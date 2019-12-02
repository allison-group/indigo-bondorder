import argparse
from collections import defaultdict
from pathlib import Path
import sys
from time import perf_counter

from indigox.config import DEFAULT_METHOD
from indigox.exception import IndigoMissingFunctionality, IndigoUnfeasibleComputation
from indigox.formalbonds import FormalBondOrders
from indigox.periodictable import PeriodicTable as PT

import openbabel as ob


mol_fcs = defaultdict(int)

#with Path("/Users/iwelsh/Dropbox/PhD/BondOrderValidationMolecules/MMFF94/MMFF94.fc_hypervalent").open('r') as f:
for line in Path("/Users/iwelsh/Dropbox/PhD/BondOrderValidationMolecules/MMFF94/MMFF94.fc_hypervalent").open('r'):
    if line.startswith("Molecule"):
        current = line.split()[1]
    elif line.startswith("Atom"):
        mol_fcs[current] += int(line.split()[-1])

def print_optimised(G):
    methods = {'fpt':'Fixed Parameter Optimisation',
               'dynamic':'Fixed Parameter Optimisation',
               'lo':'Local Minimisation',
               'local':'Local Minimisation',
               'a*':'A* search',
               'astar':'A* search',
               'ga':'Genetic Optimisation',
               'genetic':'Genetic Optimisation',
               'ball':'BALL A* search'}
    print_str = """Optimisation of {name} using {method} method.{atoms}{bonds}
Optimised energy was {energy:.5f} hartrees.
""" 
    formatting = {'name':G.graph['name'],
                  'method':methods[G.graph['bond_order_method']],
                  'atoms':[],
                  'bonds':[],
                  'energy':G.graph['bond_order_energy']}
    
    for atom in G:
        if G.node[atom]['formal_charge']:
            G.node[atom]['atom'] = atom
            if not len(formatting['atoms']):
                formatting['atoms'].append('\nAtoms with non-zero formal charges:')
            formatting['atoms'].append('{atom:^3} {element:<2} {formal_charge:>2}'
                      ''.format(**G.node[atom]))
    for a, b in G.edges():
        if G[a][b]['order'] > 1:
            if not len(formatting['bonds']):
                formatting['bonds'].append('\nBonds with order greater than 1:')
            formatting['bonds'].append('{:^3} {:^3} {:>1}'
                                       ''.format(a,b,G[a][b]['order']))
    formatting['atoms'] = '\n'.join(formatting['atoms'])
    formatting['bonds'] = '\n'.join(formatting['bonds'])
    print(print_str.format(**formatting))

def obMolToCXXBuildFunction(obMol):
    all_elements = set()
    atom_ads = []
    bond_ads = []
    atom_idx_to_pos = dict()
    formatting = dict(name=obMol.GetTitle())
    formatting['q'] = mol_fcs[formatting['name']]
    
    build_string = """void Build{name}(std::shared_ptr<indigo-bondorder::graph::MolecularGraph> G) {{
  using namespace indigo_bondorder::graph;
  using namespace indigo_bondorder;
  
  IXPeriodicTable pt = IXPeriodicTable();
  G->Clear();
  MolVertProp {element_list};
  std::vector<MolVertex> a;
  G->SetTotalCharge({q});
  
  {setting_elements}
  {adding_vertices}
  {adding_edges}
}}
"""
    for atom in ob.OBMolAtomIter(obMol):
        a = atom.GetIdx()
        element = PT[atom.GetAtomicNum()].symbol
        atom_idx_to_pos[a] = len(atom_ads)
        atom_ads.append("a.push_back(G->AddVertex({}));".format(element))
        all_elements.add(element)
    for bond in ob.OBMolBondIter(obMol):
        a = bond.GetBeginAtomIdx();
        b = bond.GetEndAtomIdx();
        bond_ads.append("G->AddEdge(a[{}], a[{}]);".format(atom_idx_to_pos[a],
                                                           atom_idx_to_pos[b]))
    
    formatting["element_list"] = ", ".join(sorted(all_elements))
    setting_elements = []
    for e in sorted(all_elements):
        setting_elements.append("{0}.element = std::make_shared<data::IXElement>(*pt.GetElement(\"{0}\"));".format(e))
    formatting['setting_elements'] = "\n  ".join(setting_elements)
    formatting['adding_vertices'] = "\n  ".join(atom_ads)
    formatting['adding_edges'] = "\n  ".join(bond_ads)
    
    return build_string.format(**formatting), "void Build{name}(std::shared_ptr<indigo-bondorder::graph::MolecularGraph> G);".format(**formatting)

def main():
    cl_parser = argparse.ArgumentParser(
        description='Determine bond orders and formal charges of molecules.')
    cl_parser.add_argument('file', type=Path,
                           help='molecule structure file to optimise')
    cl_parser.add_argument('-i', '--input', nargs='?', default='pdb', type=str,
                           help='format of the input files')
    cl_parser.add_argument('-m', '--method', nargs='?', default=DEFAULT_METHOD,
                           help='optimisation method to use')
    cl_parser.add_argument('-q', '--charge', nargs='?', default=0, type=int,
                           help='total molecular charge of the molecule')
    cl_parser.add_argument('--no-printout', action='store_true')
    
    args = cl_parser.parse_args()
    
    convert = ob.OBConversion()
    convert.SetInAndOutFormats(args.input, 'pdb')
    obMol = ob.OBMol()
    if not args.file.is_file():
        print('Unable to read from file: {}'.format(args.file), file=sys.stderr)
        sys.exit(1)
    not_at_end = convert.ReadFile(obMol, str(args.file))
    start_time = perf_counter()
    count = 1
    header = [];
    implementation = [];
    seenbodkou = False
    while not_at_end:
        count += 1
        if obMol.GetData("COMPND").GetValue().strip() != "FASJIS":
            not_at_end = convert.Read(obMol)
            continue
        obMol.SetTitle(obMol.GetData("COMPND").GetValue().strip())
        #print(obMol.GetData("COMPND").GetValue().strip())
        try:
            opt = FormalBondOrders.determine_bond_orders(obMol, 'a*', 
                                                         mol_fcs[obMol.GetData("COMPND").GetValue().strip()])
            opt.graph['name'] = obMol.GetData("COMPND").GetValue().strip()
        except IndigoUnfeasibleComputation:
            not_at_end = convert.Read(obMol)
            continue
        
        if obMol.GetTitle() == "BODKOU":
            if seenbodkou:
                obMol.SetTitle("BODKOU_")
            else:
                seenbodkou = True;
        
        i, h = obMolToCXXBuildFunction(obMol)
        header.append(h)
        implementation.append(i)
#         main_str = """Build{0}(G);
#   fpt.SetMolecularGraph(G);
#   astar.SetMolecularGraph(G);
#   ElectronOpt::Configuration::ALGORITHM = ElectronOpt::Algorithm::FPT;
#   fpt_count = fpt.Run();
#   fpt_ene = fpt.GetMinimisedEnergy();
#   ElectronOpt::Configuration::ALGORITHM = ElectronOpt::Algorithm::ASTAR;
#   try {{
#     astar_count = astar.Run();
#     astar_ene = astar.GetMinimisedEnergy();
#   }} catch (std::exception& e) {{
#     astar_count = 0;
#     astar_ene = 0;
#   }}
#   if (fpt_ene != astar_ene || fpt_count != astar_count)
#     std::cout << "{0} A* does not match FPT. " << astar_count << ", " << astar_ene << " vs. " << fpt_count << ", " << fpt_ene << std::endl;
#   else std::cout << "{0} A* matches FPT. " << astar_count << ", " << astar_ene << std::endl;
#         """
#         print(main_str.format(obMol.GetTitle()))
        
        if not args.no_printout:
            print_optimised(opt)
        not_at_end = convert.Read(obMol)

    print("{:.3f} seconds run time".format(perf_counter()-start_time))
#     with Path("/Users/iwelsh/GitHub/indigox/indigox/mmff94.hpp").open('w') as f:
#         h_str = """#include "maths/molecular_graph.hpp"
#  
# {}
#         """.format("\n".join(header))
#         f.write(h_str)
#          
#     with Path("/Users/iwelsh/GitHub/indigox/indigox/mmff94.cpp").open('w') as f:
#         i_str = """#include <vector>
#  
# #include "classes/periodictable.hpp"
# #include "mmff94.hpp"
#  
# {}
#         """.format("\n\n".join(implementation))
#         f.write(i_str)
     
        
            

if __name__ == "__main__":
    main()
