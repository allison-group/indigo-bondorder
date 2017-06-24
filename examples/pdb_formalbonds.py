import argparse
from pathlib import Path
import sys

from indigox.config import DEFAULT_METHOD
from indigox.exception import IndigoMissingFunctionality
from indigox.formalbonds import FormalBondOrders
import openbabel as ob

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
    while not_at_end:
        opt = FormalBondOrders.determine_bond_orders(obMol, args.method, 
                                                 args.charge)
    
        if not args.no_printout:
            print_optimised(opt)
        not_at_end = convert.Read(obMol)
        
            

if __name__ == "__main__":
    main()
