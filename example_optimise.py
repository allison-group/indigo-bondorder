from indigox.exception import IndigoMissingFunctionality
from indigox.formalbonds import FormalBondOrders

import openbabel as ob


def build_mol(molecule='benzene'):
    if molecule == 'benzene':
        molecule = {'atoms':[(1,6),(2,6),(3,6),(4,6),(5,6),(6,6),
                             (7,1),(8,1),(9,1),(10,1),(11,1),(12,1),],
                    'bonds':[(1,2),(1,6),(1,7),(2,3),(2,8),(3,4),(3,9),
                             (4,5),(4,10),(5,6),(5,11),(6,12)],}
        charge = 0
    elif molecule == 'diacid-diol':
        molecule = {'atoms':[(1,8),(2,8),(3,8),(4,8),(5,8),(6,8),
                             (7,6),(8,6),(9,6),(10,6),(11,1),(12,1),
                             (13,1),(14,1),(15,1),],
                    'bonds':[(1,7),(1,15),(2,7),(3,8),(3,13),(4,9),(4,14),
                             (5,10),(6,10),(7,8),(8,9),(8,11),(9,10),(9,12)]}
        charge = -1
    else:
        raise IndigoMissingFunctionality('Unknown molecule: {}'.format(molecule))
    
    
    mol = ob.OBMol()
    for idx, z in molecule['atoms']:
        obAtom = mol.NewAtom(idx)
        obAtom.SetAtomicNum(z)
        
    for a, b in molecule['bonds']:
        mol.AddBond(a,b, 1)
        
    return mol, charge

def print_optimised(G):
    print('Optimisation method used: {}'.format(G.graph['bond_order_method'].upper()))
    print('Optimised energy: {:.5f}'.format(G.graph['bond_order_energy']))
    if G.graph['bond_order_method'] != 'ball':
        print('Optimised formal charges:')
        for atom in G:
            if G.node[atom]['formal_charge']:
                G.node[atom]['atom'] = atom
                print('{atom:>2} {element:<2} {formal_charge:>2}'
                      ''.format(**G.node[atom]))
    print('Optimised bond orders:')
    for a, b in G.edges():
        if G[a][b]['order'] > 1:
            print('{:>2} {:>2} {:>1}'.format(a,b,G[a][b]['order']))
        

def main():
    for mol in ['benzene','diacid-diol']:
        test_mol, charge = build_mol(mol)
        print('Optimisting {}'.format(mol))
        for method in ['fpt','lo','a*','ga','ball']:
            opt_bofc = FormalBondOrders.determine_bond_orders(test_mol, method, charge)
            print_optimised(opt_bofc)
            print('\n')
    
if __name__ == "__main__":
    main()