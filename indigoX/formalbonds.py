from indigo.algorithm.bondorder.dynamic import DynamicFPT
from indigo.algorithm.bondorder.localoptimisation import LocalOptimisation
from indigo.config import SUPPORTED_ELEMENTS
from indigo.data.basedata import BaseAlgorithm
from indigo.data.periodictable import PeriodicTable as PT
from indigo.exception import IndigoUnfeasibleComputation, IndigoMissingParameters
from indigo.graph.modgraph import nx_graph
import openbabel as ob


class FormalBondOrders(BaseAlgorithm):
    def __init__(self, system, method='fpt'):
        self.sys = system
        self.G = nx_graph(system)
        self.method = method.lower()
        all_e = set(self.G.node[n]['element'] for n in self.G)
        if not all_e.issubset(SUPPORTED_ELEMENTS):
            raise IndigoUnfeasibleComputation('Cannot calculate bond orders and '
                                              'formal charges with {} elements'
                                        ''.format(all_e - SUPPORTED_ELEMENTS))
                          
    def aromatics(self, assignment):
        mol = ob.OBMol()
        for atom in assignment:
            obAtom = mol.NewAtom(atom)
            obAtom.SetAtomicNum(PT[assignment.node[atom]['element']].number)
            obAtom.SetFormalCharge(assignment.node[atom]['formal_charge'])
            
        for a, b in assignment.edges():
            mol.AddBond(a, b, assignment[a][b]['order'])
            
        for bond in ob.OBMolBondIter(mol):
            a = bond.GetBeginAtomIdx()
            b = bond.GetEndAtomIdx()
            assignment[a][b]['aromatic'] = bond.IsAromatic()
    
    def optimise(self):
        if self.method == 'fpt':
            processor = DynamicFPT(self.G)
        elif self.method == 'lo':
            processor = LocalOptimisation(self.G)
        else:
            raise IndigoMissingParameters('Unknown optimisation method: {}'
                                          ''.format(self.method))
        assignment, ene = processor.run()
        self.aromatics(assignment)
        assignment.graph['bond_order_energy'] = ene
        assignment.graph['bond_order_method'] = self.method
        return assignment
            
    
    @classmethod
    def determine_bond_orders(cls, system, method='fpt'):
        return cls(system, method).optimise()
        
    