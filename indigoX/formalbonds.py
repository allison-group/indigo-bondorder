from indigoX.astar import AStar
from indigoX.config import SUPPORTED_ELEMENTS, BALL_AVAILABLE
from indigoX.exception import (IndigoUnfeasibleComputation,
                               IndigoMissingParameters)
from indigoX.fpt import FPT
from indigoX.ga import GeneticAlogrithm
from indigoX.lopt import LocalOptimisation
from indigoX.misc import obmol_to_graph, _get_logger
from indigoX.periodictable import PeriodicTable as PT
import openbabel as ob
if BALL_AVAILABLE:
    from indigoX.ball import BallOpt



class FormalBondOrders(object):
    log = _get_logger()
    
    def __init__(self, obMol, method='fpt', total_charge=0):
        self.sys = obMol
        self.G = obmol_to_graph(obMol)
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
        if self.method in ['fpt','dynamic']:
            processor = FPT(self.G)
        elif self.method in ['lo','local']:
            processor = LocalOptimisation(self.G)
        elif self.method in ['a*','astar']:
            processor = AStar(self.G)
        elif self.method in ['ga','genetic']:
            processor = GeneticAlogrithm(self.G)
        elif BALL_AVAILABLE and self.method == 'ball':
            processor = BallOpt(self.G)
            self.log.warning("BALL method selected. Formal charges will not be "
                             "optimised.")
        elif not BALL_AVAILABLE and self.method == 'ball':
            raise IndigoUnfeasibleComputation('The BALL library cannot be loaded.')
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
        
    