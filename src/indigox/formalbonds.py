import warnings

from indigox.astar import AStar
from indigox.config import SUPPORTED_ELEMENTS, BALL_AVAILABLE, DEFAULT_METHOD
from indigox.exception import (IndigoUnfeasibleComputation,
                               IndigoMissingParameters, IndigoWarning)
from indigox.fpt import FPT
from indigox.ga import GeneticAlogrithm
from indigox.lopt import LocalOptimisation
from indigox.misc import obmol_to_graph, BondOrderAssignment
from indigox.periodictable import PeriodicTable as PT
from networkx import is_connected

import openbabel as ob


if BALL_AVAILABLE:
    from indigox.ball import BallOpt



class FormalBondOrders(BondOrderAssignment):
    
    def __init__(self, obMol, method, total_charge):
        self.sys = obMol
        self.G = obmol_to_graph(obMol, total_charge)
        self.method = method.lower()
        all_e = set(self.G.node[n]['element'] for n in self.G)
        if not is_connected(self.G):
            raise IndigoUnfeasibleComputation('Cannot calculate bond orders and '
                                              'formal charges on disconnected '
                                              'molecules.')
        if not all_e.issubset(SUPPORTED_ELEMENTS):
            raise IndigoUnfeasibleComputation('Cannot calculate bond orders and '
                                              'formal charges with {} elements'
                                        ''.format(all_e - SUPPORTED_ELEMENTS))
        if DEFAULT_METHOD not in ['fpt','a*']:
            warnings.warn("Default optimisation method is not exact or does "
                          "not optimsed formal charges and bond orders.",
                          IndigoWarning)
        
            
                          
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
    
    def initialise(self):
        pass
    
    def run(self):
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
    def determine_bond_orders(cls, system, method=DEFAULT_METHOD, total_charge=0):
        return cls(system, method, total_charge).run()
        
    