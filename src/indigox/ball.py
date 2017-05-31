from copy import deepcopy
from indigox.config import BALL_DATA_FILE, INFINITY, MAX_SOLUTIONS
from indigox.misc import BondOrderAssignment, graph_to_dist_graph, node_energy
try:
    import BALLCore as BALL
    BALL_AVAILABLE = True
    BALL_ELEMENTS = dict(
        H=BALL.PTE['H'],  He=BALL.PTE['HE'], Li=BALL.PTE['LI'],
        Be=BALL.PTE['BE'], B=BALL.PTE['B'],  C=BALL.PTE['C'],
        N=BALL.PTE['N'],  O=BALL.PTE['O'],  F=BALL.PTE['F'],
        Ne=BALL.PTE['NE'], Na=BALL.PTE['NA'], Mg=BALL.PTE['MG'],
        Al=BALL.PTE['AL'], Si=BALL.PTE['SI'], P=BALL.PTE['P'],
        S=BALL.PTE['S'],  Cl=BALL.PTE['CL'], Ar=BALL.PTE['AR'],
        K=BALL.PTE['K'],  Ca=BALL.PTE['CA'], Sc=BALL.PTE['SC'],
        Ti=BALL.PTE['TI'], V=BALL.PTE['V'],  Cr=BALL.PTE['CR'],
        Mn=BALL.PTE['MN'], Fe=BALL.PTE['FE'], Co=BALL.PTE['CO'],
        Ni=BALL.PTE['NI'], Cu=BALL.PTE['CU'], Zn=BALL.PTE['ZN'],
        Ga=BALL.PTE['GA'], Ge=BALL.PTE['GE'], As=BALL.PTE['AS'],
        Se=BALL.PTE['SE'], Br=BALL.PTE['BR'], Kr=BALL.PTE['KR'],
        Rb=BALL.PTE['RB'], Sr=BALL.PTE['SR'], Y=BALL.PTE['Y'],
        Zr=BALL.PTE['ZR'], Nb=BALL.PTE['NB'], Mo=BALL.PTE['MO'],
        Tc=BALL.PTE['TC'], Ru=BALL.PTE['RU'], Rh=BALL.PTE['RH'],
        Pd=BALL.PTE['PD'], Ag=BALL.PTE['AG'], Cd=BALL.PTE['CD'],
        In=BALL.PTE['IN'], Sn=BALL.PTE['SN'], Sb=BALL.PTE['SB'],
        Te=BALL.PTE['TE'], I=BALL.PTE['I'],  Xe=BALL.PTE['XE'],
        Cs=BALL.PTE['CS'], Ba=BALL.PTE['BA'], La=BALL.PTE['LA'],
        Ce=BALL.PTE['CE'], Pr=BALL.PTE['PR'], Nd=BALL.PTE['ND'],
        Pm=BALL.PTE['PM'], Sm=BALL.PTE['SM'], Eu=BALL.PTE['EU'],
        Gd=BALL.PTE['GD'], Tb=BALL.PTE['TB'], Dy=BALL.PTE['DY'],
        Ho=BALL.PTE['HO'], Er=BALL.PTE['ER'], Tm=BALL.PTE['TM'],
        Yb=BALL.PTE['YB'], Lu=BALL.PTE['LU'], Hf=BALL.PTE['HF'],
        Ta=BALL.PTE['TA'], W=BALL.PTE['W'],  Re=BALL.PTE['RE'],
        Os=BALL.PTE['OS'], Ir=BALL.PTE['IR'], Pt=BALL.PTE['PT'],
        Au=BALL.PTE['AU'], Hg=BALL.PTE['HG'], Tl=BALL.PTE['TL'],
        Pb=BALL.PTE['PB'], Bi=BALL.PTE['BI'], At=BALL.PTE['AT'],
        Rn=BALL.PTE['RN'], Fr=BALL.PTE['FR'], Ra=BALL.PTE['RA'],
        Ac=BALL.PTE['AC'], Th=BALL.PTE['TH'], Pa=BALL.PTE['PA'],
        U=BALL.PTE['U'],  Np=BALL.PTE['NP'], Pu=BALL.PTE['PU'],
        Po=BALL.PTE['PO'], Am=BALL.PTE['AM'], Cm=BALL.PTE['CM'],
        Bk=BALL.PTE['BK'], Cf=BALL.PTE['CF'], Es=BALL.PTE['ES'],
        Fm=BALL.PTE['FM'], Md=BALL.PTE['MD'], No=BALL.PTE['NO'],
        Lr=BALL.PTE['LR'], Rf=BALL.PTE['RF'], Db=BALL.PTE['DB'],
        Sg=BALL.PTE['SG'], Bh=BALL.PTE['BH'], Hs=BALL.PTE['HS'],
        Mt=BALL.PTE['MT'],)
    
    # setup the bond order processor
    bop = BALL.AssignBondOrderProcessor()
    # alias' for long name
    opts = BALL.AssignBondOrderProcessor.Option
    algo = BALL.AssignBondOrderProcessor.Algorithm
    
    bop.options.setBool(opts.KEKULIZE_RINGS, True)
    bop.options.setBool(opts.OVERWRITE_SINGLE_BOND_ORDERS, True)
    bop.options.setBool(opts.OVERWRITE_DOUBLE_BOND_ORDERS, True)
    bop.options.setBool(opts.OVERWRITE_TRIPLE_BOND_ORDERS, True)
    bop.options.set(opts.ALGORITHM, algo.A_STAR)
    bop.options.setReal(opts.BOND_LENGTH_WEIGHTING, 0)
    bop.options.setInteger(opts.MAX_NUMBER_OF_SOLUTIONS, MAX_SOLUTIONS)
    bop.options.setBool(opts.COMPUTE_ALSO_NON_OPTIMAL_SOLUTIONS, False)
    bop.options.setBool(opts.ADD_HYDROGENS, False)
    bop.options.set(opts.INIFile, str(BALL_DATA_FILE))
except ImportError:
    BALL_AVAILABLE = False

class BallOpt(BondOrderAssignment):
    def __init__(self, G):
        self.init_G = G
        
    def initialise(self):
        self.G = graph_to_dist_graph(self.init_G)
        self.system = BALL.System()
        self.mol = BALL.Molecule()
        self.atoms = {}
        self.bonds = []
        
        for a, d in self.init_G.nodes(True):
            ball_e = BALL_ELEMENTS[d['element']]
            atom = BALL.Atom()
            atom.setName(str(a))
            atom.setElement(ball_e)
            self.atoms[a] = atom
            
        for a, b, d in self.init_G.edges(data=True):
            bond = self.atoms[a].createBond(self.atoms[b])
            bond.setOrder(1)
            self.bonds.append(bond)
        
        for atom in self.atoms.values():
            self.mol.insert(atom)
            
        self.system.insert(self.mol)

    def run(self):
        if not BALL_AVAILABLE:
            self.log.warning('BALL method is unavailable as BALLCore could not '
                             'be loaded.')
            for x in self.init_G:
                self.init_G.node[x]['formal_charge'] = 0
                for y in self.init_G[x]:
                    self.init_G[x][y]['order'] = 1
            return self.init_G, INFINITY
        else:
            self.log.warning("BALL method selected. Formal charges will not be "
                             "optimised.")
        best_ene = INFINITY * INFINITY
        best_g = None
        
        self.initialise()
        self.system.apply(bop)
        for i in range(bop.getNumberOfComputedSolutions()):
            bop.apply(i)
            for atom in BALL.atoms(self.system):
                a = int(str(atom.getName()))
                fc = int(atom.getFormalCharge())
                self.G.node[(a,)]['fc'] = fc
                
            for bond in BALL.bonds(self.system):
                a = int(str(bond.getFirstAtom().getName()))
                b = int(str(bond.getSecondAtom().getName()))
                bo = int(bond.getOrder())
                if a > b:
                    a, b = b, a
                self.G.node[(a, b)]['e-'] = bo * 2
            
            i_ene = round(sum(node_energy(self.G, n) for n in self.G),5)
            if i_ene < best_ene:
                best_ene = i_ene
                best_g = self.assignment_to_graph()
                
        return best_g, best_ene
    
    def assignment_to_graph(self):
        G = deepcopy(self.init_G)
        for v in self.G:
            if len(v) == 1:
                G.node[v[0]]['formal_charge'] = 0
            if len(v) == 2:
                G[v[0]][v[1]]['order'] = self.G.node[v]['e-'] // 2
        return G
        
       
        
