import BALLCore as BALL
from indigo.config import SOURCE_DIR, INFINITY, MAX_SOLUTIONS
from indigo.misc import BondOrderAssignment, graph_to_dist_graph, node_energy
                   
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
bop.options.set(opts.INIFile, str(SOURCE_DIR / 'data' / 'OriginalBO.xml'))

class BallOpt(BondOrderAssignment):
    def __init__(self, G, include_fc=False, ref_ene=None, bond_ene=None):
        self.init_G = G
        self.fc = include_fc
        self.ref = ref_ene
        self.bond_ene = bond_ene
        
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
            if self.fc:
                atom.setFormalCharge(d['fc'])
            else:
                atom.setFormalCharge(0)
            self.atoms[a] = atom
            
        for a, b, d in self.init_G.edges(data=True):
            bond = self.atoms[a].createBond(self.atoms[b])
            bond.setOrder(1)
            self.bonds.append(bond)
        
        for atom in self.atoms.values():
            self.mol.insert(atom)
            
        self.system.insert(self.mol)

    def run(self):
        best_ref = INFINITY * INFINITY
        best_bond = INFINITY * INFINITY
        matched_ref = matched_bond = False
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
            b_only = round(sum(node_energy(self.G, n) for n in self.G if len(n) == 2),5)
            if not matched_ref and -1e-10 < self.ref - i_ene < 1e-10:
                matched_ref = True
                best_ref = i_ene
            if not matched_bond and -1e-10 < self.bond_ene - b_only < 1e-10:
                matched_bond = True
                best_bond = b_only
            if not matched_ref and i_ene < best_ref:
                best_ref = i_ene
            if not matched_bond and b_only < best_bond:
                best_bond = b_only
                
            if matched_bond and matched_ref:
                break
                
        return None, (best_ref, best_bond)
       
        
