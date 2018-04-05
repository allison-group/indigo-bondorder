import indigox as ix

PT = ix.PeriodicTable()
H = PT["H"]
C = PT["C"]
O = PT["O"]
N = PT["N"]

# Build a molecule
mol = ix.Molecule()
mol.SetTotalCharge(-1)
# Add all the atoms
c1 = mol.NewAtom(C) ; c1.SetName("C1")
c2 = mol.NewAtom(C) ; c2.SetName("C2")
c3 = mol.NewAtom(C) ; c3.SetName("C3")
c4 = mol.NewAtom(C) ; c4.SetName("C4")
c5 = mol.NewAtom(C) ; c5.SetName("C5")
c6 = mol.NewAtom(C) ; c6.SetName("C6")
h7 = mol.NewAtom(H) ; h7.SetName("H7")
h8 = mol.NewAtom(H) ; h8.SetName("H8")
n9 = mol.NewAtom(N) ; n9.SetName("N9")
h10 = mol.NewAtom(H) ; h10.SetName("H10")
h11 = mol.NewAtom(H) ; h11.SetName("H11")
c12 = mol.NewAtom(C) ; c12.SetName("C12")
o13 = mol.NewAtom(O) ; o13.SetName("O13")
o14 = mol.NewAtom(O) ; o14.SetName("O14")
o15 = mol.NewAtom(O) ; o15.SetName("O15")
o16 = mol.NewAtom(O) ; o16.SetName("O16")
# Add all the bonds
for a, b in [(c1,c2),(c1,c6),(c1,h7),(c2,c3),(c2,h8),(c3,c4),(c3,n9),
             (c4,c5),(c4,h10),(c5,c6),(c5,h11),(c6,c12),(n9,o13),(n9,o14),
             (c12,o15),(c12,o16)]:
    mol.NewBond(a,b)

# Print out some information on the molecule
print("Number of atoms: {}, number of bonds: {}".format(mol.NumAtoms(), mol.NumBonds()))

# Setup to use the FPT algorithm with single electrons without preplacing
# to calculate bond orders and formal charges
opts = ix.Options.AssignElectrons
opts.ALGORITHM = opts.Algorithm.FPT
opts.FPT.ADD_EDGES_TO_TD = False
opts.FPT.MINIMUM_PROPAGATION_DEPTH = 1
opts.USE_ELECTRON_PAIRS = False

# Calculate bond orders and formal charges
count = mol.AssignElectrons()
print("{} resonace structure(s) calculated with a score of {}.".format(count, mol.GetMinimumElectronAssignmentScore()))

# Print out each of the structures
for i in range(count):
    mol.ApplyElectronAssignment(i)
    print("Structure {}:".format(i))
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0:
            print("Atom {} has a formal charge of {}.".format(atom.GetName(), atom.GetFormalCharge()))
    for bond in mol.GetBonds():
        if bond.GetOrder() != ix.BondOrder.SINGLE_BOND:
            print("Bond between {} and {} is a {}".format(bond.GetSourceAtom().GetName(), bond.GetTargetAtom().GetName(), bond.GetOrder()))
