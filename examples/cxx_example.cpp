#include <iostream>
#include <indigo-bondorder/indigo-bondorder.hpp>
#include <chrono>

/**
 *
 * A simple example of Bond Order and Formal Charge assignment, using the
 * small molecule 4-Nitrobenzoate (https://pubchem.ncbi.nlm.nih.gov/compound/4419940).
 *
 * This returns 8 resonance structures because the aromatic ring, the nitro group
 * and the benzoate group can have two formal configurations each (depending on
 * which oxygens are charged and which set of bonds to represent in the aromatic
 * ring).
 *
 * Note: If this example has trouble reading file locations, you may need to set
 * Options::DATA_DIRECTORY to a the full path of your .data folder
 *
 * @return
 */
int main() {
    using namespace indigo_bondorder;

    // Set the options to use
    Options::AssignElectrons::ALGORITHM = Options::AssignElectrons::Algorithm::FPT;
    Options::AssignElectrons::FPT::ADD_EDGES_TO_TD = false;
    Options::AssignElectrons::FPT::MINIMUM_PROPAGATION_DEPTH = 1;
    Options::AssignElectrons::USE_ELECTRON_PAIRS = false;

    auto before_PT = std::chrono::high_resolution_clock::now();

    // Prepare elements
    PeriodicTable_p PT = PeriodicTable::GetInstance();
    Element_p H = PT->GetElement("H");
    Element_p C = PT->GetElement("C");
    Element_p O = PT->GetElement("O");
    Element_p N = PT->GetElement("N");

    auto after_PT = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = after_PT - before_PT;
    std::cout << "PT Generation: \t\t" << elapsed.count() << " s\n";

    // Build the molecule
    Molecule_p m = CreateMolecule();
    m->SetTotalCharge(-1);

    auto after_CreateMolecule = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_2 = after_CreateMolecule - after_PT;
    std::cout << "Create Molecule: \t" << elapsed_2.count() << " s\n";

    // Add the atoms
    Atom_p c1  = m->NewAtom(C);  c1->SetName("C1") ; 
    Atom_p c2  = m->NewAtom(C);  c2->SetName("C2") ;
    Atom_p c3  = m->NewAtom(C);  c3->SetName("C3") ;
    Atom_p c4  = m->NewAtom(C);  c4->SetName("C4") ;
    Atom_p c5  = m->NewAtom(C);  c5->SetName("C5") ;
    Atom_p c6  = m->NewAtom(C);  c6->SetName("C6") ;
    Atom_p h7  = m->NewAtom(H);  h7->SetName("H7") ;
    Atom_p h8  = m->NewAtom(H);  h8->SetName("H8") ;
    Atom_p n9  = m->NewAtom(N);  n9->SetName("N9") ;
    Atom_p h10 = m->NewAtom(H); h10->SetName("H10");
    Atom_p h11 = m->NewAtom(H); h11->SetName("H11");
    Atom_p c12 = m->NewAtom(C); c12->SetName("C12");
    Atom_p o13 = m->NewAtom(O); o13->SetName("O13");
    Atom_p o14 = m->NewAtom(O); o14->SetName("O14");
    Atom_p o15 = m->NewAtom(O); o15->SetName("O15");
    Atom_p o16 = m->NewAtom(O); o16->SetName("O16");

    auto after_AddAtoms = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_3 = after_AddAtoms - after_CreateMolecule;
    std::cout << "Add Atoms: \t\t\t" << elapsed_3.count() << " s\n";

    // Add the bonds
    m->NewBond(c1,c2);
    m->NewBond(c1,c6);
    m->NewBond(c1,h7);
    m->NewBond(c2,c3);
    m->NewBond(c2,h8);
    m->NewBond(c3,c4);
    m->NewBond(c3,n9);
    m->NewBond(c4,c5);
    m->NewBond(c4,h10);
    m->NewBond(c5,c6);
    m->NewBond(c5,h11);
    m->NewBond(c6,c12);
    m->NewBond(n9,o13);
    m->NewBond(n9,o14);
    m->NewBond(c12,o15);
    m->NewBond(c12,o16);

    auto after_AddBonds = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_4 = after_AddBonds - after_AddAtoms;
    std::cout << "Add Bonds: \t\t\t" << elapsed_4.count() << " s\n\n";
    
    // Print out some information
    std::cout << "Number of atoms: " << m->NumAtoms() << ", number of bonds: " << m->NumBonds() << "\n";

    auto before_calc = std::chrono::high_resolution_clock::now();

    // Calculate bond orders and formal charges
    Uint count = m->AssignElectrons();

    auto after_calc = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_5 = after_calc - before_calc;
    std::cout << "Calculate BO and FC: \t" << elapsed_5.count() << " s\n";

    std::cout << count << " resonance structure(s) calculated with a score of " << m->GetMinimumElectronAssignmentScore() << ".\n";

    std::vector<double> times;

    // Print out each of the structures
    for (Uint i = 0; i < count; ++i) {
        auto before_assignment = std::chrono::high_resolution_clock::now();
        m->ApplyElectronAssignment(i);
        std::chrono::duration<double> to_add = std::chrono::high_resolution_clock::now() - before_assignment;
        times.push_back(to_add.count());

        for (MolAtomIterator it = m->BeginAtom(); it != m->EndAtom(); ++it) {
            Atom_p at = *it;
            if (at->GetFormalCharge() != 0)
                std::cout << "Atom " << at->GetName() << " has a formal charge of " << at->GetFormalCharge() << ".\n";
        }
        for (MolBondIterator it = m->BeginBond(); it != m->EndBond(); ++it) {
            Bond_p bt = *it;
            if (bt->GetOrder() != 1)
                std::cout << "Bond between " << bt->GetSourceAtom()->GetName() << " and " << bt->GetTargetAtom()->GetName() << " has an order of " << bt->GetOrder() << ".\n";
        }
        std::cout << std::endl;
    }

    std::cout << "\nElectron assignment times:" << std::endl;

    for (auto time : times) {
        std::cout << "\t" << time << " s\n";
    }

    return 0;
}

