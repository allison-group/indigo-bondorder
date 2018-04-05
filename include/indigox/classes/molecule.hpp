/** @file molecule.hpp
 *  @brief Declaration of Molecule class
 *  @author Ivan Welsh
 *  @date 6 January 2018
 *  @lastmodify 8 January 2018
 *  @version 0.1
 *  @copyright The MIT License
 */
#include <map>
#include <unordered_map>
#include <vector>

#include "../api.hpp"

#include "../algorithm/electron_optimisation.hpp"

#include "molecular_graph.hpp"
#include "../utils/counter.hpp"

#ifndef INDIGOX_MOLECULE_HPP
#define INDIGOX_MOLECULE_HPP

namespace indigox {
  
  // Related typedefs
  typedef std::vector<Atom_p> MolAtoms;
  typedef std::vector<Bond_p> MolBonds;
  typedef MolAtoms::iterator MolAtomIterator;
  typedef MolBonds::iterator MolBondIterator;
  typedef MolAtoms::const_iterator const_MolAtomIterator;
  typedef MolBonds::const_iterator const_MolBondIterator;
  
  /** @class Molecule molecule.hpp classes/molecule.hpp
   *  @brief Class for storing and manipulating molecules.
   */
  class Molecule : public utils::CountableObject<Molecule>,
  public std::enable_shared_from_this<Molecule> {
    
  public:
    /// @name Initialisation methods
    
    /// @brief Default constructor
    Molecule();
    
    /// @brief Constructor for pre-named molecule
    Molecule(String name);
    
    /// @brief Destructor
    ~Molecule();
    
    /// @name Data retrival methods
    
    /// @returns a shared pointer to an atom with the given index.
    Atom_p GetAtomIndex(uid_t);
    Atom_p GetAtomUniqueID(uid_t);
    
    /// @returns a shared pointer to a bond with the given index.
    Bond_p GetBondIndex(uid_t);
    Bond_p GetBondUniqueID(uid_t);
    
    /// @returns a shared pointer to the bond between two atoms, if it exists.
    Bond_p GetBond(Atom_p, Atom_p) const;
    
    /// @returns the chemical formula of the molecule.
    String GetFormula() const;
    
    /// @returns a shared pointer to a graph representation of the molecule
    MolecularGraph_p GetMolecularGraph();
    
    /// @returns the name of the molecule.
    String GetName() const;
    
    /// @returns the molecular charge of the molecule
    Int GetTotalCharge() const;
    
    bool IsModified() const;
    
    /// @returns the number of atoms in the molecule.
    Uint NumAtoms() const;
    
    /// @returns the number of bonds in the molecule.
    Uint NumBonds() const;
    
    /// @name Data modification methods
    
    /// @brief Sets the name of the molecule
    void SetName(String);
    
    /// @brief Sets the molecular charge of the molecule.
    void SetTotalCharge(Int);
    
    /// @brief Resets the indices of contained atoms and bonds.
    void ResetIndices();
    
    /// @name Molecule modification methods
    
    /// @brief Adds a new atom to the molecule.
    /// @returns a shared pointer to the new atom.
    Atom_p NewAtom();
    Atom_p NewAtom(Element_p);
    Atom_p NewAtom(uid_t, Element_p);
    
    /// @brief Adds a new bond between two atoms to the molecule.
    /// @returns a shared pointer to the new bond.
    Bond_p NewBond(Atom_p, Atom_p);
    
    /// @brief Removes a given atom from the molecule.
    void RemoveAtom(Atom_p);
    
    /// @brief Removes a given bond from the molecule.
    void RemoveBond(Bond_p);
    
    Uint AssignElectrons();
    bool ApplyElectronAssignment(Uint);
    Score GetMinimumElectronAssignmentScore();
    
    /// @name Iterator methods
    Atom_p Begin(MolAtomIterator&);
    Atom_p Next(MolAtomIterator&);
    
    MolAtomIterator BeginAtom();
    MolAtomIterator EndAtom();
    MolBondIterator BeginBond();
    MolBondIterator EndBond();
    
    
  private:
    String name_;
    Int q_ = 0;
    MolAtoms atoms_;  // own atoms that are part of me
    MolBonds bonds_;  // own bonds that are part of me
    std::unordered_map<uid_t, Atom_wp> idx_to_atom_;
    std::unordered_map<uid_t, Bond_wp> idx_to_bond_;
    std::map<Atom_p, MolVertex> atom_to_vertex_;
    std::map<Bond_p, MolEdge> bond_to_edge_;
    MolecularGraph_p graph_;
    std::unique_ptr<ElectronOpt> elnopt_;
    bool modified_ = true;
  };
  
}

#endif /* INDIGOX_MOLECULE_HPP */
