/** @file bond.hpp
 *  @brief Bond declaration
 *  @author Ivan Welsh
 *  @date 5 January 2018
 *  @lastmodify 8 January 2018
 *  @version 0.1
 *  @copyright The MIT License
 */
#include <array>

#include "../api.hpp"
#include "../utils/counter.hpp"

#ifndef INDIGOX_CLASSES_BOND_HPP
#define INDIGOX_CLASSES_BOND_HPP

namespace indigox {
  
  // Related typedefs
  typedef std::array<Atom_wp, 2> BondAtoms;
  typedef BondAtoms::iterator BondAtomIterator;
  typedef BondAtoms::const_iterator const_BondAtomIterator;
  
  enum BondOrder {
    SINGLE_BOND = 1,
    DOUBLE_BOND,
    TRIPLE_BOND,
    QUADRUPLE_BOND,
    AROMATIC_BOND,
    ONEANDAHALF_BOND,
    TWOANDAHALF_BOND,
    UNDEFINED_BOND
  };
  
  class Bond
  : public utils::CountableObject<Bond>,
  public std::enable_shared_from_this<Bond> {
  
  public:
    
    /// @name Initialisation methods.
    
    /// @brief Default constructor.
    Bond();
    
    /// @brief Normal use constructor.
    Bond(Atom_p, Atom_p);
    
    /// @brief Destructor
    ~Bond();
    
    /// @name Data retrival methods
    
    /// @returns the index of this bond.
    uid_t GetIndex() const;
    
    /// @returns the molecule this bond is part of.
    Molecule_p GetMolecule() const;
    
    /// @returns the bond order of this bond.
    BondOrder GetOrder() const;
    
    /// @returns the source atom of this bond.
    Atom_p GetSourceAtom() const;
    
    /// @returns the target atom of this bond.
    Atom_p GetTargetAtom() const;
    
    /// @returns a string representation of this bond.
    String ToString() const;
    
    /// @name Data modification methods
    
    /// @brief Sets the index of this bond.
    void SetIndex(uid_t);
    
    /// @brief Sets the molecule this bond is part of.
    void SetMolecule(Molecule_p);
    
    /// @brief Sets the bond order of this bond.
    void SetOrder(BondOrder);
    
    /// @brief Sets the source atom of this bond.
    void SetSourceAtom(Atom_p);
    
    /// @brief Sets the target atom of this bond.
    void SetTargetAtom(Atom_p);
    
    void Clear();
    
    /// @name Iterator methods
    
    /// @brief Sets the given iterator to the next position in atoms_.
    Atom_p Next(BondAtomIterator&);
    
    /// @brief Sets the given iterator to the start of atoms_.
    Atom_p Begin(BondAtomIterator&);
    
    /// @brief Sets the given iterator to the end of atoms_.
    BondAtomIterator BeginAtom();
    BondAtomIterator EndAtom();
    
  private:
    BondAtoms atoms_;
    Molecule_wp mol_;
    BondOrder order_ = UNDEFINED_BOND;
    uid_t idx_;
    
  };
  
}

#endif /* INDIGOX_CLASSES_BOND_HPP */
