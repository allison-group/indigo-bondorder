/** @file atom.hpp
 *  @brief Atom declaration
 *  @author Ivan Welsh
 *  @date 5 January 2018
 *  @lastmodify 8 January 2018
 *  @version 0.1
 *  @copyright The MIT License
 */

#ifndef INDIGOX_CLASSES_ATOM_HPP
#define INDIGOX_CLASSES_ATOM_HPP

#include <array>
#include <cstdint>
#include <iostream>
#include <memory>
#include <vector>

#include "../api.hpp"
#include "../utils/counter.hpp"

namespace indigox {

  // Related typedefs
  typedef std::vector<Bond_wp> AtomBonds;
  typedef AtomBonds::iterator AtomBondIterator;
  typedef AtomBonds::const_iterator const_AtomBondIterator;
  
  class Atom
  : public utils::CountableObject<Atom>,
  public std::enable_shared_from_this<Atom> {
  
  public:
    /// @name Initialisation methods.
    
    /// @brief Default constructor.
    Atom();
    
    /// @brief Normal constructor, links Atom to a Molecule.
    Atom(Molecule_p m);
    
    /// @brief Destructor
    ~Atom();
    
    /// @name Data retrival methods.
    
    /// @returns the element of this atom.
    Element_p GetElement() const;
    
    /// @returns the formal charge on this atom.
    Int GetFormalCharge() const;
    
    /// @returns the index of this atom.
    uid_t GetIndex() const;
    
    /// @returns the molecule this atom is part of.
    Molecule_p GetMolecule() const;
    
    /// @returns the name of this atom.
    String GetName() const;
    
    /// @returns the x coordinate of this atom.
    Float GetX() const;
    
    /// @returns the y coordinate of this atom.
    Float GetY() const;
    
    /// @returns the z coordinate of this atom.
    Float GetZ() const;
    
    /// @returns a string representation of this atom.
    String ToString();
    
    /// @name Data modification methods.
    
    /// @brief Sets the element of this atom.
    void SetElement(Element_p);
    void SetElement(String);
    void SetElement(Int);
    
    /// @brief Sets the formal charge on this atom.
    void SetFormalCharge(Int);
    
    /// @brief Sets the index of this atom.
    void SetIndex(uid_t);
    
    /// @brief Sets the molecule this atom is part of. Performs no bookkeeping.
    void SetMolecule(Molecule_p);
    
    /// @brief Sets the name of this atom.
    void SetName(String);
    
    /// @brief Sets the position of this atom
    void SetPosition(Float, Float, Float);
    
    /// @name Atom modification methods.
    
    /// @brief Adds a bond to this atom. Performs no bookkeeping.
    void AddBond(Bond_p);
    
    /// @brief Removes a bond from this atom. Performs no bookkeeping.
    void RemoveBond(Bond_p);
    
    void Clear();
    
    /// @name Iterator methods.
    
    /// @brief Sets the given iterator to the next position in bonds_.
    Bond_p Next(AtomBondIterator&);
    
    /// @brief Sets the given iterator to the start of bonds_.
    Bond_p Begin(AtomBondIterator&);
    
    /// @brief Sets the given iterator to the end of bonds_.
    AtomBondIterator BeginBond();
    AtomBondIterator EndBond();
    
  private:
    Molecule_wp mol_;
    Element_p element_;
    Int fc_ = 0;
    uid_t idx_ = 0;
    String name_ = "None";
    std::array<Float, 3> pos_ = {{0.0, 0.0, 0.0}};
    
    AtomBonds bonds_;
    
  };
}

#endif /* INDIGOX_CLASSES_ATOM_HPP */
