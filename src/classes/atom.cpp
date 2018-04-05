/** @file atom.cpp
 *  @brief Atom implementation
 *  @author Ivan Welsh
 *  @date 5 January 2018
 *  @lastmodify 8 January 2018
 *  @version 0.1
 *  @copyright The MIT License
 */

#include <algorithm>
#include <cstdint>
#include <memory>
#include <sstream>

#include "api.hpp"
#include "classes/atom.hpp"
#include "classes/periodictable.hpp"
#include "utils/counter.hpp"

namespace indigox {
  
  /*
   *  Atom implementation
   */
  
  // Initalisation methods
  Atom::Atom() : utils::CountableObject<Atom>(), mol_(), element_() {
    idx_ = GetUniqueID();
  }
  
  /** @param m the molecule to assign this atom to. */
  Atom::Atom(Molecule_p m)
  : utils::CountableObject<Atom>(), mol_(m), element_() {
    idx_ = GetUniqueID();
  }
  
  Atom::~Atom() { }

  // Data retrival methods
  Element_p Atom::GetElement() const { return element_; }
  Int Atom::GetFormalCharge() const { return (Int)fc_; }
  uid_t Atom::GetIndex() const { return idx_; }
  Molecule_p Atom::GetMolecule() const {
    if (!mol_.expired()) return mol_.lock();
    else return Molecule_p();
  }
  String Atom::GetName() const { return name_; }
  Float Atom::GetX() const { return pos_[0]; }
  Float Atom::GetY() const { return pos_[1]; }
  Float Atom::GetZ() const { return pos_[2]; }
  String Atom::ToString() {
    std::stringstream ss;
    ss << idx_ << " \"" << name_ << "\" @ (" << pos_[0] << "," << pos_[1] ;
    ss << "," << pos_[2] << ")";
    return ss.str();
  }
  
  // Data modification methods
  /** @param e the element to set. */
  void Atom::SetElement(Element_p e) { element_ = e; }
  
  /** @param e the name or symbol of the element to set. */
  void Atom::SetElement(String e) {
    PeriodicTable_p PT = PeriodicTable::GetInstance();
    element_ = PT->GetElement(e);
  }
  
  /** @param e the atomic number of the element to set. */
  void Atom::SetElement(Int e) {
    PeriodicTable_p PT = PeriodicTable::GetInstance();
    element_ = PT->GetElement((uint8_t)e);
  }
  
  /** @param q the value to give this atom's formal charge. */
  void Atom::SetFormalCharge(Int q) { fc_ = q; }
  
  /** @details Index is used for internal bookkeeping and should not be
   *  considered stable.
   *  @param i the index to set.
   */
  void Atom::SetIndex(uid_t i) { idx_ = i; }
  
  /** @param m the molecule to set. */
  void Atom::SetMolecule(Molecule_p m) { mol_ = m; }
  
  /** @param n the name to set. */
  void Atom::SetName(String n) { name_ = n; }
  
  /** @param x the x coordinate to set.
   *  @param y the y coordinate to set.
   *  @param z the z coordinate to set.
   */
  void Atom::SetPosition(Float x, Float y, Float z) {
    pos_[0] = x;
    pos_[1] = y;
    pos_[2] = z;
  }
  
  // Atom modification methods
  void Atom::AddBond(Bond_p b) {
    bonds_.emplace_back(b);
  }
  
  void Atom::RemoveBond(Bond_p b) {
    AtomBondIterator it = bonds_.begin();
    for (; it != bonds_.end(); ++it) {
      if (it->lock() == b) break;
    }
    if (it != bonds_.end()) bonds_.erase(it);
  }
  
  void Atom::Clear() {
    mol_.reset();
    element_.reset();
    bonds_.clear();
    pos_.fill(0.0);
  }
  
  // Iterator methods
  Bond_p Atom::Begin(AtomBondIterator &it) {
    for (it = bonds_.begin(); it != bonds_.end(); ++it) {
      if (!it->expired()) break;
    }
    return (it == bonds_.end()) ? Bond_p() : it->lock();
  }
  
  Bond_p Atom::Next(AtomBondIterator &it) {
    for (++it; it != bonds_.end(); ++it) {
      if (!it->expired()) break;
    }
    return (it == bonds_.end()) ? Bond_p() : it->lock();
  }
  
  AtomBondIterator Atom::BeginBond() { return bonds_.begin(); }
  AtomBondIterator Atom::EndBond() { return bonds_.end(); }
  
}
