/** @file bond.cpp
 *  @brief Bond implementation
 *  @author Ivan Welsh
 *  @date 5 January 2018
 *  @lastmodify 8 January 2018
 *  @version 0.1
 *  @copyright The MIT License
 */
#include <sstream>

#include "api.hpp"
#include "classes/atom.hpp"
#include "classes/bond.hpp"
#include "utils/counter.hpp"

namespace indigox {
  
  /*
   *  Bond implementation
   */
  
  // Initalisation methods
  Bond::Bond() : utils::CountableObject<Bond>(), mol_() {
    atoms_[0] = Atom_wp();
    atoms_[1] = Atom_wp();
    idx_ = GetUniqueID();
  }
  
  /** @param a the source atom for this bond.
   *  @param b the target atom for this bond.
   */
  Bond::Bond(Atom_p a, Atom_p b) : utils::CountableObject<Bond>(), mol_() {
    atoms_[0] = Atom_wp(a);
    atoms_[1] = Atom_wp(b);
    idx_ = GetUniqueID();
  }
  
  Bond::~Bond() { }
  
  // Data retrival methods
  uid_t Bond::GetIndex() const { return idx_; }
  Molecule_p Bond::GetMolecule() const {
    if (!mol_.expired()) return mol_.lock();
    else return Molecule_p();
  }
  BondOrder Bond::GetOrder() const { return order_; }
  Atom_p Bond::GetSourceAtom() const {
    if (!atoms_[0].expired()) return atoms_[0].lock();
    else return Atom_p();
  }
  Atom_p Bond::GetTargetAtom() const {
    if (!atoms_[1].expired()) return atoms_[1].lock();
    else return Atom_p();
  }
  String Bond::ToString() const {
    std::stringstream ss;
    ss << "Bond(";
    if (!atoms_[0].expired()) {
      Atom_p tmp = atoms_[0].lock();
      ss << tmp->GetIndex() << ":\"" << tmp->GetName() << "\"";
    } else {
      ss << "\"None\"";
    }
    ss << ",";
    if (!atoms_[1].expired()) {
      Atom_p tmp = atoms_[1].lock();
      ss << tmp->GetIndex() << ":\"" << tmp->GetName() << "\"";
    } else {
      ss << "\"None\"";
    }
    ss << ")";
    return ss.str();
  }
  
  // Data modification methods
  /** @details Index is used for internal bookkeeping and should not be
   *  considered stable.
   *  @param i the index to set.
   */
  void Bond::SetIndex(uid_t i) { idx_ = i; }
  
  /** @param m the molecule to set. */
  void Bond::SetMolecule(Molecule_p m) { mol_ = m; }
  
  /** @param o the bond order to set this bond to. */
  void Bond::SetOrder(BondOrder o) { order_ = o; }
  
  /** @param a the atom to set this bond's source atom to. */
  void Bond::SetSourceAtom(Atom_p a) { atoms_[0] = Atom_wp(a); }
  
  /** @param a the atom to set this bond's target atom to. */
  void Bond::SetTargetAtom(Atom_p a) { atoms_[1] = Atom_wp(a); }
  
  void Bond::Clear() {
    mol_.reset();
    atoms_.fill(Atom_wp());
  }
  
  // Iterator methods
  Atom_p Bond::Begin(BondAtomIterator &it) {
    it = atoms_.begin();
    return (it == atoms_.end()) ? Atom_p() : it->lock();
  }
  
  Atom_p Bond::Next(BondAtomIterator &it) {
    ++it;
    return (it == atoms_.end()) ? Atom_p() : it->lock();
  }
  
  BondAtomIterator Bond::BeginAtom() { return atoms_.begin(); }
  BondAtomIterator Bond::EndAtom() { return atoms_.end(); }
  
}

