//
//  iters.cpp
//  indigox
//
//  Created by Ivan Welsh on 7/01/18.
//  Copyright © 2018 Hermes Productions. All rights reserved.
//
#include <iostream>

#include "api.hpp"
#include "classes/iters.hpp"
#include "classes/molecule.hpp"

using namespace indigox;

AtomIterator::AtomIterator() : type_(UNDEFINED), ptr_() { }

AtomIterator::AtomIterator(Molecule_p m) : type_(MOLECULE_ATOM) {
  parentMol_ = m;
  ptr_ = parentMol_->Begin(itMol_);
}

AtomIterator::AtomIterator(Bond_p b) : type_(BOND_ATOM) {
  parentBond_ = b;
  ptr_ = parentBond_->Begin(itBond_);
}

AtomIterator::AtomIterator(Atom_p a) : type_(ATOM_NEIGHBOUR) {
  parentAtom_ = a;
  Bond_p tmp = parentAtom_->Begin(itAtom_);
  ptr_ = tmp->GetSourceAtom();
  if (!ptr_ || ptr_->GetUniqueID() != parentAtom_->GetUniqueID()) {
    ptr_ = tmp->GetTargetAtom();
  }
}

AtomIterator::AtomIterator(const AtomIterator& it) {
  ResetToOther(it);
}

AtomIterator::~AtomIterator() { }

AtomIterator& AtomIterator::operator=(const AtomIterator& it) {
  if (this != &it) ResetToOther(it);
  return *this;
}

void AtomIterator::ResetToOther(const AtomIterator& it) {
  type_ = it.type_;
  itMol_ = it.itMol_;
  parentMol_ = it.parentMol_;
  ptr_ = it.ptr_;
}

AtomIterator::operator bool() const { return bool(ptr_); }

AtomIterator& AtomIterator::operator++() {
  switch (type_) {
    case BOND_ATOM:
      ptr_ = parentBond_->Next(itBond_);
      break;
    case MOLECULE_ATOM:
      ptr_ = parentMol_->Next(itMol_);
      break;
    case ATOM_NEIGHBOUR:
      {
        Bond_p tmp = parentAtom_->Next(itAtom_);
        ptr_ = tmp->GetSourceAtom();
        if (!ptr_ || ptr_->GetUniqueID() != parentAtom_->GetUniqueID()) {
          ptr_ = tmp->GetTargetAtom();
        }
      }
      break;
      
    default:
      break;
  }
  return *this;
}

AtomIterator AtomIterator::operator++(int) {
  AtomIterator tmp(*this);
  operator++();
  return tmp;
}

Atom_p AtomIterator::operator->() const { return ptr_; }
Atom_p AtomIterator::operator*() const { return ptr_; }
