//
//  helpers.cpp
//  indigox
//
//  Created by Welsh, Ivan on 8/01/18.
//  Copyright © 2018 Allison Group. All rights reserved.
//

#include "api.hpp"
#include "classes/atom.hpp"
#include "classes/bond.hpp"
#include "classes/molecule.hpp"
#include "utils/helpers.hpp"

namespace indigox {
  /// @returns a shared pointer to a new Atom instance.
  Atom_p CreateAtom() { return Atom_p(new Atom()); }
  
  /// @returns a shared pointer to a new Bond instance.
  Bond_p CreateBond() { return Bond_p(new Bond()); }
  
  /// @returns a shared pointer to a new Molecule instance.
  Molecule_p CreateMolecule() { return Molecule_p(new Molecule()); }
  
}
