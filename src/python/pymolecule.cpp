//
//  pymolecule.cpp
//  indigox
//
//  Created by Welsh, Ivan on 8/01/18.
//  Copyright © 2018 Allison Group. All rights reserved.
//

#include "api.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "python/interface.hpp"

#include "classes/atom.hpp"
#include "classes/bond.hpp"
#include "classes/molecule.hpp"
#include "classes/periodictable.hpp"
#include "utils/helpers.hpp"

namespace py = pybind11;
namespace indigox {
  void GeneratePyMolecule(py::module& m) {
    py::class_<Molecule, Molecule_p>(m, "Molecule")
    // Hidden methods
    .def(py::init<>(&CreateMolecule))
    // Data retrival
    .def("GetAtomIndex", &Molecule::GetAtomIndex)
    .def("GetAtomUniqueID", &Molecule::GetAtomUniqueID)
    .def("GetBondIndex", &Molecule::GetBondIndex)
    .def("GetBond", &Molecule::GetBond)
    .def("GetFormula", &Molecule::GetFormula)
    .def("GetName", &Molecule::GetName)
    .def("GetTotalCharge", &Molecule::GetTotalCharge)
    .def("NumAtoms", &Molecule::NumAtoms)
    .def("NumBonds", &Molecule::NumBonds)
    // Data modification
    .def("SetName", &Molecule::SetName)
    .def("SetTotalCharge", &Molecule::SetTotalCharge)
    // Molecule modification
    .def("NewAtom", py::overload_cast<>(&Molecule::NewAtom))
    .def("NewAtom", py::overload_cast<Element_p>(&Molecule::NewAtom))
    .def("NewAtom", py::overload_cast<uid_t,Element_p>(&Molecule::NewAtom))
    .def("NewBond", &Molecule::NewBond)
    .def("RemoveAtom", &Molecule::RemoveAtom)
    .def("RemoveBond", &Molecule::RemoveBond)
    // Structural determination
    .def("AssignElectrons", &Molecule::AssignElectrons)
    .def("ApplyElectronAssignment", &Molecule::ApplyElectronAssignment)
    .def("GetMinimumElectronAssignmentScore", &Molecule::GetMinimumElectronAssignmentScore)
    // Iterators
    .def("GetAtoms", [](const Molecule_p m){
      return py::make_iterator(m->BeginAtom(), m->EndAtom());
    }, py::keep_alive<0, 1>())
    .def("GetBonds", [](const Molecule_p m){
      return py::make_iterator(m->BeginBond(), m->EndBond());
    }, py::keep_alive<0, 1>())
    ;
  }
}

