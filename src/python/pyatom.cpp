//
//  pyatom.cpp
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
#include "classes/molecule.hpp"
#include "classes/periodictable.hpp"
#include "utils/helpers.hpp"

namespace py = pybind11;

namespace indigox {
  void GeneratePyAtom(pybind11::module& m) {
    py::class_<Atom, Atom_p>(m, "Atom")
    // Hidden methods
    .def(py::init<>(&CreateAtom))
    .def("__str__", &Atom::ToString)
    .def("__repr__", &Atom::ToString)
    // Getters
    .def("GetElement", &Atom::GetElement)
    .def("GetFormalCharge", &Atom::GetFormalCharge)
    .def("GetIndex", &Atom::GetIndex)
    .def("GetMolecule", &Atom::GetMolecule)
    .def("GetName", &Atom::GetName)
    .def("GetX", &Atom::GetX)
    .def("GetY", &Atom::GetY)
    .def("GetZ", &Atom::GetZ)
    .def("GetUniqueID", &Atom::GetUniqueID)
    // Setters
    .def("SetElement", py::overload_cast<Element_p>(&Atom::SetElement))
    .def("SetElement", py::overload_cast<String>(&Atom::SetElement))
    .def("SetElement", py::overload_cast<Int>(&Atom::SetElement))
    .def("SetFormalCharge", &Atom::SetFormalCharge)
    .def("SetIndex", &Atom::SetIndex)
    .def("SetMolecule", &Atom::SetMolecule)
    .def("SetName", &Atom::SetName)
    .def("SetPosition", &Atom::SetPosition)
    ;
  }
}
