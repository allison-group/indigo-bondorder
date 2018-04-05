//
//  pybond.cpp
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
#include "utils/helpers.hpp"

namespace py = pybind11;

namespace indigox {
  void GeneratePyBond(py::module& m) {
    // Bond class
    py::class_<Bond, Bond_p>(m, "Bond")
    // Hidden methods
    .def(py::init<>(&CreateBond))
    .def("__str__", &Bond::ToString)
    .def("__repr__", &Bond::ToString)
    // Getters
    .def("GetIndex", &Bond::GetIndex)
    .def("GetMolecule", &Bond::GetMolecule)
    .def("GetOrder", &Bond::GetOrder)
    .def("GetSourceAtom", &Bond::GetSourceAtom)
    .def("GetTargetAtom", &Bond::GetTargetAtom)
    // Setters
    .def("SetIndex", &Bond::SetIndex)
    .def("SetMolecule", &Bond::SetMolecule)
    .def("SetOrder", &Bond::SetOrder)
    .def("SetSourceAtom", &Bond::SetSourceAtom)
    .def("SetTargetAtom", &Bond::SetTargetAtom)
    ;
    
    // BondOrder enum
    py::enum_<BondOrder>(m, "BondOrder")
    .value("UNDEFINED_BOND", BondOrder::UNDEFINED_BOND)
    .value("SINGLE_BOND", BondOrder::SINGLE_BOND)
    .value("DOUBLE_BOND", BondOrder::DOUBLE_BOND)
    .value("TRIPLE_BOND", BondOrder::TRIPLE_BOND)
    .value("QUADRUPLE_BOND", BondOrder::QUADRUPLE_BOND)
    .value("AROMATIC_BOND", BondOrder::AROMATIC_BOND)
    .value("ONEANDAHALF_BOND", BondOrder::ONEANDAHALF_BOND)
    .value("TWOANDAHALF_BOND", BondOrder::TWOANDAHALF_BOND)
    ;
  }
}
