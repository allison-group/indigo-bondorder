//
//  interface.hpp
//  indigo-bondorder
//
//  Created by Welsh, Ivan on 8/01/18.
//  Copyright © 2018 Allison Group. All rights reserved.
//
#include "../api.hpp"
#include <pybind11/pybind11.h>

#ifndef INDIGO_BONDORDER_PYTHON_INTERFACE_HPP
#define INDIGO_BONDORDER_PYTHON_INTERFACE_HPP

/// @todo add opaque set<string> stuff so can add to string
// PYBIND11_MAKE_OPAQUE(std::set<indigo-bondorder::String>)

namespace indigox {
  void GenerateOptions(pybind11::module& m);
  void GeneratePyAtom(pybind11::module& m);
  void GeneratePyBond(pybind11::module& m);
  void GeneratePyMolecule(pybind11::module& m);
  void GeneratePyPeriodicTable(pybind11::module& m);
  void GeneratePyElement(pybind11::module& m);
}

#endif /* INDIGO_BONDORDER_PYTHON_INTERFACE_HPP */
