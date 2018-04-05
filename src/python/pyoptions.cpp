#include "api.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "python/interface.hpp"

#include "utils/options.hpp"

namespace py = pybind11;
using namespace indigox;

void indigox::GenerateOptions(py::module& m) {
  // General options
  typedef Options o_;
  py::class_<o_, std::shared_ptr<o_>> PyOptions(m, "Options");
  PyOptions.def_static("Reset", &o_::Reset)
  .def_readwrite_static("DATA_DIRECTORY", &o_::DATA_DIRECTORY)
  .def_readwrite_static("PERIODIC_TABLE_FILE", &o_::PERIODIC_TABLE_FILE)
  ;
  
  // AssignElectrons general options
  typedef o_::AssignElectrons e_;
  py::class_<e_, std::shared_ptr<e_>> PyAssElec(PyOptions, "AssignElectrons");
  py::enum_<e_::Algorithm>(PyAssElec, "Algorithm")
  .value("LOCAL_OPTIMISATION", e_::Algorithm::LOCAL_OPTIMISATION)
  .value("ASTAR", e_::Algorithm::ASTAR)
  .value("FPT", e_::Algorithm::FPT)
  ;
  PyAssElec.def_static("Reset", &e_::Reset)
  .def_readwrite_static("ALGORITHM", &e_::ALGORITHM)
  .def_readwrite_static("ATOM_ENERGY_FILE", &e_::ATOM_ENERGY_FILE)
  .def_readwrite_static("BOND_ENERGY_FILE", &e_::BOND_ENERGY_FILE)
  .def_readwrite_static("INF", &e_::INF)
  .def_readwrite_static("MAXIMUM_BOND_ORDER", &e_::MAXIMUM_BOND_ORDER)
  .def_readwrite_static("USE_ELECTRON_PAIRS", &e_::USE_ELECTRON_PAIRS)
  .def_readwrite_static("AUTO_USE_ELECTRON_PAIRS",&e_::AUTO_USE_ELECTRON_PAIRS)
  .def_readwrite_static("USE_CHARGED_BOND_ENERGIES", &e_::USE_CHARGED_BOND_ENERGIES)
  .def_readwrite_static("ALLOWED_ELEMENTS", &e_::ALLOWED_ELEMENTS)
  .def_readwrite_static("HIGHEST_MAGNITUDE_CHARGE", &e_::HIGHEST_MAGNITUDE_CHARGE)
  .def_readwrite_static("ALLOW_CHARGED_CARBON", &e_::ALLOW_CHARGED_CARBON)
  .def_readwrite_static("MAXIMUM_RESULT_COUNT", &e_::MAXIMUM_RESULT_COUNT)
  .def_readwrite_static("PREPLACE_ELECTRONS", &e_::PREPLACE_ELECTRONS)
  ;
  
  // AStar specific options
  typedef e_::AStar a_;
  py::class_<a_, std::shared_ptr<a_>> PyAStar(PyAssElec, "AStar");
  py::enum_<a_::Heuristic>(PyAStar, "Heuristic")
  .value("PROMISCUOUS", a_::Heuristic::PROMISCUOUS)
  .value("ABSTEMIOUS", a_::Heuristic::ABSTEMIOUS)
  ;
  PyAStar.def_static("Reset", &a_::Reset)
  .def_readwrite_static("HEURISTIC", &a_::HEURISTIC)
  .def_readwrite_static("MEGABYTE_LIMIT", &a_::MEGABYTE_LIMIT)
  ;
  
  // FPT specific options
  typedef e_::FPT f_;
  py::class_<f_, std::shared_ptr<f_>> PyFPT(PyAssElec, "FPT");
  py::enum_<f_::PermAlgo>(PyFPT, "PermutationAlgorithm")
  .value("RANDOM", f_::PermAlgo::RANDOM)
  .value("QUICKBB", f_::PermAlgo::QUICKBB)
  .value("MINDEGREE", f_::PermAlgo::MINDEGREE)
  .value("MINADDEDGES", f_::PermAlgo::MINADDEDGES)
  ;
  PyFPT.def_static("Reset", &f_::Reset)
  .def_readwrite_static("LIBTW_JAR_FILE", &f_::LIBTW_JAR_FILE)
  .def_readwrite_static("ADD_EDGES_TO_TD", &f_::ADD_EDGES_TO_TD)
  .def_readwrite_static("PERMUTATION_ALGORITM", &f_::PERM_ALGO)
  .def_readwrite_static("MINIMUM_PROPAGATION_DEPTH", &f_::MINIMUM_PROPAGATION_DEPTH)
  ;
  
  // LO specific options
  typedef e_::LocalOptimisation l_;
  py::class_<l_, std::shared_ptr<l_>> PyLO(PyAssElec, "LocalOptimisation");
  PyLO.def_static("Reset", &l_::Reset)
  .def_readwrite_static("OPTIMISE_ALL_MINIMUMS", &l_::OPTIMISE_ALL_MINIMUMS)
  .def_readwrite_static("CACHE_RESULTS", &l_::CACHE_RESULTS)
  .def_readwrite_static("CACHE_INFINITIES", &l_::CACHE_INFINITIES)
  .def_readwrite_static("TIMEOUT_LIMIT", &l_::TIMEOUT_LIMIT)
  ;
}
