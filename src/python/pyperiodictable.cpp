#include "api.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "python/interface.hpp"

#include "classes/periodictable.hpp"

namespace py = pybind11;
namespace indigox {
  void GeneratePyPeriodicTable(py::module& m) {
    py::class_<PeriodicTable, PeriodicTable_p>(m, "PeriodicTable")
    .def(py::init<>(&PeriodicTable::GetInstance), py::return_value_policy::reference)
    .def("GetInstance", &PeriodicTable::GetInstance, py::return_value_policy::reference)
    .def("GetElement", py::overload_cast<uint8_t>(&PeriodicTable::GetElement, py::const_))
    .def("GetElement", py::overload_cast<String>(&PeriodicTable::GetElement, py::const_))
    .def("__getitem__", py::overload_cast<uint8_t>(&PeriodicTable::GetElement, py::const_))
    .def("__getitem__", py::overload_cast<String>(&PeriodicTable::GetElement, py::const_))
    .def("NumElements", &PeriodicTable::NumElements)
    ;
  }
}

