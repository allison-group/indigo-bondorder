#include "api.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "python/interface.hpp"

#include "classes/periodictable.hpp"

namespace py = pybind11;
namespace indigox {
  void GeneratePyElement(py::module& m) {
    py::class_<Element, Element_p>(m, "Element")
    .def("__eq__", py::overload_cast<Element_p, uint8_t>(&operator==))
    .def("__eq__", py::overload_cast<Element_p, String>(&operator==))
    .def("__eq__", py::overload_cast<Element_p, Element_p>(&operator==))
    .def("__ne__", py::overload_cast<Element_p, uint8_t>(&operator!=))
    .def("__ne__", py::overload_cast<Element_p, String>(&operator!=))
    .def("__ne__", py::overload_cast<Element_p, Element_p>(&operator!=))
    .def("__str__", &Element::ToString)
    .def("__repr__", &Element::ToString)
    .def("GetAtomicMass", &Element::GetAtomicMass)
    .def("GetAtomicNumber", &Element::GetAtomicNumber)
    .def("GetAtomicRadius", &Element::GetAtomicRadius)
    .def("GetCovalentRadius", &Element::GetCovalentRadius)
    .def("GetVanDerWaalsRadius", &Element::GetVanDerWaalsRadius)
    .def("GetName", &Element::GetName)
    .def("GetSymbol", &Element::GetSymbol)
    .def("GetGroup", &Element::GetGroup)
    .def("GetPeriod", &Element::GetPeriod)
    .def("GetValenceElectronCount", &Element::GetValenceElectronCount)
    .def("GetOctet", &Element::GetOctet)
    .def("GetHypervalentOctet", &Element::GetHypervalentOctet)
    .def("GetElectronegativity", &Element::GetElectronegativity)
    ;
  }
}
