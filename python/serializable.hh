//
// Copyright (c) 2017-2024 CNRS INRIA
// This file was borrowed from the Pinocchio library:
// https://github.com/stack-of-tasks/pinocchio
//

#ifndef HPP_FCL_PYTHON_SERIALIZABLE_H
#define HPP_FCL_PYTHON_SERIALIZABLE_H

#include "serialization.hh"
#include "hpp/fcl/serialization/serializable.h"

namespace hpp {
namespace fcl {
namespace python {

namespace bp = boost::python;

template <typename Derived>
struct SerializableVisitor
    : public bp::def_visitor<SerializableVisitor<Derived> > {
  template <class PyClass>
  void visit(PyClass& cl) const {
    cl.def("saveToText", &Derived::saveToText, bp::arg("filename"),
           "Saves *this inside a text file.")
        .def("loadFromText", &Derived::loadFromText, bp::arg("filename"),
             "Loads *this from a text file.")

        .def("saveToString", &Derived::saveToString, bp::arg("self"),
             "Parses the current object to a string.")
        .def("loadFromString", &Derived::loadFromString,
             bp::args("self", "string"),
             "Parses from the input string the content of the current object.")

        .def("saveToXML", &Derived::saveToXML, bp::args("filename", "tag_name"),
             "Saves *this inside a XML file.")
        .def("loadFromXML", &Derived::loadFromXML,
             bp::args("self", "filename", "tag_name"),
             "Loads *this from a XML file.")

        .def("saveToBinary",
             (void(Derived::*)(const std::string&) const) &
                 Derived::saveToBinary,
             bp::args("self", "filename"), "Saves *this inside a binary file.")
        .def("loadFromBinary",
             (void(Derived::*)(const std::string&)) & Derived::loadFromBinary,
             bp::args("self", "filename"), "Loads *this from a binary file.")

        .def("saveToBinary",
             (void(Derived::*)(boost::asio::streambuf&) const) &
                 Derived::saveToBinary,
             bp::args("self", "buffer"), "Saves *this inside a binary buffer.")
        .def("loadFromBinary",
             (void(Derived::*)(boost::asio::streambuf&)) &
                 Derived::loadFromBinary,
             bp::args("self", "buffer"), "Loads *this from a binary buffer.")

        .def("saveToBinary",
             (void(Derived::*)(serialization::StaticBuffer&) const) &
                 Derived::saveToBinary,
             bp::args("self", "buffer"),
             "Saves *this inside a static binary buffer.")
        .def("loadFromBinary",
             (void(Derived::*)(serialization::StaticBuffer&)) &
                 Derived::loadFromBinary,
             bp::args("self", "buffer"),
             "Loads *this from a static binary buffer.");
    serialize<Derived>();
  }
};
}  // namespace python
}  // namespace fcl
}  // namespace hpp

#endif  // ifndef HPP_FCL_PYTHON_SERIALIZABLE_H
