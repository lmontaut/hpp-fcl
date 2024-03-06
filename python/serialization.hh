//
// Copyright (c) 2017-2024 CNRS INRIA
// This file was borrowed from the Pinocchio library:
// https://github.com/stack-of-tasks/pinocchio
//

#ifndef HPP_FCL_PYTHON_SERIALIZATION_H
#define HPP_FCL_PYTHON_SERIALIZATION_H

#include "hpp/fcl/serialization/archive.h"
#include "namespace.hh"

namespace hpp {
namespace fcl {
namespace python {

template <typename T>
void serialize() {
  namespace bp = boost::python;

  bp::scope current_scope = getOrCreatePythonNamespace("serialization");

  bp::def("loadFromBinary",
          (void (*)(T&, boost::asio::streambuf&))
              hpp::fcl::serialization::loadFromBinary<T>,
          bp::args("object", "stream_buffer"),
          "Load an object from a binary buffer.");

  bp::def("saveToBinary",
          (void (*)(const T&, boost::asio::streambuf&))
              hpp::fcl::serialization::saveToBinary<T>,
          bp::args("object", "stream_buffer"),
          "Save an object to a binary buffer.");

  bp::def("loadFromBinary",
          (void (*)(T&, serialization::StaticBuffer&))
              hpp::fcl::serialization::loadFromBinary<T>,
          bp::args("object", "static_buffer"),
          "Load an object from a static binary buffer.");

  bp::def("saveToBinary",
          (void (*)(const T&, serialization::StaticBuffer&))
              hpp::fcl::serialization::saveToBinary<T>,
          bp::args("object", "static_buffer"),
          "Save an object to a static binary buffer.");
}

}  // namespace python
}  // namespace fcl
}  // namespace hpp

#endif
