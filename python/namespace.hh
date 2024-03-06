//
// Copyright (c) 2017-2024 CNRS INRIA
// This file was borrowed from the Pinocchio library:
// https://github.com/stack-of-tasks/pinocchio
//

#ifndef HPP_FCL_PYTHON_NAMESPACE_H
#define HPP_FCL_PYTHON_NAMESPACE_H

#include <boost/python.hpp>
#include <string>

namespace hpp {
namespace fcl {
namespace python {

///
///  \returns the name of the current Python scope.
///
inline std::string getCurrentScopeName() {
  namespace bp = boost::python;
  bp::scope current_scope;

  return std::string(bp::extract<const char*>(current_scope.attr("__name__")));
}

///
///  \brief Helper to create or simply return an existing namespace in Python
///
/// \param[in] submodule_name name of the submodule
///
/// \returns The submodule related to the namespace name.
///
inline boost::python::object getOrCreatePythonNamespace(
    const std::string& submodule_name) {
  namespace bp = boost::python;

  const std::string complete_submodule_name =
      getCurrentScopeName() + "." + submodule_name;

  bp::object submodule(
      bp::borrowed(PyImport_AddModule(complete_submodule_name.c_str())));
  bp::scope current_scope;
  current_scope.attr(submodule_name.c_str()) = submodule;

  return submodule;
}

}  // namespace python
}  // namespace fcl
}  // namespace hpp

#endif  // HPP_FCL_PYTHON_NAMESPACE_H
