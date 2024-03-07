//
// Copyright (c) 2017-2024 CNRS INRIA
// This file was borrowed from the Pinocchio library:
// https://github.com/stack-of-tasks/pinocchio
//

#ifndef HPP_FCL_SERIALIZATION_SERIALIZABLE_H
#define HPP_FCL_SERIALIZATION_SERIALIZABLE_H

#include "hpp/fcl/serialization/archive.h"

namespace hpp {
namespace fcl {
namespace serialization {

template <class Derived>
struct Serializable {
 private:
  Derived& derived() { return *static_cast<Derived*>(this); }
  const Derived& derived() const { return *static_cast<const Derived*>(this); }

 public:
  /// \brief Loads a Derived object from a text file.
  void loadFromText(const std::string& filename) {
    ::hpp::fcl::serialization::loadFromText(this->derived(), filename);
  }

  /// \brief Saves a Derived object as a text file.
  void saveToText(const std::string& filename) const {
    ::hpp::fcl::serialization::saveToText(this->derived(), filename);
  }

  /// \brief Loads a Derived object from a stream string.
  void loadFromStringStream(std::istringstream& is) {
    ::hpp::fcl::serialization::loadFromStringStream(this->derived(), is);
  }

  /// \brief Saves a Derived object to a string stream.
  void saveToStringStream(std::stringstream& ss) const {
    ::hpp::fcl::serialization::saveToStringStream(this->derived(), ss);
  }

  /// \brief Loads a Derived object from a  string.
  void loadFromString(const std::string& str) {
    ::hpp::fcl::serialization::loadFromString(this->derived(), str);
  }

  /// \brief Saves a Derived object to a string.
  std::string saveToString() const {
    return ::hpp::fcl::serialization::saveToString(this->derived());
  }

  /// \brief Loads a Derived object from an XML file.
  void loadFromXML(const std::string& filename, const std::string& tag_name) {
    ::hpp::fcl::serialization::loadFromXML(this->derived(), filename, tag_name);
  }

  /// \brief Saves a Derived object as an XML file.
  void saveToXML(const std::string& filename,
                 const std::string& tag_name) const {
    ::hpp::fcl::serialization::saveToXML(this->derived(), filename, tag_name);
  }

  /// \brief Loads a Derived object from an binary file.
  void loadFromBinary(const std::string& filename) {
    ::hpp::fcl::serialization::loadFromBinary(this->derived(), filename);
  }

  /// \brief Saves a Derived object as an binary file.
  void saveToBinary(const std::string& filename) const {
    ::hpp::fcl::serialization::saveToBinary(this->derived(), filename);
  }

  /// \brief Loads a Derived object from a binary container.
  void loadFromBinary(boost::asio::streambuf& container) {
    ::hpp::fcl::serialization::loadFromBinary(this->derived(), container);
  }

  /// \brief Saves a Derived object as a binary container.
  void saveToBinary(boost::asio::streambuf& container) const {
    ::hpp::fcl::serialization::saveToBinary(this->derived(), container);
  }
};

}  // namespace serialization
}  // namespace fcl
}  // namespace hpp

#endif  // ifndef HPP_FCL_SERIALIZATION_SERIALIZABLE_H
