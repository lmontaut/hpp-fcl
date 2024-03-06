//
// Copyright (c) 2021-2024 CNRS INRIA
// This file was borrowed from the Pinocchio library:
// https://github.com/stack-of-tasks/pinocchio
//

#ifndef HPP_FCL_SERIALIZATION_STATIC_BUFFER_H
#define HPP_FCL_SERIALIZATION_STATIC_BUFFER_H

#include <vector>

namespace hpp {
namespace fcl {
namespace serialization {

/// \brief Static buffer with pre-allocated memory.
struct StaticBuffer {
  ///  \brief Defautl constructor from a given size
  explicit StaticBuffer(const size_t n) : m_size(n) { this->m_data.reserve(n); }

  /// \brief Returns the current size of the buffer
  size_t size() const { return this->m_size; }

  /// \brief Returns the pointer on the data
  char* data() { return this->m_data.data(); }

  /// \brief Returns the pointer on the data (const version)
  const char* data() const { return this->m_data.data(); }

  /// \brief Increase the capacity of the vector to a value that's greater or
  /// equal to new_size.
  ///
  /// \param[in] new_size New capacity of the buffer.
  ///
  void resize(const size_t new_size) {
    this->m_size = new_size;
    this->m_data.reserve(new_size);
  }

 protected:
  size_t m_size;
  std::vector<char> m_data;
};

}  // namespace serialization
}  // namespace fcl
}  // namespace hpp

#endif  // ifndef HPP_FCL_SERIALIZATION_STATIC_BUFFER_H
