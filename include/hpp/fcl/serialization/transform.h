//
// Copyright (c) 2022 INRIA
//

#ifndef HPP_FCL_SERIALIZATION_TRANSFORM_H
#define HPP_FCL_SERIALIZATION_TRANSFORM_H

#include "hpp/fcl/math/transform.h"
#include "hpp/fcl/serialization/fwd.h"

namespace boost {
namespace serialization {

template <class Archive>
void save(Archive& ar, const hpp::fcl::Transform3f& tf,
          const unsigned int /*version*/) {
  ar& make_nvp("translation", make_array(tf.translation().data(), 3));
  ar& make_nvp("rotation", make_array(tf.rotation().data(), 9));
}

template <class Archive>
void load(Archive& ar, hpp::fcl::Transform3f& tf,
          const unsigned int /*version*/) {
  ar >> make_nvp("translation", make_array(tf.translation().data(), 3));
  ar >> make_nvp("rotation", make_array(tf.rotation().data(), 9));
}

template <class Archive>
void serialize(Archive& ar, hpp::fcl::Transform3f& tf,
               const unsigned int version) {
  split_free(ar, tf, version);
}

}  // namespace serialization
}  // namespace boost

#endif  // ifndef HPP_FCL_SERIALIZATION_TRANSFORM_H
