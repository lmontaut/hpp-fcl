/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2011-2014, Willow Garage, Inc.
 *  Copyright (c) 2014-2015, Open Source Robotics Foundation
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Open Source Robotics Foundation nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

/** \author Jia Pan */

#ifndef COAL_BV_SPLITTER_H
#define COAL_BV_SPLITTER_H

#include "coal/BVH/BVH_internal.h"
#include "coal/BV/kIOS.h"
#include "coal/BV/OBBRSS.h"
#include <vector>
#include <iostream>

namespace coal {

/// @brief Three types of split algorithms are provided in FCL as default
enum SplitMethodType {
  SPLIT_METHOD_MEAN,
  SPLIT_METHOD_MEDIAN,
  SPLIT_METHOD_BV_CENTER
};

/// @brief A class describing the split rule that splits each BV node
template <typename BV>
class BVSplitter {
 public:
  BVSplitter(SplitMethodType method)
      : split_vector(0, 0, 0), split_method(method) {}

  /// @brief Default deconstructor
  virtual ~BVSplitter() {}

  /// @brief Set the geometry data needed by the split rule
  void set(Vec3s* vertices_, Triangle32* tri_indices_, BVHModelType type_) {
    vertices = vertices_;
    tri_indices = tri_indices_;
    type = type_;
  }

  /// @brief Compute the split rule according to a subset of geometry and the
  /// corresponding BV node
  void computeRule(const BV& bv, unsigned int* primitive_indices,
                   unsigned int num_primitives) {
    switch (split_method) {
      case SPLIT_METHOD_MEAN:
        computeRule_mean(bv, primitive_indices, num_primitives);
        break;
      case SPLIT_METHOD_MEDIAN:
        computeRule_median(bv, primitive_indices, num_primitives);
        break;
      case SPLIT_METHOD_BV_CENTER:
        computeRule_bvcenter(bv, primitive_indices, num_primitives);
        break;
      default:
        std::cerr << "Split method not supported" << std::endl;
    }
  }

  /// @brief Apply the split rule on a given point
  bool apply(const Vec3s& q) const { return q[split_axis] > split_value; }

  /// @brief Clear the geometry data set before
  void clear() {
    vertices = NULL;
    tri_indices = NULL;
    type = BVH_MODEL_UNKNOWN;
  }

 protected:
  /// @brief The axis based on which the split decision is made. For most BV,
  /// the axis is aligned with one of the world coordinate, so only split_axis
  /// is needed. For oriented node, we can use a vector to make a better split
  /// decision.
  int split_axis;
  Vec3s split_vector;

  /// @brief The split threshold, different primitives are splitted according
  /// whether their projection on the split_axis is larger or smaller than the
  /// threshold
  Scalar split_value;

  /// @brief The mesh vertices or points handled by the splitter
  Vec3s* vertices;

  /// @brief The triangles handled by the splitter
  Triangle32* tri_indices;

  /// @brief Whether the geometry is mesh or point cloud
  BVHModelType type;

  /// @brief The split algorithm used
  SplitMethodType split_method;

  /// @brief Split algorithm 1: Split the node from center
  void computeRule_bvcenter(const BV& bv, unsigned int*, unsigned int) {
    Vec3s center = bv.center();
    int axis = 2;

    if (bv.width() >= bv.height() && bv.width() >= bv.depth())
      axis = 0;
    else if (bv.height() >= bv.width() && bv.height() >= bv.depth())
      axis = 1;

    split_axis = axis;
    split_value = center[axis];
  }

  /// @brief Split algorithm 2: Split the node according to the mean of the data
  /// contained
  void computeRule_mean(const BV& bv, unsigned int* primitive_indices,
                        unsigned int num_primitives) {
    int axis = 2;

    if (bv.width() >= bv.height() && bv.width() >= bv.depth())
      axis = 0;
    else if (bv.height() >= bv.width() && bv.height() >= bv.depth())
      axis = 1;

    split_axis = axis;
    Scalar sum = 0;

    if (type == BVH_MODEL_TRIANGLES) {
      for (unsigned int i = 0; i < num_primitives; ++i) {
        const Triangle32& t = tri_indices[primitive_indices[i]];
        sum += (vertices[t[0]][split_axis] + vertices[t[1]][split_axis] +
                vertices[t[2]][split_axis]);
      }

      sum /= 3;
    } else if (type == BVH_MODEL_POINTCLOUD) {
      for (unsigned int i = 0; i < num_primitives; ++i) {
        sum += vertices[primitive_indices[i]][split_axis];
      }
    }

    split_value = sum / Scalar(num_primitives);
  }

  /// @brief Split algorithm 3: Split the node according to the median of the
  /// data contained
  void computeRule_median(const BV& bv, unsigned int* primitive_indices,
                          unsigned int num_primitives) {
    int axis = 2;

    if (bv.width() >= bv.height() && bv.width() >= bv.depth())
      axis = 0;
    else if (bv.height() >= bv.width() && bv.height() >= bv.depth())
      axis = 1;

    split_axis = axis;
    std::vector<Scalar> proj((size_t)num_primitives);

    if (type == BVH_MODEL_TRIANGLES) {
      for (unsigned int i = 0; i < num_primitives; ++i) {
        const Triangle32& t = tri_indices[primitive_indices[i]];
        proj[i] = (vertices[t[0]][split_axis] + vertices[t[1]][split_axis] +
                   vertices[t[2]][split_axis]) /
                  3;
      }
    } else if (type == BVH_MODEL_POINTCLOUD) {
      for (unsigned int i = 0; i < num_primitives; ++i)
        proj[i] = vertices[primitive_indices[i]][split_axis];
    }

    std::sort(proj.begin(), proj.end());

    if (num_primitives % 2 == 1) {
      split_value = proj[(num_primitives - 1) / 2];
    } else {
      split_value =
          (proj[num_primitives / 2] + proj[num_primitives / 2 - 1]) / 2;
    }
  }
};

template <>
bool COAL_DLLAPI BVSplitter<OBB>::apply(const Vec3s& q) const;

template <>
bool COAL_DLLAPI BVSplitter<RSS>::apply(const Vec3s& q) const;

template <>
bool COAL_DLLAPI BVSplitter<kIOS>::apply(const Vec3s& q) const;

template <>
bool COAL_DLLAPI BVSplitter<OBBRSS>::apply(const Vec3s& q) const;

template <>
void COAL_DLLAPI BVSplitter<OBB>::computeRule_bvcenter(
    const OBB& bv, unsigned int* primitive_indices,
    unsigned int num_primitives);

template <>
void COAL_DLLAPI BVSplitter<OBB>::computeRule_mean(
    const OBB& bv, unsigned int* primitive_indices,
    unsigned int num_primitives);

template <>
void COAL_DLLAPI BVSplitter<OBB>::computeRule_median(
    const OBB& bv, unsigned int* primitive_indices,
    unsigned int num_primitives);

template <>
void COAL_DLLAPI BVSplitter<RSS>::computeRule_bvcenter(
    const RSS& bv, unsigned int* primitive_indices,
    unsigned int num_primitives);

template <>
void COAL_DLLAPI BVSplitter<RSS>::computeRule_mean(
    const RSS& bv, unsigned int* primitive_indices,
    unsigned int num_primitives);

template <>
void COAL_DLLAPI BVSplitter<RSS>::computeRule_median(
    const RSS& bv, unsigned int* primitive_indices,
    unsigned int num_primitives);

template <>
void COAL_DLLAPI BVSplitter<kIOS>::computeRule_bvcenter(
    const kIOS& bv, unsigned int* primitive_indices,
    unsigned int num_primitives);

template <>
void COAL_DLLAPI BVSplitter<kIOS>::computeRule_mean(
    const kIOS& bv, unsigned int* primitive_indices,
    unsigned int num_primitives);

template <>
void COAL_DLLAPI BVSplitter<kIOS>::computeRule_median(
    const kIOS& bv, unsigned int* primitive_indices,
    unsigned int num_primitives);

template <>
void COAL_DLLAPI BVSplitter<OBBRSS>::computeRule_bvcenter(
    const OBBRSS& bv, unsigned int* primitive_indices,
    unsigned int num_primitives);

template <>
void COAL_DLLAPI BVSplitter<OBBRSS>::computeRule_mean(
    const OBBRSS& bv, unsigned int* primitive_indices,
    unsigned int num_primitives);

template <>
void COAL_DLLAPI BVSplitter<OBBRSS>::computeRule_median(
    const OBBRSS& bv, unsigned int* primitive_indices,
    unsigned int num_primitives);

}  // namespace coal

#endif
