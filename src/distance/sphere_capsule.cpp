/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2021-2022, CNRS
 *  Author: Louis Montaut
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

#include "coal/math/transform.h"
#include "coal/shape/geometric_shapes.h"

#include "coal/internal/shape_shape_func.h"
#include "../narrowphase/details.h"

#include "coal/tracy.hh"

namespace coal {
struct GJKSolver;

namespace internal {
template <>
Scalar ShapeShapeDistance<Sphere, Capsule>(
    const CollisionGeometry* o1, const Transform3s& tf1,
    const CollisionGeometry* o2, const Transform3s& tf2, const GJKSolver*,
    const bool, Vec3s& p1, Vec3s& p2, Vec3s& normal) {
  COAL_TRACY_ZONE_SCOPED_N(
      "coal::internal::ShapeShapeDistance<Sphere, Capsule>");
  const Sphere& s1 = static_cast<const Sphere&>(*o1);
  const Capsule& s2 = static_cast<const Capsule&>(*o2);
  return details::sphereCapsuleDistance(s1, tf1, s2, tf2, p1, p2, normal);
}

template <>
Scalar ShapeShapeDistance<Capsule, Sphere>(
    const CollisionGeometry* o1, const Transform3s& tf1,
    const CollisionGeometry* o2, const Transform3s& tf2, const GJKSolver*,
    const bool, Vec3s& p1, Vec3s& p2, Vec3s& normal) {
  COAL_TRACY_ZONE_SCOPED_N(
      "coal::internal::ShapeShapeDistance<Capsule, Sphere>");
  const Capsule& s1 = static_cast<const Capsule&>(*o1);
  const Sphere& s2 = static_cast<const Sphere&>(*o2);
  const Scalar distance =
      details::sphereCapsuleDistance(s2, tf2, s1, tf1, p2, p1, normal);
  normal = -normal;
  return distance;
}

}  // namespace internal

}  // namespace coal
