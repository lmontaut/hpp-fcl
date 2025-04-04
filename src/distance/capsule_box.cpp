#include "coal/math/transform.h"
#include "coal/shape/geometric_shapes.h"

#include "coal/internal/shape_shape_func.h"
#include "../narrowphase/details.h"

#include "coal/tracy.hh"

namespace coal {

namespace internal {

template <>
Scalar ShapeShapeDistance<Box, Capsule>(const CollisionGeometry* o1,
                                        const Transform3s& tf1,
                                        const CollisionGeometry* o2,
                                        const Transform3s& tf2,
                                        const GJKSolver*, const bool, Vec3s& p1,
                                        Vec3s& p2, Vec3s& normal) {
  const Box& s1 = static_cast<const Box&>(*o1);
  const Capsule& s2 = static_cast<const Capsule&>(*o2);
  const Scalar dist =
      details::capsuleBoxDistance(s2, tf2, s1, tf1, p2, p1, normal);
  normal *= -1;
  return dist;
}

template <>
Scalar ShapeShapeDistance<Capsule, Box>(const CollisionGeometry* o1,
                                        const Transform3s& tf1,
                                        const CollisionGeometry* o2,
                                        const Transform3s& tf2,
                                        const GJKSolver*, const bool, Vec3s& p1,
                                        Vec3s& p2, Vec3s& normal) {
  const Capsule& s1 = static_cast<const Capsule&>(*o1);
  const Box& s2 = static_cast<const Box&>(*o2);
  const Scalar dist =
      details::capsuleBoxDistance(s1, tf1, s2, tf2, p1, p2, normal);
  return dist;
}

}  // namespace internal
}  // namespace coal
