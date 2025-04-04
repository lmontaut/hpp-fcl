/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2011-2014, Willow Garage, Inc.
 *  Copyright (c) 2014-2015, Open Source Robotics Foundation
 *  Copyright (c) 2018-2019, Centre National de la Recherche Scientifique
 *  Copyright (c) 2021-2024, INRIA
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
/** \author Jia Pan, Florent Lamiraux */

#ifndef COAL_SRC_NARROWPHASE_DETAILS_H
#define COAL_SRC_NARROWPHASE_DETAILS_H

#include "coal/internal/traversal_node_setup.h"
#include "coal/narrowphase/narrowphase.h"

namespace coal {
namespace details {
// Compute the point on a line segment that is the closest point on the
// segment to to another point. The code is inspired by the explanation
// given by Dan Sunday's page:
//   http://geomalgorithms.com/a02-_lines.html
static inline void lineSegmentPointClosestToPoint(const Vec3s& p,
                                                  const Vec3s& s1,
                                                  const Vec3s& s2, Vec3s& sp) {
  Vec3s v = s2 - s1;
  Vec3s w = p - s1;

  Scalar c1 = w.dot(v);
  Scalar c2 = v.dot(v);

  if (c1 <= 0) {
    sp = s1;
  } else if (c2 <= c1) {
    sp = s2;
  } else {
    Scalar b = c1 / c2;
    Vec3s Pb = s1 + v * b;
    sp = Pb;
  }
}

/// @param p1 witness point on the Sphere.
/// @param p2 witness point on the Capsule.
/// @param normal pointing from shape 1 to shape 2 (sphere to capsule).
/// @return the distance between the two shapes (negative if penetration).
inline Scalar sphereCapsuleDistance(const Sphere& s1, const Transform3s& tf1,
                                    const Capsule& s2, const Transform3s& tf2,
                                    Vec3s& p1, Vec3s& p2, Vec3s& normal) {
  Vec3s pos1(tf2.transform(Vec3s(0., 0., s2.halfLength)));
  Vec3s pos2(tf2.transform(Vec3s(0., 0., -s2.halfLength)));
  Vec3s s_c = tf1.getTranslation();

  Vec3s segment_point;

  lineSegmentPointClosestToPoint(s_c, pos1, pos2, segment_point);
  normal = segment_point - s_c;
  Scalar norm(normal.norm());
  Scalar r1 = s1.radius + s1.getSweptSphereRadius();
  Scalar r2 = s2.radius + s2.getSweptSphereRadius();
  Scalar dist = norm - r1 - r2;

  static const Scalar eps(std::numeric_limits<Scalar>::epsilon());
  if (norm > eps) {
    normal.normalize();
  } else {
    normal << 1, 0, 0;
  }
  p1 = s_c + normal * r1;
  p2 = segment_point - normal * r2;
  return dist;
}

/// @param p1 witness point on the Sphere.
/// @param p2 witness point on the Cylinder.
/// @param normal pointing from shape 1 to shape 2 (sphere to cylinder).
/// @return the distance between the two shapes (negative if penetration).
inline Scalar sphereCylinderDistance(const Sphere& s1, const Transform3s& tf1,
                                     const Cylinder& s2, const Transform3s& tf2,
                                     Vec3s& p1, Vec3s& p2, Vec3s& normal) {
  static const Scalar eps(sqrt(std::numeric_limits<Scalar>::epsilon()));
  Scalar r1(s1.radius);
  Scalar r2(s2.radius);
  Scalar lz2(s2.halfLength);
  // boundaries of the cylinder axis
  Vec3s A(tf2.transform(Vec3s(0, 0, -lz2)));
  Vec3s B(tf2.transform(Vec3s(0, 0, lz2)));
  // Position of the center of the sphere
  Vec3s S(tf1.getTranslation());
  // axis of the cylinder
  Vec3s u(tf2.getRotation().col(2));
  /// @todo a tiny performance improvement could be achieved using the abscissa
  /// with S as the origin
  assert((B - A - (s2.halfLength * 2) * u).norm() < eps);
  Vec3s AS(S - A);
  // abscissa of S on cylinder axis with A as the origin
  Scalar s(u.dot(AS));
  Vec3s P(A + s * u);
  Vec3s PS(S - P);
  Scalar dPS = PS.norm();
  // Normal to cylinder axis such that plane (A, u, v) contains sphere
  // center
  Vec3s v(0, 0, 0);
  Scalar dist;
  if (dPS > eps) {
    // S is not on cylinder axis
    v = (1 / dPS) * PS;
  }
  if (s <= 0) {
    if (dPS <= r2) {
      // closest point on cylinder is on cylinder disc basis
      dist = -s - r1;
      p1 = S + r1 * u;
      p2 = A + dPS * v;
      normal = u;
    } else {
      // closest point on cylinder is on cylinder circle basis
      p2 = A + r2 * v;
      Vec3s Sp2(p2 - S);
      Scalar dSp2 = Sp2.norm();
      if (dSp2 > eps) {
        normal = (1 / dSp2) * Sp2;
        p1 = S + r1 * normal;
        dist = dSp2 - r1;
        assert(fabs(dist) - (p1 - p2).norm() < eps);
      } else {
        // Center of sphere is on cylinder boundary
        normal = p2 - .5 * (A + B);
        assert(u.dot(normal) >= 0);
        normal.normalize();
        dist = -r1;
        p1 = S + r1 * normal;
      }
    }
  } else if (s <= (s2.halfLength * 2)) {
    // 0 < s <= s2.lz
    normal = -v;
    dist = dPS - r1 - r2;
    p2 = P + r2 * v;
    p1 = S - r1 * v;
  } else {
    // lz < s
    if (dPS <= r2) {
      // closest point on cylinder is on cylinder disc basis
      dist = s - (s2.halfLength * 2) - r1;
      p1 = S - r1 * u;
      p2 = B + dPS * v;
      normal = -u;
    } else {
      // closest point on cylinder is on cylinder circle basis
      p2 = B + r2 * v;
      Vec3s Sp2(p2 - S);
      Scalar dSp2 = Sp2.norm();
      if (dSp2 > eps) {
        normal = (1 / dSp2) * Sp2;
        p1 = S + r1 * normal;
        dist = dSp2 - r1;
        assert(fabs(dist) - (p1 - p2).norm() < eps);
      } else {
        // Center of sphere is on cylinder boundary
        normal = p2 - .5 * (A + B);
        normal.normalize();
        p1 = S + r1 * normal;
        dist = -r1;
      }
    }
  }

  // Take swept-sphere radius into account
  const Scalar ssr1 = s1.getSweptSphereRadius();
  const Scalar ssr2 = s2.getSweptSphereRadius();
  if (ssr1 > 0 || ssr2 > 0) {
    p1 += ssr1 * normal;
    p2 -= ssr2 * normal;
    dist -= (ssr1 + ssr2);
  }

  return dist;
}

/// @param p1 witness point on the first Sphere.
/// @param p2 witness point on the second Sphere.
/// @param normal pointing from shape 1 to shape 2 (sphere1 to sphere2).
/// @return the distance between the two spheres (negative if penetration).
inline Scalar sphereSphereDistance(const Sphere& s1, const Transform3s& tf1,
                                   const Sphere& s2, const Transform3s& tf2,
                                   Vec3s& p1, Vec3s& p2, Vec3s& normal) {
  const coal::Vec3s& center1 = tf1.getTranslation();
  const coal::Vec3s& center2 = tf2.getTranslation();
  Scalar r1 = (s1.radius + s1.getSweptSphereRadius());
  Scalar r2 = (s2.radius + s2.getSweptSphereRadius());

  Vec3s c1c2 = center2 - center1;
  Scalar cdist = c1c2.norm();
  Vec3s unit(1, 0, 0);
  if (cdist > Eigen::NumTraits<Scalar>::epsilon()) unit = c1c2 / cdist;
  Scalar dist = cdist - r1 - r2;
  normal = unit;
  p1.noalias() = center1 + r1 * unit;
  p2.noalias() = center2 - r2 * unit;
  return dist;
}

/** @brief the minimum distance from a point to a line */
inline Scalar segmentSqrDistance(const Vec3s& from, const Vec3s& to,
                                 const Vec3s& p, Vec3s& nearest) {
  Vec3s diff = p - from;
  Vec3s v = to - from;
  Scalar t = v.dot(diff);

  if (t > 0) {
    Scalar dotVV = v.squaredNorm();
    if (t < dotVV) {
      t /= dotVV;
      diff -= v * t;
    } else {
      t = 1;
      diff -= v;
    }
  } else
    t = 0;

  nearest.noalias() = from + v * t;
  return diff.squaredNorm();
}

/// @brief Whether a point's projection is in a triangle
inline bool projectInTriangle(const Vec3s& p1, const Vec3s& p2, const Vec3s& p3,
                              const Vec3s& normal, const Vec3s& p) {
  Vec3s edge1(p2 - p1);
  Vec3s edge2(p3 - p2);
  Vec3s edge3(p1 - p3);

  Vec3s p1_to_p(p - p1);
  Vec3s p2_to_p(p - p2);
  Vec3s p3_to_p(p - p3);

  Vec3s edge1_normal(edge1.cross(normal));
  Vec3s edge2_normal(edge2.cross(normal));
  Vec3s edge3_normal(edge3.cross(normal));

  Scalar r1, r2, r3;
  r1 = edge1_normal.dot(p1_to_p);
  r2 = edge2_normal.dot(p2_to_p);
  r3 = edge3_normal.dot(p3_to_p);
  if ((r1 > 0 && r2 > 0 && r3 > 0) || (r1 <= 0 && r2 <= 0 && r3 <= 0)) {
    return true;
  }
  return false;
}

/// @param p1 witness point on the first Sphere.
/// @param p2 witness point on the second Sphere.
/// @param normal pointing from shape 1 to shape 2 (sphere1 to sphere2).
/// @return the distance between the two shapes (negative if penetration).
inline Scalar sphereTriangleDistance(const Sphere& s, const Transform3s& tf1,
                                     const TriangleP& tri,
                                     const Transform3s& tf2, Vec3s& p1,
                                     Vec3s& p2, Vec3s& normal) {
  const Vec3s& P1 = tf2.transform(tri.a);
  const Vec3s& P2 = tf2.transform(tri.b);
  const Vec3s& P3 = tf2.transform(tri.c);

  Vec3s tri_normal = (P2 - P1).cross(P3 - P1);
  tri_normal.normalize();
  const Vec3s& center = tf1.getTranslation();
  // Note: comparing an object with a swept-sphere radius of r1 against another
  // object with a swept-sphere radius of r2 is equivalent to comparing the
  // first object with a swept-sphere radius of r1 + r2 against the second
  // object with a swept-sphere radius of 0.
  const Scalar& radius =
      s.radius + s.getSweptSphereRadius() + tri.getSweptSphereRadius();
  assert(radius >= 0);
  assert(s.radius >= 0);
  Vec3s p1_to_center = center - P1;
  Scalar distance_from_plane = p1_to_center.dot(tri_normal);
  Vec3s closest_point(
      Vec3s::Constant(std::numeric_limits<Scalar>::quiet_NaN()));
  Scalar min_distance_sqr, distance_sqr;

  if (distance_from_plane < 0) {
    distance_from_plane *= -1;
    tri_normal *= -1;
  }

  if (projectInTriangle(P1, P2, P3, tri_normal, center)) {
    closest_point = center - tri_normal * distance_from_plane;
    min_distance_sqr = distance_from_plane * distance_from_plane;
  } else {
    // Compute distance to each edge and take minimal distance
    Vec3s nearest_on_edge;
    min_distance_sqr = segmentSqrDistance(P1, P2, center, closest_point);

    distance_sqr = segmentSqrDistance(P2, P3, center, nearest_on_edge);
    if (distance_sqr < min_distance_sqr) {
      min_distance_sqr = distance_sqr;
      closest_point = nearest_on_edge;
    }
    distance_sqr = segmentSqrDistance(P3, P1, center, nearest_on_edge);
    if (distance_sqr < min_distance_sqr) {
      min_distance_sqr = distance_sqr;
      closest_point = nearest_on_edge;
    }
  }

  normal = (closest_point - center).normalized();
  p1 = center + normal * (s.radius + s.getSweptSphereRadius());
  p2 = closest_point - normal * tri.getSweptSphereRadius();
  const Scalar distance = std::sqrt(min_distance_sqr) - radius;
  return distance;
}

/// @param p1 closest (or most penetrating) point on the Halfspace,
/// @param p2 closest (or most penetrating) point on the shape,
/// @param normal the halfspace normal.
/// @return the distance between the two shapes (negative if penetration).
inline Scalar halfspaceDistance(const Halfspace& h, const Transform3s& tf1,
                                const ShapeBase& s, const Transform3s& tf2,
                                Vec3s& p1, Vec3s& p2, Vec3s& normal) {
  // TODO(louis): handle multiple contact points when the halfspace normal is
  // parallel to the shape's surface (every primitive except sphere and
  // ellipsoid).

  // Express halfspace in world frame
  Halfspace new_h = transform(h, tf1);

  // Express halfspace normal in shape frame
  Vec3s n_2(tf2.getRotation().transpose() * new_h.n);

  // Compute support of shape in direction of halfspace normal
  int hint = 0;
  p2.noalias() =
      getSupport<details::SupportOptions::WithSweptSphere>(&s, -n_2, hint);
  p2 = tf2.transform(p2);

  const Scalar dist = new_h.signedDistance(p2);
  p1.noalias() = p2 - dist * new_h.n;
  normal.noalias() = new_h.n;

  const Scalar dummy_precision =
      std::sqrt(Eigen::NumTraits<Scalar>::dummy_precision());
  COAL_UNUSED_VARIABLE(dummy_precision);
  assert(new_h.distance(p1) <= dummy_precision);
  return dist;
}

/// @param p1 closest (or most penetrating) point on the Plane,
/// @param p2 closest (or most penetrating) point on the shape,
/// @param normal the halfspace normal.
/// @return the distance between the two shapes (negative if penetration).
inline Scalar planeDistance(const Plane& plane, const Transform3s& tf1,
                            const ShapeBase& s, const Transform3s& tf2,
                            Vec3s& p1, Vec3s& p2, Vec3s& normal) {
  // TODO(louis): handle multiple contact points when the plane normal is
  // parallel to the shape's surface (every primitive except sphere and
  // ellipsoid).

  // Express plane as two halfspaces in world frame
  std::array<Halfspace, 2> new_h = transformToHalfspaces(plane, tf1);

  // Express halfspace normals in shape frame
  Vec3s n_h1(tf2.getRotation().transpose() * new_h[0].n);
  Vec3s n_h2(tf2.getRotation().transpose() * new_h[1].n);

  // Compute support of shape in direction of halfspace normal and its opposite
  int hint = 0;
  Vec3s p2h1 =
      getSupport<details::SupportOptions::WithSweptSphere>(&s, -n_h1, hint);
  p2h1 = tf2.transform(p2h1);

  hint = 0;
  Vec3s p2h2 =
      getSupport<details::SupportOptions::WithSweptSphere>(&s, -n_h2, hint);
  p2h2 = tf2.transform(p2h2);

  Scalar dist1 = new_h[0].signedDistance(p2h1);
  Scalar dist2 = new_h[1].signedDistance(p2h2);

  const Scalar dummy_precision =
      std::sqrt(Eigen::NumTraits<Scalar>::dummy_precision());
  COAL_UNUSED_VARIABLE(dummy_precision);

  Scalar dist;
  if (dist1 >= dist2) {
    dist = dist1;
    p2.noalias() = p2h1;
    p1.noalias() = p2 - dist * new_h[0].n;
    normal.noalias() = new_h[0].n;
    assert(new_h[0].distance(p1) <= dummy_precision);
  } else {
    dist = dist2;
    p2.noalias() = p2h2;
    p1.noalias() = p2 - dist * new_h[1].n;
    normal.noalias() = new_h[1].n;
    assert(new_h[1].distance(p1) <= dummy_precision);
  }
  return dist;
}

/// Taken from book Real Time Collision Detection, from Christer Ericson
/// @param pb the witness point on the box surface
/// @param ps the witness point on the sphere.
/// @param normal pointing from box to sphere
/// @return the distance between the two shapes (negative if penetration).
inline Scalar boxSphereDistance(const Box& b, const Transform3s& tfb,
                                const Sphere& s, const Transform3s& tfs,
                                Vec3s& pb, Vec3s& ps, Vec3s& normal) {
  const Vec3s& os = tfs.getTranslation();
  const Vec3s& ob = tfb.getTranslation();
  const Matrix3s& Rb = tfb.getRotation();

  pb = ob;

  bool outside = false;
  const Vec3s os_in_b_frame(Rb.transpose() * (os - ob));
  int axis = -1;
  Scalar min_d = (std::numeric_limits<Scalar>::max)();
  for (int i = 0; i < 3; ++i) {
    Scalar facedist;
    if (os_in_b_frame(i) < -b.halfSide(i)) {  // outside
      pb.noalias() -= b.halfSide(i) * Rb.col(i);
      outside = true;
    } else if (os_in_b_frame(i) > b.halfSide(i)) {  // outside
      pb.noalias() += b.halfSide(i) * Rb.col(i);
      outside = true;
    } else {
      pb.noalias() += os_in_b_frame(i) * Rb.col(i);
      if (!outside &&
          (facedist = b.halfSide(i) - std::fabs(os_in_b_frame(i))) < min_d) {
        axis = i;
        min_d = facedist;
      }
    }
  }
  normal = pb - os;
  Scalar pdist = normal.norm();
  Scalar dist;    // distance between sphere and box
  if (outside) {  // pb is on the box
    dist = pdist - s.radius;
    normal /= -pdist;
  } else {  // pb is inside the box
    if (os_in_b_frame(axis) >= 0) {
      normal = Rb.col(axis);
    } else {
      normal = -Rb.col(axis);
    }
    dist = -min_d - s.radius;
  }
  ps = os - s.radius * normal;
  if (!outside || dist <= 0) {
    // project point pb onto the box's surface
    pb = ps - dist * normal;
  }

  // Take swept-sphere radius into account
  const Scalar ssrb = b.getSweptSphereRadius();
  const Scalar ssrs = s.getSweptSphereRadius();
  if (ssrb > 0 || ssrs > 0) {
    pb += ssrb * normal;
    ps -= ssrs * normal;
    dist -= (ssrb + ssrs);
  }

  return dist;
}

/// @brief return distance between two halfspaces
/// @param p1 the witness point on the first halfspace.
/// @param p2 the witness point on the second halfspace.
/// @param normal pointing from first to second halfspace.
/// @return the distance between the two shapes (negative if penetration).
///
/// @note If the two halfspaces don't have the same normal (or opposed
/// normals), they collide and their distance is set to -infinity as there is no
/// translation that can separate them; they have infinite penetration depth.
/// The points p1 and p2 are the same point and represent the origin of the
/// intersection line between the objects. The normal is the direction of this
/// line.
inline Scalar halfspaceHalfspaceDistance(const Halfspace& s1,
                                         const Transform3s& tf1,
                                         const Halfspace& s2,
                                         const Transform3s& tf2, Vec3s& p1,
                                         Vec3s& p2, Vec3s& normal) {
  Halfspace new_s1 = transform(s1, tf1);
  Halfspace new_s2 = transform(s2, tf2);

  Scalar distance;
  Vec3s dir = (new_s1.n).cross(new_s2.n);
  Scalar dir_sq_norm = dir.squaredNorm();

  if (dir_sq_norm < std::numeric_limits<Scalar>::epsilon())  // parallel
  {
    if (new_s1.n.dot(new_s2.n) > 0) {
      // If the two halfspaces have the same normal, one is inside the other
      // and they can't be separated. They have inifinte penetration depth.
      distance = -(std::numeric_limits<Scalar>::max)();
      if (new_s1.d <= new_s2.d) {
        normal = new_s1.n;
        p1 = normal * distance;
        p2 = new_s2.n * new_s2.d;
        assert(new_s2.distance(p2) <=
               Eigen::NumTraits<Scalar>::dummy_precision());
      } else {
        normal = -new_s1.n;
        p1 << new_s1.n * new_s1.d;
        p2 = -(normal * distance);
        assert(new_s1.distance(p1) <=
               Eigen::NumTraits<Scalar>::dummy_precision());
      }
    } else {
      distance = -(new_s1.d + new_s2.d);
      normal = new_s1.n;
      p1 = new_s1.n * new_s1.d;
      p2 = new_s2.n * new_s2.d;
    }
  } else {
    // If the halfspaces are not parallel, they are in collision.
    // Their distance, in the sens of the norm of separation vector, is infinite
    // (it's impossible to find a translation which separates them)
    distance = -(std::numeric_limits<Scalar>::max)();
    // p1 and p2 are the same point, corresponding to a point on the
    // intersection line between the two objects. Normal is the direction of
    // that line.
    normal = dir;
    p1 = p2 =
        ((new_s2.n * new_s1.d - new_s1.n * new_s2.d).cross(dir)) / dir_sq_norm;
    // Sources: https://en.wikipedia.org/wiki/Plane%E2%80%93plane_intersection
    // and      https://en.wikipedia.org/wiki/Cross_product
  }

  // Take swept-sphere radius into account
  const Scalar ssr1 = s1.getSweptSphereRadius();
  const Scalar ssr2 = s2.getSweptSphereRadius();
  if (ssr1 > 0 || ssr2 > 0) {
    p1 += ssr1 * normal;
    p2 -= ssr2 * normal;
    distance -= (ssr1 + ssr2);
  }

  return distance;
}

/// @brief return distance between plane and halfspace.
/// @param p1 the witness point on the halfspace.
/// @param p2 the witness point on the plane.
/// @param normal pointing from halfspace to plane.
/// @return the distance between the two shapes (negative if penetration).
///
/// @note If plane and halfspace don't have the same normal (or opposed
/// normals), they collide and their distance is set to -infinity as there is no
/// translation that can separate them; they have infinite penetration depth.
/// The points p1 and p2 are the same point and represent the origin of the
/// intersection line between the objects. The normal is the direction of this
/// line.
inline Scalar halfspacePlaneDistance(const Halfspace& s1,
                                     const Transform3s& tf1, const Plane& s2,
                                     const Transform3s& tf2, Vec3s& p1,
                                     Vec3s& p2, Vec3s& normal) {
  Halfspace new_s1 = transform(s1, tf1);
  Plane new_s2 = transform(s2, tf2);

  Scalar distance;
  Vec3s dir = (new_s1.n).cross(new_s2.n);
  Scalar dir_sq_norm = dir.squaredNorm();

  if (dir_sq_norm < std::numeric_limits<Scalar>::epsilon())  // parallel
  {
    normal = new_s1.n;
    distance = new_s1.n.dot(new_s2.n) > 0 ? (new_s2.d - new_s1.d)
                                          : -(new_s1.d + new_s2.d);
    p1 = new_s1.n * new_s1.d;
    p2 = new_s2.n * new_s2.d;
    assert(new_s1.distance(p1) <= Eigen::NumTraits<Scalar>::dummy_precision());
    assert(new_s2.distance(p2) <= Eigen::NumTraits<Scalar>::dummy_precision());
  } else {
    // If the halfspace and plane are not parallel, they are in collision.
    // Their distance, in the sens of the norm of separation vector, is infinite
    // (it's impossible to find a translation which separates them)
    distance = -(std::numeric_limits<Scalar>::max)();
    // p1 and p2 are the same point, corresponding to a point on the
    // intersection line between the two objects. Normal is the direction of
    // that line.
    normal = dir;
    p1 = p2 =
        ((new_s2.n * new_s1.d - new_s1.n * new_s2.d).cross(dir)) / dir_sq_norm;
    // Sources: https://en.wikipedia.org/wiki/Plane%E2%80%93plane_intersection
    // and      https://en.wikipedia.org/wiki/Cross_product
  }

  // Take swept-sphere radius into account
  const Scalar ssr1 = s1.getSweptSphereRadius();
  const Scalar ssr2 = s2.getSweptSphereRadius();
  if (ssr1 > 0 || ssr2 > 0) {
    p1 += ssr1 * normal;
    p2 -= ssr2 * normal;
    distance -= (ssr1 + ssr2);
  }

  return distance;
}

/// @brief return distance between two planes
/// @param p1 the witness point on the first plane.
/// @param p2 the witness point on the second plane.
/// @param normal pointing from first to second plane.
/// @return the distance between the two shapes (negative if penetration).
///
/// @note If the two planes don't have the same normal (or opposed
/// normals), they collide and their distance is set to -infinity as there is no
/// translation that can separate them; they have infinite penetration depth.
/// The points p1 and p2 are the same point and represent the origin of the
/// intersection line between the objects. The normal is the direction of this
/// line.
inline Scalar planePlaneDistance(const Plane& s1, const Transform3s& tf1,
                                 const Plane& s2, const Transform3s& tf2,
                                 Vec3s& p1, Vec3s& p2, Vec3s& normal) {
  Plane new_s1 = transform(s1, tf1);
  Plane new_s2 = transform(s2, tf2);

  Scalar distance;
  Vec3s dir = (new_s1.n).cross(new_s2.n);
  Scalar dir_sq_norm = dir.squaredNorm();

  if (dir_sq_norm < std::numeric_limits<Scalar>::epsilon())  // parallel
  {
    p1 = new_s1.n * new_s1.d;
    p2 = new_s2.n * new_s2.d;
    assert(new_s1.distance(p1) <= Eigen::NumTraits<Scalar>::dummy_precision());
    assert(new_s2.distance(p2) <= Eigen::NumTraits<Scalar>::dummy_precision());
    distance = (p1 - p2).norm();

    if (distance > Eigen::NumTraits<Scalar>::dummy_precision()) {
      normal = (p2 - p1).normalized();
    } else {
      normal = new_s1.n;
    }
  } else {
    // If the planes are not parallel, they are in collision.
    // Their distance, in the sens of the norm of separation vector, is infinite
    // (it's impossible to find a translation which separates them)
    distance = -(std::numeric_limits<Scalar>::max)();
    // p1 and p2 are the same point, corresponding to a point on the
    // intersection line between the two objects. Normal is the direction of
    // that line.
    normal = dir;
    p1 = p2 =
        ((new_s2.n * new_s1.d - new_s1.n * new_s2.d).cross(dir)) / dir_sq_norm;
    // Sources: https://en.wikipedia.org/wiki/Plane%E2%80%93plane_intersection
    // and      https://en.wikipedia.org/wiki/Cross_product
  }

  // Take swept-sphere radius into account
  const Scalar ssr1 = s1.getSweptSphereRadius();
  const Scalar ssr2 = s2.getSweptSphereRadius();
  if (ssr1 > 0 || ssr2 > 0) {
    p1 += ssr1 * normal;
    p2 -= ssr2 * normal;
    distance -= (ssr1 + ssr2);
  }

  return distance;
}

/// See the prototype below
inline Scalar computePenetration(const Vec3s& P1, const Vec3s& P2,
                                 const Vec3s& P3, const Vec3s& Q1,
                                 const Vec3s& Q2, const Vec3s& Q3,
                                 Vec3s& normal) {
  Vec3s u((P2 - P1).cross(P3 - P1));
  normal = u.normalized();
  Scalar depth1((P1 - Q1).dot(normal));
  Scalar depth2((P1 - Q2).dot(normal));
  Scalar depth3((P1 - Q3).dot(normal));
  return std::max(depth1, std::max(depth2, depth3));
}

// Compute penetration distance and normal of two triangles in collision
// Normal is normal of triangle 1 (P1, P2, P3), penetration depth is the
// minimal distance (Q1, Q2, Q3) should be translated along the normal so
// that the triangles are collision free.
//
// Note that we compute here an upper bound of the penetration distance,
// not the exact value.
inline Scalar computePenetration(const Vec3s& P1, const Vec3s& P2,
                                 const Vec3s& P3, const Vec3s& Q1,
                                 const Vec3s& Q2, const Vec3s& Q3,
                                 const Transform3s& tf1, const Transform3s& tf2,
                                 Vec3s& normal) {
  Vec3s globalP1(tf1.transform(P1));
  Vec3s globalP2(tf1.transform(P2));
  Vec3s globalP3(tf1.transform(P3));
  Vec3s globalQ1(tf2.transform(Q1));
  Vec3s globalQ2(tf2.transform(Q2));
  Vec3s globalQ3(tf2.transform(Q3));
  return computePenetration(globalP1, globalP2, globalP3, globalQ1, globalQ2,
                            globalQ3, normal);
}

// inline Scalar capsuleBox(Contact* con, double margin, const Eigen::Vector3d&
// pos1,
//                const Eigen::Matrix3d& mat1, const Eigen::Vector3d& size1,
//                const Eigen::Vector3d& pos2, const Eigen::Matrix3d& mat2,
//                const Eigen::Vector3d& size2) {

inline Scalar capsuleBoxDistance(const Capsule& capsule, const Transform3s& tf1,
                                 const Box& box, const Transform3s& tf2,
                                 Vec3s& p1, Vec3s& p2, Vec3s& normal) {
  // Temporary variables
  Vec3s tmp1, tmp2, tmp3, halfaxis, axis, dif;
  Vec3s pos;  // position of capsule in box-local frame

  Scalar halflength;      // half of capsule's length
  Scalar bestdist;        // closest contact point distance
  Scalar bestdistmax;     // init value for bestdist
  Scalar bestsegmentpos;  // between -1 and 1 : which point on the segment is
                          // closest to the box
  Scalar secondpos;  // distance of 2nd contact position on capsule segment from
                     // the first
  Scalar dist;
  Scalar bestboxpos;  // closest contact point, position on the box's edge
  Scalar mul, e1, e2, dp, de;

  Scalar ma, mb, mc, u, v, det, x1, x2, idet;  // linelinedist temps

  int s1, s2;        // hold linelinedist info
  int i, j, c1, c2;  // temporary variables
  int cltype = -4;   // closest type
  int clface;        // closest face
  int clcorner = 0;  // closest corner (0..7 in binary)
  int cledge;        // closest edge axis
  int axisdir;       // direction of capsule axis in relation to the box
  int ax1, ax2, ax;  // axis temporaries

  halflength = capsule.halfLength;
  secondpos =
      -4;  // initialize to no 2nd contact (valid values are between -1 and 1)

  // Bring capsule to box-local frame (center's box is at (0,0,0))
  tmp1.noalias() = tf1.translation() - tf2.translation();
  // And axis parallel to world
  pos.noalias() = tf2.rotation().transpose() * tmp1;

  // Capsule's axis
  tmp1 = tf1.rotation().col(2);

  // Do the same for the capsule axis
  axis.noalias() = tf2.rotation().transpose() * tmp1;
  // Scale to get actual capsule half-axis
  halfaxis = axis * halflength;

  axisdir = 0;
  if (halfaxis(0) > 0) axisdir += 1;
  if (halfaxis(1) > 0) axisdir += 2;
  if (halfaxis(2) > 0) axisdir += 4;

  // Under this notion "axisdir" and "7-axisdir" point in opposite directions,
  // essentially the same for a capsule

  bestdistmax = 2 * (capsule.radius + halflength + box.halfSide(0) +
                     box.halfSide(1) + box.halfSide(2));  // initialize bestdist
  bestdist = bestdistmax;
  bestsegmentpos = 0;

  tmp2.setZero();

  // Test to see if maybe a face of the box is closest to the capsule
  for (i = -1; i <= 1; i += 2) {
    tmp1 = pos + halfaxis * i;
    tmp2 = tmp1;

    for (c1 = 0, j = 0, c2 = -1; j < 3; j++) {
      if (tmp1(j) < -box.halfSide(j)) {
        c1++;
        c2 = j;
        tmp1(j) = -box.halfSide(j);
      } else if (tmp1(j) > box.halfSide(j)) {
        c1++;
        c2 = j;
        tmp1(j) = box.halfSide(j);
      }
    }

    if (c1 > 1) continue;

    tmp1 = tmp1 - tmp2;
    dist = tmp1.dot(tmp1);

    if (dist < bestdist) {
      bestdist = dist;
      bestsegmentpos = i;
      cltype = -2 + i;
      clface = c2;
    }
  }

  tmp2.setZero();

  // Check for corners and edges
  for (j = 0; j < 3; j++) {
    for (i = 0; i < 8; i++) {
      if ((i & (1 << j)) == 0) {
        // Trick to get a corner
        tmp3(0) = ((i & 1) ? 1 : -1) * box.halfSide(0);
        tmp3(1) = ((i & 2) ? 1 : -1) * box.halfSide(1);
        tmp3(2) = ((i & 4) ? 1 : -1) * box.halfSide(2);
        tmp3(j) = 0;

        // tmp3 is the starting point on the box
        // tmp2 is the direction along the "j"-th axis
        // pos is the capsule's center
        // halfaxis is the capsule direction

        // Find closest point between capsule and the edge
        dif = tmp3 - pos;

        ma = box.halfSide(j) * box.halfSide(j);
        mb = -box.halfSide(j) * halfaxis(j);
        mc = capsule.halfLength * capsule.halfLength;

        u = -box.halfSide(j) * dif(j);
        v = halfaxis.dot(dif);

        det = ma * mc - mb * mb;
        if (std::abs(det) < Eigen::NumTraits<Scalar>::dummy_precision())
          continue;
        idet = 1 / det;

        // sX : X=1 means middle of segment. X=0 or 2 one or the other end
        x1 = (mc * u - mb * v) * idet;
        x2 = (ma * v - mb * u) * idet;

        s1 = s2 = 1;

        if (x1 > 1) {
          x1 = 1;
          s1 = 2;
          x2 = (v - mb) * (1 / mc);
        } else if (x1 < -1) {
          x1 = -1;
          s1 = 0;
          x2 = (v + mb) * (1 / mc);
        }

        if (x2 > 1) {
          x2 = 1;
          s2 = 2;
          x1 = (u - mb) * (1 / ma);
          if (x1 > 1)
            x1 = 1, s1 = 2;
          else if (x1 < -1)
            x1 = -1, s1 = 0;
        } else if (x2 < -1) {
          x2 = -1;
          s2 = 0;
          x1 = (u + mb) * (1 / ma);
          if (x1 > 1)
            x1 = 1, s1 = 2;
          else if (x1 < -1)
            x1 = -1, s1 = 0;
        }

        dif = tmp3 - pos;
        dif += halfaxis * (-x2);
        dif(j) += box.halfSide(j) * x1;

        tmp1(2) = dif.dot(dif);

        c1 = s1 * 3 + s2;

        // The -MINVAL might not be necessary. Fixes numerical problem when axis
        // is numerically parallel to the box
        if (tmp1(2) < bestdist - Eigen::NumTraits<Scalar>::dummy_precision()) {
          bestdist = tmp1(2);
          bestsegmentpos = x2;
          bestboxpos = x1;

          // c1<6 means that closest point on the box is at the lower end
          // or in the middle of the edge
          c2 = c1 / 6;

          clcorner = i + (1 << j) * c2;  // which corner is the closest
          cledge = j;                    // which axis
          cltype = c1;                   // save clamped info
        }
      }
    }
  }

  // Special case for 2D plane checks - this is complicated math for specific
  // case handling I've retained the structure but simplified the notation where
  // possible
  for (j = 0; j < 3; j++) {
    if (j == 2) {
      struct Vec2D {
        Scalar x, y;

        Vec2D operator-(const Vec2D& other) const {
          return {x - other.x, y - other.y};
        }

        Scalar dot(const Vec2D& other) const {
          return x * other.x + y * other.y;
        }
      };

      Vec2D p, s, dd;
      Scalar uu, vv, w, ee1, best, l;

      bestdist = bestdistmax;

      p = {pos(0), pos(1)};
      dd = {halfaxis(0), halfaxis(1)};
      s = {box.halfSide(0), box.halfSide(1)};

      l = std::sqrt(dd.x * dd.x + dd.y * dd.y);

      uu = dd.x * s.y;
      vv = dd.y * s.x;
      w = dd.x * p.y - dd.y * p.x;

      best = -1;

      ee1 = +uu - vv;
      if ((ee1 < 0) == (w < 0)) {
        if (best < std::abs(ee1)) {
          best = std::abs(ee1);
          c1 = 0;
        }
      }
      ee1 = -uu - vv;
      if ((ee1 < 0) == (w < 0)) {
        if (best < std::abs(ee1)) {
          best = std::abs(ee1);
          c1 = 1;
        }
      }
      ee1 = +uu + vv;
      if ((ee1 < 0) == (w < 0)) {
        if (best < std::abs(ee1)) {
          best = std::abs(ee1);
          c1 = 2;
        }
      }
      ee1 = -uu + vv;
      if ((ee1 < 0) == (w < 0)) {
        if (best < std::abs(ee1)) {
          best = std::abs(ee1);
          c1 = 3;
        }
      }

      ee1 = std::abs(w) / l;
      ee1 = dd.x * dd.x + dd.y * dd.y;
      ee1 = p.x + (+s.y - p.y) / dd.y * dd.x;
    }
  }

  // Invalid type
  if (cltype == -4) return 0;

  if (cltype >= 0 && cltype / 3 != 1) {  // closest to a corner of the box
    c1 = axisdir ^ clcorner;

    // Hack to find the relative orientation of capsule and corner
    // there are 2 cases:
    //    1: pointing to or away from the corner
    //    2: oriented along a face or an edge

    if (c1 == 0 || c1 == 7)
      goto skip;  // case 1: no chance of additional contact

    if (c1 == 1 || c1 == 2 || c1 == 4) {
      mul = 1;
      de = 1 - bestsegmentpos;
      dp = 1 + bestsegmentpos;
    }

    if (c1 == 3 || c1 == 5 || c1 == 6) {
      mul = -1;
      c1 = 7 - c1;
      dp = 1 - bestsegmentpos;
      de = 1 + bestsegmentpos;
    }

    // "de" and "dp" distance from first closest point on the capsule to both
    // ends of it mul is a direction along the capsule's axis

    if (c1 == 1) ax = 0, ax1 = 1, ax2 = 2;
    if (c1 == 2) ax = 1, ax1 = 2, ax2 = 0;
    if (c1 == 4) ax = 2, ax1 = 0, ax2 = 1;

    if (axis(ax) * axis(ax) > 0.5) {  // second point along the edge of the box
      secondpos = de;                 // initial position from the
      e1 = 2 * box.halfSide(ax) / std::abs(halfaxis(ax));

      if (e1 < secondpos) {
        secondpos =
            e1;  // we overshoot, move back to the other corner of the edge
      }
      secondpos *= mul;
    } else {  // second point along a face of the box
      secondpos = dp;

      // check for overshoot again
      e1 = 2 * box.halfSide(ax1) / std::abs(halfaxis(ax1));
      if (e1 < secondpos) secondpos = e1;

      e1 = 2 * box.halfSide(ax2) / std::abs(halfaxis(ax2));
      if (e1 < secondpos) secondpos = e1;

      secondpos *= -mul;
    }
  } else if (cltype >= 0 && cltype / 3 == 1) {  // we are on box's edge
    // Hacks to find the relative orientation of capsule and edge
    // there are 2 cases:
    //    c1= 2^n: edge and capsule are oriented in a T configuration (no more
    //    contacts) c1!=2^n: oriented in a cross X configuration

    c1 = axisdir ^ clcorner;  // same trick
    c1 &= 7 - (1 << cledge);  // even more hacks

    if (c1 != 1 && c1 != 2 && c1 != 4) goto skip;

    if (cledge == 0) ax1 = 1, ax2 = 2;
    if (cledge == 1) ax1 = 2, ax2 = 0;
    if (cledge == 2) ax1 = 0, ax2 = 1;
    ax = cledge;

    // Then it finds with which face the capsule has a lower angle and switches
    // the axis names
    if (std::abs(axis(ax1)) > std::abs(axis(ax2))) ax1 = ax2;
    ax2 = 3 - ax - ax1;

    // Keep track of the axis orientation (mul will tell us which direction
    // along the capsule to find the second point) you can notice all other
    // references to the axis "halfaxis" are with absolute value
    if (c1 & (1 << ax2)) {
      mul = 1;
      secondpos = 1 - bestsegmentpos;
    } else {
      mul = -1;
      secondpos = 1 + bestsegmentpos;
    }

    // Now we have to find out whether we point towards the opposite side or
    // towards one of the sides and also find the farthest point along the
    // capsule that is above the box
    e1 = 2 * box.halfSide(ax2) / std::abs(halfaxis(ax2));
    if (e1 < secondpos) secondpos = e1;

    if (((axisdir & (1 << ax)) != 0) ==
        ((c1 & (1 << ax2)) != 0))  // that is insane
      e2 = 1 - bestboxpos;
    else
      e2 = 1 + bestboxpos;

    e1 = box.halfSide(ax) * e2 / std::abs(halfaxis(ax));

    if (e1 < secondpos) secondpos = e1;

    secondpos *= mul;
  } else if (cltype < 0) {
    // Similarly we handle the case when one capsule's end is closest to a face
    // of the box and find where is the other end pointing to and clamping to
    // the farthest point of the capsule that's above the box
    if (clface == -1)
      goto skip;  // here the closest point is inside the box, no need for a
                  // second point
    if (cltype == -3)
      mul = 1;
    else
      mul = -1;

    secondpos = 2;

    tmp1 = pos;
    tmp1 += halfaxis * (-mul);

    for (i = 0; i < 3; i++) {
      if (i != clface) {
        e1 = (box.halfSide(i) - tmp1(i)) / halfaxis(i) * mul;
        if (e1 > 0)
          if (e1 < secondpos) secondpos = e1;

        e1 = (-box.halfSide(i) - tmp1(i)) / halfaxis(i) * mul;
        if (e1 > 0)
          if (e1 < secondpos) secondpos = e1;
      }
    }
    secondpos *= mul;
  }

skip:
  // Create sphere in original orientation at first contact point
  tmp1 = pos + halfaxis * bestsegmentpos;
  tmp2.noalias() = tf2.rotation() * tmp1 + tf2.translation();

  Sphere sphere(tmp2(0));
  const Scalar res = boxSphereDistance(box, tf2, sphere, tf1, p2, p1, normal);
  normal *= -1;

  return res;
}

}  // namespace details
}  // namespace coal

#endif  // COAL_SRC_NARROWPHASE_DETAILS_H
