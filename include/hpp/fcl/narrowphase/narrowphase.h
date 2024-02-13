/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2011-2014, Willow Garage, Inc.
 *  Copyright (c) 2014-2015, Open Source Robotics Foundation
 *  Copyright (c) 2018-2019, Centre National de la Recherche Scientifique
 *  All rights reserved.
 *  Copyright (c) 2021-2023, INRIA
 *  All rights reserved.
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

#ifndef HPP_FCL_NARROWPHASE_H
#define HPP_FCL_NARROWPHASE_H

#include <limits>
#include <iostream>

#include <hpp/fcl/narrowphase/gjk.h>
#include <hpp/fcl/collision_data.h>
#include <hpp/fcl/narrowphase/narrowphase_defaults.h>

namespace hpp {
namespace fcl {

/// @brief collision and distance solver based on GJK algorithm implemented in
/// fcl (rewritten the code from the GJK in bullet)
struct HPP_FCL_DLLAPI GJKSolver {
  typedef Eigen::Array<FCL_REAL, 1, 2> Array2d;

  /// @brief initialize GJK
  template <typename S1, typename S2>
  void getGJKInitialGuess(const details::MinkowskiDiff& shape, const S1& s1,
                          const S2& s2, Vec3f& guess,
                          support_func_guess_t& support_hint) const {
    switch (gjk_initial_guess) {
      case GJKInitialGuess::DefaultGuess:
        guess = Vec3f(1, 0, 0);
        support_hint.setZero();
        break;
      case GJKInitialGuess::CachedGuess:
        guess = cached_guess;
        support_hint = support_func_cached_guess;
        break;
      case GJKInitialGuess::BoundingVolumeGuess:
        if (s1.aabb_local.volume() < 0 || s2.aabb_local.volume() < 0) {
          HPP_FCL_THROW_PRETTY(
              "computeLocalAABB must have been called on the shapes before "
              "using "
              "GJKInitialGuess::BoundingVolumeGuess.",
              std::logic_error);
        }
        guess.noalias() = s1.aabb_local.center() -
                          (shape.oR1 * s2.aabb_local.center() + shape.ot1);
        support_hint.setZero();
        break;
      default:
        HPP_FCL_THROW_PRETTY("Wrong initial guess for GJK.", std::logic_error);
    }
    // TODO: use gjk_initial_guess instead
    HPP_FCL_COMPILER_DIAGNOSTIC_PUSH
    HPP_FCL_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS
    if (enable_cached_guess) {
      guess = cached_guess;
      support_hint = support_func_cached_guess;
    }
    HPP_FCL_COMPILER_DIAGNOSTIC_POP
  }

  /// @brief intersection checking between two shapes.
  /// @return true if the shapes are colliding.
  /// @note the variables `gjk_status` and `epa_status` can be used to
  /// check if GJK or EPA ran successfully.
  /// @note (12 feb, 2024): this function is never called in the
  /// hpp-fcl library.
  /// Either the `shapeDistance` in this file is called or a specialized
  /// function. However, this function is still tested in the test suite,
  /// notably in `test/collision.cpp` and `test/geometric_shapes.cpp`.
  template <typename S1, typename S2>
  bool shapeIntersect(const S1& s1, const Transform3f& tf1, const S2& s2,
                      const Transform3f& tf2, FCL_REAL& distance_lower_bound,
                      bool enable_penetration, Vec3f* contact_points,
                      Vec3f* normal) const {
    details::MinkowskiDiff shape;
    shape.set(&s1, &s2, tf1, tf2);

    // Reset GJK algorithm
    gjk.reset();

    Vec3f guess;
    support_func_guess_t support_hint;
    getGJKInitialGuess(shape, s1, s2, guess, support_hint);

    gjk.evaluate(shape, guess, support_hint);
    HPP_FCL_COMPILER_DIAGNOSTIC_PUSH
    HPP_FCL_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS
    if (gjk_initial_guess == GJKInitialGuess::CachedGuess ||
        enable_cached_guess) {
      cached_guess = gjk.getGuessFromSimplex();
      support_func_cached_guess = gjk.support_hint;
    }
    HPP_FCL_COMPILER_DIAGNOSTIC_POP

    Vec3f w0, w1;
    switch (gjk.status) {
      case details::GJK::Inside:
        if (!enable_penetration && contact_points == NULL && normal == NULL)
          return true;
        if (gjk.hasPenetrationInformation(shape)) {
          gjk.getClosestPoints(shape, w0, w1);
          distance_lower_bound = gjk.distance;
          if (normal)
            (*normal).noalias() = tf1.getRotation() * (w1 - w0).normalized();
          if (contact_points) *contact_points = tf1.transform((w0 + w1) / 2);
          return true;
        } else {
          // Reset EPA algorithm
          epa.reset(epa_max_face_num, epa_max_vertex_num);

          epa.evaluate(gjk, -guess);
          if (epa.status & details::EPA::Valid ||
              epa.status == details::EPA::OutOfFaces        // Warnings
              || epa.status == details::EPA::OutOfVertices  // Warnings
          ) {
            epa.getClosestPoints(shape, w0, w1);
            distance_lower_bound = -epa.depth;
            if (normal) (*normal).noalias() = tf1.getRotation() * epa.normal;
            if (contact_points)
              *contact_points =
                  tf1.transform(w0 - epa.normal * (epa.depth * 0.5));
            return true;
          } else if (epa.status == details::EPA::FallBack) {
            epa.getClosestPoints(shape, w0, w1);
            distance_lower_bound = -epa.depth;  // Should be zero
            if (normal) (*normal).noalias() = tf1.getRotation() * epa.normal;
            if (contact_points) *contact_points = tf1.transform(w0);
            return true;
          }
          distance_lower_bound = -(std::numeric_limits<FCL_REAL>::max)();
          // EPA failed but we know there is a collision so we should
          return true;
        }
        break;
      case details::GJK::Valid:
        distance_lower_bound = gjk.distance;
        break;
      default:;
    }

    return false;
  }

  //// @brief intersection checking between one shape and a triangle with
  /// transformation
  /// @return true if the shape are colliding.
  /// The variables `gjk_status` and `epa_status` can be used to
  /// check if GJK or EPA ran successfully.
  template <typename S>
  bool shapeTriangleInteraction(const S& s, const Transform3f& tf1,
                                const Vec3f& P1, const Vec3f& P2,
                                const Vec3f& P3, const Transform3f& tf2,
                                FCL_REAL& distance, Vec3f& p1, Vec3f& p2,
                                Vec3f& normal) const {
    bool col = false;
    // Express everything in frame 1
    const Transform3f tf_1M2(tf1.inverseTimes(tf2));
    TriangleP tri(tf_1M2.transform(P1), tf_1M2.transform(P2),
                  tf_1M2.transform(P3));

    details::MinkowskiDiff shape;
    shape.set(&s, &tri);

    // Reset GJK algorithm
    gjk.reset();

    Vec3f guess;
    support_func_guess_t support_hint;
    getGJKInitialGuess(shape, s, tri, guess, support_hint);

    gjk.evaluate(shape, guess, support_hint);

    HPP_FCL_COMPILER_DIAGNOSTIC_PUSH
    HPP_FCL_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS
    if (gjk_initial_guess == GJKInitialGuess::CachedGuess ||
        enable_cached_guess) {
      cached_guess = gjk.getGuessFromSimplex();
      support_func_cached_guess = gjk.support_hint;
    }
    HPP_FCL_COMPILER_DIAGNOSTIC_PUSH

    Vec3f w0, w1;
    switch (gjk.status) {
      case details::GJK::Inside:
        col = true;
        if (gjk.hasPenetrationInformation(shape)) {
          gjk.getClosestPoints(shape, w0, w1);
          distance = gjk.distance;
          normal.noalias() = tf1.getRotation() * (w0 - w1).normalized();
          p1 = tf1.transform(w0);
          p2 = tf1.transform(w1);
        } else {
          // Reset EPA algorithm
          epa.reset(epa_max_face_num, epa_max_vertex_num);

          epa.evaluate(gjk, -guess);
          if (epa.status & details::EPA::Valid ||
              epa.status == details::EPA::OutOfFaces        // Warnings
              || epa.status == details::EPA::OutOfVertices  // Warnings
          ) {
            epa.getClosestPoints(shape, w0, w1);
            distance = -epa.depth;
            normal.noalias() = tf1.getRotation() * epa.normal;
            p1 = tf1.transform(w0);
            p2 = tf1.transform(w1);
            HPP_FCL_ASSERT(
                distance <= epa.tolerance,
                "Distance should be negative (at least below EPA tolerance).",
                std::logic_error);
          } else {
            distance = -(std::numeric_limits<FCL_REAL>::max)();
            gjk.getClosestPoints(shape, w0, w1);
            p1 = tf1.transform(w0);
            p2 = tf1.transform(w1);
          }
        }
        break;
      case details::GJK::Valid:
      case details::GJK::EarlyStopped:
        col = false;
        gjk.getClosestPoints(shape, w0, w1);
        distance = gjk.distance;
        normal.noalias() = -tf1.getRotation() * gjk.ray;
        normal.normalize();
        p1 = tf1.transform(w0);
        p2 = tf1.transform(w1);
        break;
      case details::GJK::Failed:
        col = false;

        gjk.getClosestPoints(shape, p1, p2);
        // TODO On degenerated case, the closest point may be wrong
        // (i.e. an object face normal is colinear to gjk.ray
        // assert (distance == (w0 - w1).norm());
        distance = gjk.distance;

        p1 = tf1.transform(p1);
        p2 = tf1.transform(p2);
        normal.setZero();
        assert(distance > 0);
        break;
      default:
        HPP_FCL_ASSERT(false, "GJKSolver: should not reach this part.",
                       std::logic_error)
    }
    return col;
  }

  /// @brief distance computation between two shapes.
  /// @return true if no error occured, false otherwise.
  /// The variables `gjk_status` and `epa_status` can be used to
  /// understand the reason of the failure.
  template <typename S1, typename S2>
  bool shapeDistance(const S1& s1, const Transform3f& tf1, const S2& s2,
                     const Transform3f& tf2, FCL_REAL& distance, Vec3f& p1,
                     Vec3f& p2, Vec3f& normal) const {
#ifndef NDEBUG
    FCL_REAL eps(sqrt(std::numeric_limits<FCL_REAL>::epsilon()));
#endif
    details::MinkowskiDiff shape;
    shape.set(&s1, &s2, tf1, tf2);

    // Reset GJK algorithm
    gjk.reset();

    Vec3f guess;
    support_func_guess_t support_hint;
    getGJKInitialGuess(shape, s1, s2, guess, support_hint);

    gjk.evaluate(shape, guess, support_hint);
    if (gjk_initial_guess == GJKInitialGuess::CachedGuess ||
        enable_cached_guess) {
      cached_guess = gjk.getGuessFromSimplex();
      support_func_cached_guess = gjk.support_hint;
    }

    if (gjk.status == details::GJK::Failed) {
      // TODO: understand why GJK fails between cylinder and box
      assert(distance * distance < sqrt(eps));
      Vec3f w0, w1;
      gjk.getClosestPoints(shape, w0, w1);
      distance = 0;
      p1 = tf1.transform(w0);
      p2 = tf1.transform(w1);
      normal.setZero();
      return false;
    } else if (gjk.status == details::GJK::Valid) {
      gjk.getClosestPoints(shape, p1, p2);
      // TODO On degenerated case, the closest point may be wrong
      // (i.e. an object face normal is colinear to gjk.ray
      // assert (distance == (w0 - w1).norm());
      distance = gjk.distance;

      normal.noalias() = -tf1.getRotation() * gjk.ray;
      normal.normalize();
      p1 = tf1.transform(p1);
      p2 = tf1.transform(p2);
      return true;
    } else if (gjk.status == details::GJK::EarlyStopped) {
      distance = gjk.distance;
      p1 = p2 = normal =
          Vec3f::Constant(std::numeric_limits<FCL_REAL>::quiet_NaN());
      return true;
    } else {
      assert(gjk.status == details::GJK::Inside);
      if (gjk.hasPenetrationInformation(shape)) {
        gjk.getClosestPoints(shape, p1, p2);
        distance = gjk.distance;
        // Return contact points in case of collision
        normal.noalias() = tf1.getRotation() * (p1 - p2);
        normal.normalize();
        p1 = tf1.transform(p1);
        p2 = tf1.transform(p2);
        return true;
      } else {
        // Reset EPA algorithm
        epa.reset(epa_max_face_num, epa_max_vertex_num);

        epa.evaluate(gjk, -guess);
        if (epa.status & details::EPA::Valid ||
            epa.status == details::EPA::OutOfFaces        // Warnings
            || epa.status == details::EPA::OutOfVertices  // Warnings
            || epa.status == details::EPA::FallBack) {
          Vec3f w0, w1;
          epa.getClosestPoints(shape, w0, w1);
          distance = (std::min)(0., -epa.depth);
          HPP_FCL_ASSERT(distance <= epa.tolerance,
                         "Distance should be negative (or at least below EPA's "
                         "tolerance).",
                         std::logic_error);
          normal.noalias() = tf1.getRotation() * epa.normal;
          p1 = tf1.transform(w0);
          p2 = tf1.transform(w1);
          return true;
        }
        distance = -(std::numeric_limits<FCL_REAL>::max)();
        gjk.getClosestPoints(shape, p1, p2);
        p1 = tf1.transform(p1);
        p2 = tf1.transform(p2);
      }
      return false;
    }
  }

  HPP_FCL_COMPILER_DIAGNOSTIC_PUSH
  HPP_FCL_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS
  /// @brief Default constructor for GJK algorithm
  /// By default, we don't want EPA to allocate memory because
  /// certain functions of the `GJKSolver` class have specializations
  /// which don't use EPA (and/or GJK).
  /// So we set the maximum number of faces and vertices of EPA to zero.
  /// Only the functions that need EPA will reset the algorithm and allocate
  /// memory if needed.
  GJKSolver()
      : gjk(GJK_DEFAULT_MAX_ITERATIONS, GJK_DEFAULT_TOLERANCE),
        epa(0, 0, EPA_DEFAULT_MAX_ITERATIONS, EPA_DEFAULT_TOLERANCE) {
    epa_max_face_num = EPA_DEFAULT_MAX_FACES;
    epa_max_vertex_num = EPA_DEFAULT_MAX_VERTICES;

    gjk_initial_guess = GJKInitialGuess::DefaultGuess;
    enable_cached_guess = false;  // TODO: use gjk_initial_guess instead
    cached_guess = Vec3f(1, 0, 0);
    support_func_cached_guess = support_func_guess_t::Zero();
    distance_upper_bound = (std::numeric_limits<FCL_REAL>::max)();

    // Default settings for GJK algorithm
    gjk.gjk_variant = GJKVariant::DefaultGJK;
    gjk.convergence_criterion = GJKConvergenceCriterion::VDB;
    gjk.convergence_criterion_type = GJKConvergenceCriterionType::Relative;
  }

  /// @brief Constructor from a DistanceRequest
  ///
  /// \param[in] request DistanceRequest input
  ///
  /// See the default constructor; by default, we don't want
  /// EPA to allocate memory so we set the maximum number of faces and vertices
  /// of EPA to zero.
  GJKSolver(const DistanceRequest& request)
      : gjk(request.gjk_max_iterations, request.gjk_tolerance),
        epa(0, 0, request.epa_max_iterations, request.epa_tolerance) {
    cached_guess = Vec3f(1, 0, 0);
    support_func_cached_guess = support_func_guess_t::Zero();
    distance_upper_bound = (std::numeric_limits<FCL_REAL>::max)();

    set(request);
  }

  /// @brief setter from a DistanceRequest
  ///
  /// \param[in] request DistanceRequest input
  ///
  void set(const DistanceRequest& request) {
    // ---------------------
    // GJK settings
    gjk_initial_guess = request.gjk_initial_guess;
    enable_cached_guess = request.enable_cached_gjk_guess;
    if (gjk_initial_guess == GJKInitialGuess::CachedGuess ||
        enable_cached_guess) {
      cached_guess = request.cached_gjk_guess;
      support_func_cached_guess = request.cached_support_func_guess;
    }
    gjk.tolerance = request.gjk_tolerance;
    gjk.max_iterations = request.gjk_max_iterations;
    // For distance computation, we don't want GJK to early stop
    gjk.setDistanceEarlyBreak(std::numeric_limits<FCL_REAL>::max());
    gjk.gjk_variant = request.gjk_variant;
    gjk.convergence_criterion = request.gjk_convergence_criterion;
    gjk.convergence_criterion_type = request.gjk_convergence_criterion_type;

    // ---------------------
    // EPA settings
    epa_max_face_num = request.epa_max_face_num;
    epa_max_vertex_num = request.epa_max_vertex_num;
    epa.max_iterations = request.epa_max_iterations;
    epa.tolerance = request.epa_tolerance;

    gjk.status = details::GJK::Status::DidNotRun;
    epa.status = details::EPA::Status::DidNotRun;

#ifndef NDEBUG
    // Only in debug mode, to warn the user
    if (request.gjk_tolerance < GJK_MINIMUM_TOLERANCE) {
      std::cout << "WARNING - GJK: using a tolerance (";
      std::cout << request.gjk_tolerance;
      std::cout << ") which is lower than the recommended lowest tolerance (";
      std::cout << GJK_DEFAULT_TOLERANCE;
      std::cout << "). Selecting this tolerance might trigger assertions.\n";
    }
    if (request.epa_tolerance < EPA_MINIMUM_TOLERANCE) {
      std::cout << "WARNING - EPA: using a tolerance (";
      std::cout << request.epa_tolerance;
      std::cout << ") which is lower than the recommended lowest tolerance (";
      std::cout << EPA_MINIMUM_TOLERANCE;
      std::cout << "). Selecting this tolerance might trigger assertions.\n";
    }
#endif
  }

  /// @brief Constructor from a CollisionRequest
  ///
  /// \param[in] request CollisionRequest input
  ///
  /// See the default constructor; by default, we don't want
  /// EPA to allocate memory so we set the maximum number of faces and vertices
  /// of EPA to zero.
  GJKSolver(const CollisionRequest& request)
      : gjk(request.gjk_max_iterations, request.gjk_tolerance),
        epa(0, 0, request.epa_max_iterations, request.epa_tolerance) {
    cached_guess = Vec3f(1, 0, 0);
    support_func_cached_guess = support_func_guess_t::Zero();
    distance_upper_bound = (std::numeric_limits<FCL_REAL>::max)();

    set(request);
  }

  /// @brief setter from a CollisionRequest
  ///
  /// \param[in] request CollisionRequest input
  ///
  void set(const CollisionRequest& request) {
    // ---------------------
    // GJK settings
    gjk_initial_guess = request.gjk_initial_guess;
    // TODO: use gjk_initial_guess instead
    enable_cached_guess = request.enable_cached_gjk_guess;
    if (gjk_initial_guess == GJKInitialGuess::CachedGuess ||
        enable_cached_guess) {
      cached_guess = request.cached_gjk_guess;
      support_func_cached_guess = request.cached_support_func_guess;
    }
    gjk.tolerance = request.gjk_tolerance;
    gjk.max_iterations = request.gjk_max_iterations;
    // The distance upper bound should be at least greater to the requested
    // security margin. Otherwise, we will likely miss some collisions.
    distance_upper_bound = (std::max)(
        0., (std::max)(request.distance_upper_bound, request.security_margin));
    gjk.setDistanceEarlyBreak(distance_upper_bound);
    gjk.gjk_variant = request.gjk_variant;
    gjk.convergence_criterion = request.gjk_convergence_criterion;
    gjk.convergence_criterion_type = request.gjk_convergence_criterion_type;

    // ---------------------
    // EPA settings
    epa_max_face_num = request.epa_max_face_num;
    epa_max_vertex_num = request.epa_max_vertex_num;
    epa.max_iterations = request.epa_max_iterations;
    epa.tolerance = request.epa_tolerance;

    // ---------------------
    // Reset GJK and EPA status
    gjk.status = details::GJK::Status::DidNotRun;
    epa.status = details::EPA::Status::DidNotRun;

#ifndef NDEBUG
    // Only in debug mode, to warn the user
    if (request.gjk_tolerance < GJK_MINIMUM_TOLERANCE) {
      std::cout << "WARNING - GJK: using a tolerance (";
      std::cout << request.gjk_tolerance;
      std::cout << ") which is lower than the recommended lowest tolerance (";
      std::cout << GJK_DEFAULT_TOLERANCE;
      std::cout << "). Selecting this tolerance might trigger assertions.\n";
    }
    if (request.epa_tolerance < EPA_MINIMUM_TOLERANCE) {
      std::cout << "WARNING - EPA: using a tolerance (";
      std::cout << request.epa_tolerance;
      std::cout << ") which is lower than the recommended lowest tolerance (";
      std::cout << EPA_MINIMUM_TOLERANCE;
      std::cout << "). Selecting this tolerance might trigger assertions.\n";
    }
#endif
  }

  /// @brief Copy constructor
  GJKSolver(const GJKSolver& other) = default;

  HPP_FCL_COMPILER_DIAGNOSTIC_PUSH
  HPP_FCL_COMPILER_DIAGNOSTIC_IGNORED_DEPRECECATED_DECLARATIONS
  bool operator==(const GJKSolver& other) const {
    return epa_max_face_num == other.epa_max_face_num &&
           epa_max_vertex_num == other.epa_max_vertex_num &&
           enable_cached_guess ==
               other.enable_cached_guess &&  // TODO: use gjk_initial_guess
                                             // instead
           cached_guess == other.cached_guess &&
           support_func_cached_guess == other.support_func_cached_guess &&
           distance_upper_bound == other.distance_upper_bound &&
           gjk_initial_guess == other.gjk_initial_guess;
  }
  HPP_FCL_COMPILER_DIAGNOSTIC_POP

  bool operator!=(const GJKSolver& other) const { return !(*this == other); }

  /// @brief Whether smart guess can be provided
  /// @Deprecated Use gjk_initial_guess instead
  HPP_FCL_DEPRECATED_MESSAGE("Use gjk_initial_guess instead")
  bool enable_cached_guess;

  /// @brief smart guess
  mutable Vec3f cached_guess;

  /// @brief which warm start to use for GJK
  GJKInitialGuess gjk_initial_guess;

  /// @brief smart guess for the support function
  mutable support_func_guess_t support_func_cached_guess;

  /// @brief Distance above which the GJK solver stops its computations and
  /// processes to an early stopping.
  ///        The two witness points are incorrect, but with the guaranty that
  ///        the two shapes have a distance greather than distance_upper_bound.
  FCL_REAL distance_upper_bound;

  /// @brief maximum number of simplex face used in EPA algorithm
  size_t epa_max_face_num;

  /// @brief maximum number of simplex vertex used in EPA algorithm
  size_t epa_max_vertex_num;

  /// @brief GJK algorithm
  mutable details::GJK gjk;

  /// @brief EPA algorithm
  mutable details::EPA epa;

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

template <>
HPP_FCL_DLLAPI bool GJKSolver::shapeTriangleInteraction(
    const Sphere& s, const Transform3f& tf1, const Vec3f& P1, const Vec3f& P2,
    const Vec3f& P3, const Transform3f& tf2, FCL_REAL& distance, Vec3f& p1,
    Vec3f& p2, Vec3f& normal) const;

template <>
HPP_FCL_DLLAPI bool GJKSolver::shapeTriangleInteraction(
    const Halfspace& s, const Transform3f& tf1, const Vec3f& P1,
    const Vec3f& P2, const Vec3f& P3, const Transform3f& tf2,
    FCL_REAL& distance, Vec3f& p1, Vec3f& p2, Vec3f& normal) const;

template <>
HPP_FCL_DLLAPI bool GJKSolver::shapeTriangleInteraction(
    const Plane& s, const Transform3f& tf1, const Vec3f& P1, const Vec3f& P2,
    const Vec3f& P3, const Transform3f& tf2, FCL_REAL& distance, Vec3f& p1,
    Vec3f& p2, Vec3f& normal) const;

#define SHAPE_INTERSECT_SPECIALIZATION_BASE(S1, S2)                 \
  template <>                                                       \
  HPP_FCL_DLLAPI bool GJKSolver::shapeIntersect<S1, S2>(            \
      const S1& s1, const Transform3f& tf1, const S2& s2,           \
      const Transform3f& tf2, FCL_REAL& distance_lower_bound, bool, \
      Vec3f* contact_points, Vec3f* normal) const

#define SHAPE_INTERSECT_SPECIALIZATION(S1, S2) \
  SHAPE_INTERSECT_SPECIALIZATION_BASE(S1, S2); \
  SHAPE_INTERSECT_SPECIALIZATION_BASE(S2, S1)

SHAPE_INTERSECT_SPECIALIZATION(Sphere, Capsule);
SHAPE_INTERSECT_SPECIALIZATION_BASE(Sphere, Sphere);
SHAPE_INTERSECT_SPECIALIZATION(Sphere, Box);
SHAPE_INTERSECT_SPECIALIZATION(Sphere, Halfspace);
SHAPE_INTERSECT_SPECIALIZATION(Sphere, Plane);

SHAPE_INTERSECT_SPECIALIZATION(Halfspace, Box);
SHAPE_INTERSECT_SPECIALIZATION(Halfspace, Capsule);
SHAPE_INTERSECT_SPECIALIZATION(Halfspace, Cylinder);
SHAPE_INTERSECT_SPECIALIZATION(Halfspace, Cone);
SHAPE_INTERSECT_SPECIALIZATION(Halfspace, Plane);

SHAPE_INTERSECT_SPECIALIZATION(Plane, Box);
SHAPE_INTERSECT_SPECIALIZATION(Plane, Capsule);
SHAPE_INTERSECT_SPECIALIZATION(Plane, Cylinder);
SHAPE_INTERSECT_SPECIALIZATION(Plane, Cone);

#undef SHAPE_INTERSECT_SPECIALIZATION
#undef SHAPE_INTERSECT_SPECIALIZATION_BASE

#define SHAPE_DISTANCE_SPECIALIZATION_BASE(S1, S2)                  \
  template <>                                                       \
  HPP_FCL_DLLAPI bool GJKSolver::shapeDistance<S1, S2>(             \
      const S1& s1, const Transform3f& tf1, const S2& s2,           \
      const Transform3f& tf2, FCL_REAL& dist, Vec3f& p1, Vec3f& p2, \
      Vec3f& normal) const

#define SHAPE_DISTANCE_SPECIALIZATION(S1, S2) \
  SHAPE_DISTANCE_SPECIALIZATION_BASE(S1, S2); \
  SHAPE_DISTANCE_SPECIALIZATION_BASE(S2, S1)

SHAPE_DISTANCE_SPECIALIZATION(Sphere, Capsule);
SHAPE_DISTANCE_SPECIALIZATION(Sphere, Box);
SHAPE_DISTANCE_SPECIALIZATION(Sphere, Cylinder);
SHAPE_DISTANCE_SPECIALIZATION_BASE(Sphere, Sphere);
SHAPE_DISTANCE_SPECIALIZATION_BASE(Capsule, Capsule);
SHAPE_DISTANCE_SPECIALIZATION_BASE(TriangleP, TriangleP);

#undef SHAPE_DISTANCE_SPECIALIZATION
#undef SHAPE_DISTANCE_SPECIALIZATION_BASE

#if !(__cplusplus >= 201103L || (defined(_MSC_VER) && _MSC_VER >= 1600))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc99-extensions"
#endif
/// \name Shape intersection specializations
/// \{

// param doc is the doxygen detailled description (should be enclosed in /** */
// and contain no dot for some obscure reasons).
#define HPP_FCL_DECLARE_SHAPE_INTERSECT(Shape1, Shape2, doc)      \
  /** @brief Fast implementation for Shape1-Shape2 collision. */  \
  doc template <>                                                 \
  HPP_FCL_DLLAPI bool GJKSolver::shapeIntersect<Shape1, Shape2>(  \
      const Shape1& s1, const Transform3f& tf1, const Shape2& s2, \
      const Transform3f& tf2, FCL_REAL& distance_lower_bound,     \
      bool enable_penetration, Vec3f* contact_points, Vec3f* normal) const
#define HPP_FCL_DECLARE_SHAPE_INTERSECT_SELF(Shape, doc) \
  HPP_FCL_DECLARE_SHAPE_INTERSECT(Shape, Shape, doc)
#define HPP_FCL_DECLARE_SHAPE_INTERSECT_PAIR(Shape1, Shape2, doc) \
  HPP_FCL_DECLARE_SHAPE_INTERSECT(Shape1, Shape2, doc);           \
  HPP_FCL_DECLARE_SHAPE_INTERSECT(Shape2, Shape1, doc)

HPP_FCL_DECLARE_SHAPE_INTERSECT_SELF(Sphere, );
HPP_FCL_DECLARE_SHAPE_INTERSECT_PAIR(Sphere, Capsule, );
HPP_FCL_DECLARE_SHAPE_INTERSECT_PAIR(Sphere, Halfspace, );
HPP_FCL_DECLARE_SHAPE_INTERSECT_PAIR(Sphere, Plane, );

template <>
HPP_FCL_DLLAPI bool GJKSolver::shapeIntersect<Box, Sphere>(
    const Box& s1, const Transform3f& tf1, const Sphere& s2,
    const Transform3f& tf2, FCL_REAL& distance_lower_bound,
    bool enable_penetration, Vec3f* contact_points, Vec3f* normal) const;

#ifdef IS_DOXYGEN  // for doxygen only
/** \todo currently disabled and to re-enable it, API of function
 *  \ref obbDisjointAndLowerBoundDistance should be modified.
 *  */
template <>
HPP_FCL_DLLAPI bool GJKSolver::shapeIntersect<Box, Box>(
    const Box& s1, const Transform3f& tf1, const Box& s2,
    const Transform3f& tf2, FCL_REAL& distance_lower_bound,
    bool enable_penetration, Vec3f* contact_points, Vec3f* normal) const;
#endif
// HPP_FCL_DECLARE_SHAPE_INTERSECT_SELF(Box,);
HPP_FCL_DECLARE_SHAPE_INTERSECT_PAIR(Box, Halfspace, );
HPP_FCL_DECLARE_SHAPE_INTERSECT_PAIR(Box, Plane, );

HPP_FCL_DECLARE_SHAPE_INTERSECT_PAIR(Capsule, Halfspace, );
HPP_FCL_DECLARE_SHAPE_INTERSECT_PAIR(Capsule, Plane, );

HPP_FCL_DECLARE_SHAPE_INTERSECT_PAIR(Cylinder, Halfspace, );
HPP_FCL_DECLARE_SHAPE_INTERSECT_PAIR(Cylinder, Plane, );

HPP_FCL_DECLARE_SHAPE_INTERSECT_PAIR(Cone, Halfspace, );
HPP_FCL_DECLARE_SHAPE_INTERSECT_PAIR(Cone, Plane, );

HPP_FCL_DECLARE_SHAPE_INTERSECT_SELF(Halfspace, );

HPP_FCL_DECLARE_SHAPE_INTERSECT_SELF(Plane, );
HPP_FCL_DECLARE_SHAPE_INTERSECT_PAIR(Plane, Halfspace, );

#undef HPP_FCL_DECLARE_SHAPE_INTERSECT
#undef HPP_FCL_DECLARE_SHAPE_INTERSECT_SELF
#undef HPP_FCL_DECLARE_SHAPE_INTERSECT_PAIR

/// \}

/// \name Shape triangle interaction specializations
/// \{

// param doc is the doxygen detailled description (should be enclosed in /** */
// and contain no dot for some obscure reasons).
#define HPP_FCL_DECLARE_SHAPE_TRIANGLE(Shape, doc)                  \
  /** @brief Fast implementation for Shape-Triangle interaction. */ \
  doc template <>                                                   \
  HPP_FCL_DLLAPI bool GJKSolver::shapeTriangleInteraction<Shape>(   \
      const Shape& s, const Transform3f& tf1, const Vec3f& P1,      \
      const Vec3f& P2, const Vec3f& P3, const Transform3f& tf2,     \
      FCL_REAL& distance, Vec3f& p1, Vec3f& p2, Vec3f& normal) const

HPP_FCL_DECLARE_SHAPE_TRIANGLE(Sphere, );
HPP_FCL_DECLARE_SHAPE_TRIANGLE(Halfspace, );
HPP_FCL_DECLARE_SHAPE_TRIANGLE(Plane, );

#undef HPP_FCL_DECLARE_SHAPE_TRIANGLE

/// \}

/// \name Shape distance specializations
/// \{

// param doc is the doxygen detailled description (should be enclosed in /** */
// and contain no dot for some obscure reasons).
#define HPP_FCL_DECLARE_SHAPE_DISTANCE(Shape1, Shape2, doc)         \
  /** @brief Fast implementation for Shape1-Shape2 distance. */     \
  doc template <>                                                   \
  bool HPP_FCL_DLLAPI GJKSolver::shapeDistance<Shape1, Shape2>(     \
      const Shape1& s1, const Transform3f& tf1, const Shape2& s2,   \
      const Transform3f& tf2, FCL_REAL& dist, Vec3f& p1, Vec3f& p2, \
      Vec3f& normal) const
#define HPP_FCL_DECLARE_SHAPE_DISTANCE_SELF(Shape, doc) \
  HPP_FCL_DECLARE_SHAPE_DISTANCE(Shape, Shape, doc)
#define HPP_FCL_DECLARE_SHAPE_DISTANCE_PAIR(Shape1, Shape2, doc) \
  HPP_FCL_DECLARE_SHAPE_DISTANCE(Shape1, Shape2, doc);           \
  HPP_FCL_DECLARE_SHAPE_DISTANCE(Shape2, Shape1, doc)

HPP_FCL_DECLARE_SHAPE_DISTANCE_PAIR(Sphere, Box, );
HPP_FCL_DECLARE_SHAPE_DISTANCE_PAIR(Sphere, Capsule, );
HPP_FCL_DECLARE_SHAPE_DISTANCE_PAIR(Sphere, Cylinder, );
HPP_FCL_DECLARE_SHAPE_DISTANCE_SELF(Sphere, );

HPP_FCL_DECLARE_SHAPE_DISTANCE_SELF(
    Capsule,
    /** Closest points are based on two line-segments. */
);

HPP_FCL_DECLARE_SHAPE_DISTANCE_SELF(
    TriangleP,
    /** Do not run EPA algorithm to compute penetration depth. Use a dedicated
       method. */
);

#undef HPP_FCL_DECLARE_SHAPE_DISTANCE
#undef HPP_FCL_DECLARE_SHAPE_DISTANCE_SELF
#undef HPP_FCL_DECLARE_SHAPE_DISTANCE_PAIR

/// \}
#if !(__cplusplus >= 201103L || (defined(_MSC_VER) && _MSC_VER >= 1600))
#pragma GCC diagnostic pop
#endif
}  // namespace fcl

}  // namespace hpp

#endif
