/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2021-2024, INRIA.
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

#ifndef HPP_FCL_TRAVERSAL_NODE_HFIELD_SHAPE_H
#define HPP_FCL_TRAVERSAL_NODE_HFIELD_SHAPE_H

/// @cond INTERNAL

#include <hpp/fcl/collision_data.h>
#include <hpp/fcl/shape/geometric_shapes.h>
#include <hpp/fcl/narrowphase/narrowphase.h>
#include <hpp/fcl/shape/geometric_shapes_utility.h>
#include <hpp/fcl/internal/traversal_node_base.h>
#include <hpp/fcl/internal/traversal.h>
#include <hpp/fcl/internal/intersect.h>
#include <hpp/fcl/hfield.h>
#include <hpp/fcl/shape/convex.h>

namespace hpp {
namespace fcl {

/// @addtogroup Traversal_For_Collision
/// @{

namespace details {
template <typename BV>
Convex<Quadrilateral> buildConvexQuadrilateral(const HFNode<BV>& node,
                                               const HeightField<BV>& model) {
  const MatrixXf& heights = model.getHeights();
  const VecXf& x_grid = model.getXGrid();
  const VecXf& y_grid = model.getYGrid();

  const FCL_REAL min_height = model.getMinHeight();

  const FCL_REAL x0 = x_grid[node.x_id], x1 = x_grid[node.x_id + 1],
                 y0 = y_grid[node.y_id], y1 = y_grid[node.y_id + 1];
  const Eigen::Block<const MatrixXf, 2, 2> cell =
      heights.block<2, 2>(node.y_id, node.x_id);

  assert(cell.maxCoeff() > min_height &&
         "max_height is lower than min_height");  // Check whether the geometry
                                                  // is degenerated

  std::shared_ptr<std::vector<Vec3f>> pts(new std::vector<Vec3f>({
      Vec3f(x0, y0, min_height),
      Vec3f(x0, y1, min_height),
      Vec3f(x1, y1, min_height),
      Vec3f(x1, y0, min_height),
      Vec3f(x0, y0, cell(0, 0)),
      Vec3f(x0, y1, cell(1, 0)),
      Vec3f(x1, y1, cell(1, 1)),
      Vec3f(x1, y0, cell(0, 1)),
  }));

  std::shared_ptr<std::vector<Quadrilateral>> polygons(
      new std::vector<Quadrilateral>(6));
  (*polygons)[0].set(0, 3, 2, 1);  // x+ side
  (*polygons)[1].set(0, 1, 5, 4);  // y- side
  (*polygons)[2].set(1, 2, 6, 5);  // x- side
  (*polygons)[3].set(2, 3, 7, 6);  // y+ side
  (*polygons)[4].set(3, 0, 4, 7);  // z- side
  (*polygons)[5].set(4, 5, 6, 7);  // z+ side

  return Convex<Quadrilateral>(pts,  // points
                               8,    // num points
                               polygons,
                               6  // number of polygons
  );
}

template <typename BV>
void buildConvexTriangles(const HFNode<BV>& node, const HeightField<BV>& model,
                          Convex<Triangle>& convex1,
                          Convex<Triangle>& convex2) {
  const MatrixXf& heights = model.getHeights();
  const VecXf& x_grid = model.getXGrid();
  const VecXf& y_grid = model.getYGrid();

  const FCL_REAL min_height = model.getMinHeight();

  const FCL_REAL x0 = x_grid[node.x_id], x1 = x_grid[node.x_id + 1],
                 y0 = y_grid[node.y_id], y1 = y_grid[node.y_id + 1];
  const FCL_REAL max_height = node.max_height;
  const Eigen::Block<const MatrixXf, 2, 2> cell =
      heights.block<2, 2>(node.y_id, node.x_id);

  assert(max_height > min_height &&
         "max_height is lower than min_height");  // Check whether the geometry
                                                  // is degenerated
  HPP_FCL_UNUSED_VARIABLE(max_height);

  {
    std::shared_ptr<std::vector<Vec3f>> pts(new std::vector<Vec3f>({
        Vec3f(x0, y0, min_height),  // A
        Vec3f(x0, y1, min_height),  // B
        Vec3f(x1, y0, min_height),  // C
        Vec3f(x0, y0, cell(0, 0)),  // D
        Vec3f(x0, y1, cell(1, 0)),  // E
        Vec3f(x1, y0, cell(0, 1)),  // F
    }));

    std::shared_ptr<std::vector<Triangle>> triangles(
        new std::vector<Triangle>(8));
    (*triangles)[0].set(0, 2, 1);  // bottom
    (*triangles)[1].set(3, 4, 5);  // top
    (*triangles)[2].set(0, 1, 3);  // West 1
    (*triangles)[3].set(3, 1, 4);  // West 2
    (*triangles)[4].set(1, 2, 5);  // South-East 1
    (*triangles)[5].set(1, 5, 4);  // South-East 1
    (*triangles)[6].set(0, 5, 2);  // North 1
    (*triangles)[7].set(5, 0, 3);  // North 2

    convex1.set(pts,  // points
                6,    // num points
                triangles,
                8  // number of polygons
    );
  }

  {
    std::shared_ptr<std::vector<Vec3f>> pts(new std::vector<Vec3f>({
        Vec3f(x0, y1, min_height),  // A
        Vec3f(x1, y1, min_height),  // B
        Vec3f(x1, y0, min_height),  // C
        Vec3f(x0, y1, cell(1, 0)),  // D
        Vec3f(x1, y1, cell(1, 1)),  // E
        Vec3f(x1, y0, cell(0, 1)),  // F
    }));

    std::shared_ptr<std::vector<Triangle>> triangles(
        new std::vector<Triangle>(8));
    (*triangles)[0].set(2, 1, 0);  // bottom
    (*triangles)[1].set(3, 4, 5);  // top
    (*triangles)[2].set(0, 1, 3);  // South 1
    (*triangles)[3].set(3, 1, 4);  // South 2
    (*triangles)[4].set(0, 5, 2);  // North West 1
    (*triangles)[5].set(0, 3, 5);  // North West 2
    (*triangles)[6].set(1, 2, 5);  // East 1
    (*triangles)[7].set(4, 1, 2);  // East 2

    convex2.set(pts,  // points
                6,    // num points
                triangles,
                8  // number of polygons
    );
  }
}

inline Vec3f projectTriangle(const Vec3f& pointA, const Vec3f& pointB,
                             const Vec3f& pointC, const Vec3f& point) {
  const Project::ProjectResult result =
      Project::projectTriangle(pointA, pointB, pointC, point);
  Vec3f res = result.parameterization[0] * pointA +
              result.parameterization[1] * pointB +
              result.parameterization[2] * pointC;

  return res;
}

inline Vec3f projectTetrahedra(const Vec3f& pointA, const Vec3f& pointB,
                               const Vec3f& pointC, const Vec3f& pointD,
                               const Vec3f& point) {
  const Project::ProjectResult result =
      Project::projectTetrahedra(pointA, pointB, pointC, pointD, point);
  Vec3f res = result.parameterization[0] * pointA +
              result.parameterization[1] * pointB +
              result.parameterization[2] * pointC +
              result.parameterization[3] * pointD;

  return res;
}

enum class FaceOrientationBinPart1 {
  BOTTOM = 0,
  TOP = 1,
  SOUTH = 2,
  NORTH_WEST = 4,
  NORTH = 6,
};

enum class FaceOrientationBinPart2 {
  BOTTOM = 0,
  TOP = 1,
  WEST = 2,
  SOUTH_EAST = 4,
  EAST = 6,
};

inline Vec3f computeTriangleNormal(const Triangle& triangle,
                                   const std::vector<Vec3f>& points) {
  const Vec3f pointA = points[triangle[0]];
  const Vec3f pointB = points[triangle[1]];
  const Vec3f pointC = points[triangle[2]];

  const Vec3f normal = (pointB - pointA)
                           .cross(pointC - pointA)
                           .normalized() assert(!normal.array().isNaN().any() &&
                                                "normal is ill-defined");

  return normal;
}

inline Vec3f projectPointOnTriangle(const Vec3f& contact_point,
                                    const Triangle& triangle,
                                    const std::vector<Vec3f>& points) {
  const Vec3f pointA = points[triangle[0]];
  const Vec3f pointB = points[triangle[1]];
  const Vec3f pointC = points[triangle[2]];

  const Vec3f contact_point_projected =
      projectTriangle(pointA, pointB, pointC, contact_point);

  return contact_point_projected;
}

inline bool isContactPointLyingOnTriangle(
    const Vec3f& contact_point, const Triangle& triangle,
    const std::vector<Vec3f>& points, const FCL_REAL prec = FCL_REAL(1e-12)) {
  const Vec3f contact_point_projected =
      projectPointOnTriangle(contact_point, triangle, points);
  return contact_point_projected.isApprox(contact_point, prec);
}

template <char face_id>
bool isContactPointLyingOnFace(const Vec3f& contact_point,
                               const Convex<Triangle>& convex,
                               const FCL_REAL prec = FCL_REAL(1e-12)) {
  const std::vector<Vec3f>& points = *(convex.points);
  if (face_id <= 1) {
    const Triangle& triangle = (*(convex.polygons))[face_id];
    return isContactPointLyingOnTriangle(contact_point, triangle, points, prec);
  } else {
    const Triangle& triangle1 = (*(convex.polygons))[face_id];
    const bool is_lying_on_triangle1 =
        isContactPointLyingOnTriangle(contact_point, triangle1, points, prec);
    if (is_lying_on_triangle1) return true;
    const Triangle& triangle2 = (*(convex.polygons))[face_id + 1];
    const bool is_lying_on_triangle2 =
        isContactPointLyingOnTriangle(contact_point, triangle2, points, prec);
    return is_lying_on_triangle2;
  }
}

template <typename Polygone, typename Shape>
bool binCorrection(const Convex<Polygone>& convex, const Shape& shape,
                   const Transform3f& shape_pose, FCL_REAL& distance,
                   Vec3f& contact_1, Vec3f& contact_2, Vec3f& normal,
                   Vec3f& normal_top, const bool is_collision) {
  const Polygone& top_triangle = (*(convex.polygons))[1];
  const std::vector<Vec3f>& points = *(convex.points);

  normal_top = computeTriangleNormal(top_triangle, points);

  bool hfield_witness_is_on_bin_side =
      isContactPointLyingOnFace<1>(contact_1, convex);

  // We correct only if there is a collision with the bin
  if (is_collision) {
    int hint = 0;
    const Vec3f _support = getSupport(
        &shape, -shape_pose.rotation().transpose() * normal_top, true, hint);
    const Vec3f support =
        shape_pose.rotation() * _support + shape_pose.translation();

    // Project support into the inclined bin having triangle
    const Vec3f pointA = points[top_triangle[0]];
    const FCL_REAL offset_plane = normal_top.dot(pointA);
    const Plane projection_plane(normal_top, offset_plane);
    const FCL_REAL distance_support_projection_plane =
        projection_plane.signedDistance(support);

    const Vec3f projected_support =
        support - distance_support_projection_plane * normal_top;

    // We need now to project the projected in the triangle shape
    contact_1 = projectPointOnTriangle(projected_support, top_triangle, points);
    contact_2 = contact_1 + distance_support_projection_plane * normal_top;
    normal = normal_top;
    distance = -distance_support_projection_plane;
  }

  return hfield_witness_is_on_bin_side;
}

template <typename Polygone, typename Shape, int Options>
bool shapeDistance(const GJKSolver* nsolver, const Convex<Polygone>& convex1,
                   const Convex<Polygone>& convex2, const Transform3f& tf1,
                   const Shape& shape, const Transform3f& tf2,
                   FCL_REAL& distance, Vec3f& c1, Vec3f& c2, Vec3f& normal,
                   Vec3f& normal_top, bool& hfield_witness_is_on_bin_side) {
  enum { RTIsIdentity = Options & RelativeTransformationIsIdentity };

  const Transform3f Id;
  Vec3f contact1_1, contact1_2, normal1, normal1_top;
  Vec3f contact2_1, contact2_2, normal2, normal2_top;
  FCL_REAL distance1, distance2;
  bool collision1, collision2;
  bool hfield_witness_is_on_bin_side1, hfield_witness_is_on_bin_side2;
  // Because we do bin comparison and correction, we always compute penetration
  // information.
  bool compute_penetration = true;
  if (RTIsIdentity)
    nsolver->shapeDistance(convex1, Id, shape, tf2, distance1,
                           compute_penetration, contact1_1, contact1_2,
                           normal1);
  else
    nsolver->shapeDistance(convex1, tf1, shape, tf2, distance1,
                           compute_penetration, contact1_1, contact1_2,
                           normal1);
  collision1 = (distance1 < 0);

  hfield_witness_is_on_bin_side1 =
      binCorrection(convex1, shape, tf2, distance1, contact1_1, contact1_2,
                    normal1, normal1_top, collision1);

  if (RTIsIdentity)
    nsolver->shapeDistance(convex2, Id, shape, tf2, distance2,
                           compute_penetration, contact2_1, contact2_2, normal);
  else
    nsolver->shapeDistance(convex2, tf1, shape, tf2, distance2,
                           compute_penetration, contact2_1, contact2_2,
                           normal2);
  collision2 = (distance2 < 0);

  hfield_witness_is_on_bin_side2 =
      binCorrection(convex2, shape, tf2, distance2, contact2_1, contact2_2,
                    normal2, normal2_top, collision2);

  if (collision1 && collision2) {
    if (distance1 > distance2)  // switch values
    {
      distance = distance2;
      c1 = contact2_1;
      c2 = contact2_2;
      normal = normal2;
      normal_top = normal2_top;
      hfield_witness_is_on_bin_side = hfield_witness_is_on_bin_side2;
    } else {
      distance = distance1;
      c1 = contact1_1;
      c2 = contact1_2;
      normal = normal1;
      normal_top = normal1_top;
      hfield_witness_is_on_bin_side = hfield_witness_is_on_bin_side1;
    }
    return true;
  } else if (collision1) {
    distance = distance1;
    c1 = contact1_1;
    c2 = contact1_2;
    normal = normal1;
    normal_top = normal1_top;
    hfield_witness_is_on_bin_side = hfield_witness_is_on_bin_side1;
    return true;
  } else if (collision2) {
    distance = distance2;
    c1 = contact2_1;
    c2 = contact2_2;
    normal = normal2;
    normal_top = normal2_top;
    hfield_witness_is_on_bin_side = hfield_witness_is_on_bin_side2;
    return true;
  }

  if (distance1 > distance2)  // switch values
  {
    distance = distance2;
    c1 = contact2_1;
    c2 = contact2_2;
    normal = normal2;
    normal_top = normal2_top;
    hfield_witness_is_on_bin_side = hfield_witness_is_on_bin_side2;
  } else {
    distance = distance1;
    c1 = contact1_1;
    c2 = contact1_2;
    normal = normal1;
    normal_top = normal1_top;
    hfield_witness_is_on_bin_side = hfield_witness_is_on_bin_side1;
  }
  return false;
}

template <typename Polygone, typename Shape, int Options>
bool shapeCollision(const GJKSolver* nsolver, const Convex<Polygone>& convex1,
                    const Convex<Polygone>& convex2, const Transform3f& tf1,
                    const Shape& shape, const Transform3f& tf2,
                    FCL_REAL& distance_lower_bound, Vec3f& contact_point,
                    Vec3f& normal) {
  enum { RTIsIdentity = Options & RelativeTransformationIsIdentity };

  const Transform3f Id;
  Vec3f contact_point2, normal2;
  FCL_REAL distance_lower_bound2;
  bool collision1, collision2;
  if (RTIsIdentity)
    collision1 =
        nsolver->shapeIntersect(convex1, Id, shape, tf2, distance_lower_bound,
                                true, &contact_point, &normal);
  else
    collision1 =
        nsolver->shapeIntersect(convex1, tf1, shape, tf2, distance_lower_bound,
                                true, &contact_point, &normal);

  if (RTIsIdentity)
    collision2 =
        nsolver->shapeIntersect(convex2, Id, shape, tf2, distance_lower_bound2,
                                true, &contact_point2, &normal2);
  else
    collision2 =
        nsolver->shapeIntersect(convex2, tf1, shape, tf2, distance_lower_bound2,
                                true, &contact_point2, &normal2);

  if (collision1 && collision2) {
    // In some case, EPA might returns something like
    // -(std::numeric_limits<FCL_REAL>::max)().
    if (distance_lower_bound != -(std::numeric_limits<FCL_REAL>::max)() &&
        distance_lower_bound2 != -(std::numeric_limits<FCL_REAL>::max)()) {
      if (distance_lower_bound > distance_lower_bound2)  // switch values
      {
        distance_lower_bound = distance_lower_bound2;
        contact_point = contact_point2;
        normal = normal2;
      }
    } else if (distance_lower_bound2 !=
               -(std::numeric_limits<FCL_REAL>::max)()) {
      distance_lower_bound = distance_lower_bound2;
      contact_point = contact_point2;
      normal = normal2;
    }
    return true;
  } else if (collision1) {
    return true;
  } else if (collision2) {
    distance_lower_bound = distance_lower_bound2;
    contact_point = contact_point2;
    normal = normal2;
    return true;
  }

  return false;
}
}  // namespace details

/// @brief Traversal node for collision between height field and shape
template <typename BV, typename S,
          int _Options = RelativeTransformationIsIdentity>
class HeightFieldShapeCollisionTraversalNode
    : public CollisionTraversalNodeBase {
 public:
  typedef CollisionTraversalNodeBase Base;
  typedef Eigen::Array<FCL_REAL, 1, 2> Array2d;

  enum {
    Options = _Options,
    RTIsIdentity = _Options & RelativeTransformationIsIdentity
  };

  HeightFieldShapeCollisionTraversalNode(const CollisionRequest& request)
      : CollisionTraversalNodeBase(request) {
    model1 = NULL;
    model2 = NULL;

    num_bv_tests = 0;
    num_leaf_tests = 0;
    query_time_seconds = 0.0;

    nsolver = NULL;
    shape_inflation.setZero();
    count = 0;
  }

  /// @brief Whether the BV node in the first BVH tree is leaf
  bool isFirstNodeLeaf(unsigned int b) const {
    return model1->getBV(b).isLeaf();
  }

  /// @brief Obtain the left child of BV node in the first BVH
  int getFirstLeftChild(unsigned int b) const {
    return static_cast<int>(model1->getBV(b).leftChild());
  }

  /// @brief Obtain the right child of BV node in the first BVH
  int getFirstRightChild(unsigned int b) const {
    return static_cast<int>(model1->getBV(b).rightChild());
  }

  /// test between BV b1 and shape
  /// @param b1 BV to test,
  /// @retval sqrDistLowerBound square of a lower bound of the minimal
  ///         distance between bounding volumes.
  /// @brief BV culling test in one BVTT node
  bool BVDisjoints(unsigned int b1, unsigned int /*b2*/,
                   FCL_REAL& sqrDistLowerBound) const {
    if (this->enable_statistics) this->num_bv_tests++;

    bool disjoint;
    if (RTIsIdentity) {
      assert(false && "must never happened");
      disjoint = !this->model1->getBV(b1).bv.overlap(
          this->model2_bv, this->request, sqrDistLowerBound);
    } else {
      disjoint = !overlap(this->tf1.getRotation(), this->tf1.getTranslation(),
                          this->model1->getBV(b1).bv, this->model2_bv,
                          this->request, sqrDistLowerBound);
    }

    if (disjoint)
      internal::updateDistanceLowerBoundFromBV(this->request, *this->result,
                                               sqrDistLowerBound);

    assert(!disjoint || sqrDistLowerBound > 0);
    return disjoint;
  }

  /// @brief Intersection testing between leaves (one Convex and one shape)
  void leafCollides(unsigned int b1, unsigned int /*b2*/,
                    FCL_REAL& sqrDistLowerBound) const {
    count++;
    if (this->enable_statistics) this->num_leaf_tests++;
    const HFNode<BV>& node = this->model1->getBV(b1);

    // Split quadrilateral primitives into two convex shapes corresponding to
    // polyhedron with triangular bases. This is essential to keep the convexity

    //    typedef Convex<Quadrilateral> ConvexQuadrilateral;
    //    const ConvexQuadrilateral convex =
    //    details::buildConvexQuadrilateral(node,*this->model1);

    typedef Convex<Triangle> ConvexTriangle;
    ConvexTriangle convex1, convex2;
    details::buildConvexTriangles(node, *this->model1, convex1, convex2);

    // Compute aabb_local for BoundingVolumeGuess case in the GJK solver
    if (nsolver->gjk_initial_guess == GJKInitialGuess::BoundingVolumeGuess) {
      convex1.computeLocalAABB();
      convex2.computeLocalAABB();
    }

    FCL_REAL distance;
    //    Vec3f contact_point, normal;
    Vec3f c1, c2, normal, normal_top;
    bool hfield_witness_is_on_bin_side;

    bool collision = details::shapeDistance<Triangle, S, Options>(
        nsolver, convex1, convex2, this->tf1, *(this->model2), this->tf2,
        distance, c1, c2, normal, normal_top, hfield_witness_is_on_bin_side);

    FCL_REAL distToCollision =
        distance - this->request.security_margin * (normal_top.dot(normal));
    if (distToCollision <= this->request.collision_distance_threshold) {
      sqrDistLowerBound = 0;
      if (this->result->numContacts() < this->request.num_max_contacts) {
        if (normal_top.isApprox(normal) &&
            (collision || !hfield_witness_is_on_bin_side)) {
          this->result->addContact(Contact(this->model1, this->model2, (int)b1,
                                           (int)Contact::NONE, c1, c2, normal,
                                           distance));
          assert(this->result->isCollision());
        }
      }
    } else
      sqrDistLowerBound = distToCollision * distToCollision;

    //    const Vec3f c1 = contact_point - distance * 0.5 * normal;
    //    const Vec3f c2 = contact_point + distance * 0.5 * normal;
    internal::updateDistanceLowerBoundFromLeaf(this->request, *this->result,
                                               distToCollision, c1, c2, normal);

    assert(this->result->isCollision() || sqrDistLowerBound > 0);
  }

  const GJKSolver* nsolver;

  const HeightField<BV>* model1;
  const S* model2;
  BV model2_bv;

  Array2d shape_inflation;

  mutable int num_bv_tests;
  mutable int num_leaf_tests;
  mutable FCL_REAL query_time_seconds;
  mutable int count;
};

/// @}

/// @addtogroup Traversal_For_Distance
/// @{

/// @brief Traversal node for distance between height field and shape
template <typename BV, typename S,
          int _Options = RelativeTransformationIsIdentity>
class HeightFieldShapeDistanceTraversalNode : public DistanceTraversalNodeBase {
 public:
  typedef DistanceTraversalNodeBase Base;

  enum {
    Options = _Options,
    RTIsIdentity = _Options & RelativeTransformationIsIdentity
  };

  HeightFieldShapeDistanceTraversalNode() : DistanceTraversalNodeBase() {
    model1 = NULL;
    model2 = NULL;

    num_leaf_tests = 0;
    query_time_seconds = 0.0;

    rel_err = 0;
    abs_err = 0;
    nsolver = NULL;
  }

  /// @brief Whether the BV node in the first BVH tree is leaf
  bool isFirstNodeLeaf(unsigned int b) const {
    return model1->getBV(b).isLeaf();
  }

  /// @brief Obtain the left child of BV node in the first BVH
  int getFirstLeftChild(unsigned int b) const {
    return model1->getBV(b).leftChild();
  }

  /// @brief Obtain the right child of BV node in the first BVH
  int getFirstRightChild(unsigned int b) const {
    return model1->getBV(b).rightChild();
  }

  /// @brief BV culling test in one BVTT node
  FCL_REAL BVDistanceLowerBound(unsigned int b1, unsigned int /*b2*/) const {
    return model1->getBV(b1).bv.distance(
        model2_bv);  // TODO(jcarpent): tf1 is not taken into account here.
  }

  /// @brief Distance testing between leaves (one triangle and one shape)
  void leafComputeDistance(unsigned int b1, unsigned int /*b2*/) const {
    if (this->enable_statistics) this->num_leaf_tests++;

    const BVNode<BV>& node = this->model1->getBV(b1);

    typedef Convex<Quadrilateral> ConvexQuadrilateral;
    const ConvexQuadrilateral convex =
        details::buildConvexQuadrilateral(node, *this->model1);

    FCL_REAL d;
    Vec3f closest_p1, closest_p2, normal;

    nsolver->shapeDistance(convex, this->tf1, *(this->model2), this->tf2, d,
                           closest_p1, closest_p2, normal);

    this->result->update(d, this->model1, this->model2, b1,
                         DistanceResult::NONE, closest_p1, closest_p2, normal);
  }

  /// @brief Whether the traversal process can stop early
  bool canStop(FCL_REAL c) const {
    if ((c >= this->result->min_distance - abs_err) &&
        (c * (1 + rel_err) >= this->result->min_distance))
      return true;
    return false;
  }

  FCL_REAL rel_err;
  FCL_REAL abs_err;

  const GJKSolver* nsolver;

  const HeightField<BV>* model1;
  const S* model2;
  BV model2_bv;

  mutable int num_bv_tests;
  mutable int num_leaf_tests;
  mutable FCL_REAL query_time_seconds;
};

/// @}

}  // namespace fcl
}  // namespace hpp

/// @endcond

#endif
