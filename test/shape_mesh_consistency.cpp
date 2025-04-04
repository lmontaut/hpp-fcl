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

#define BOOST_TEST_MODULE COAL_SHAPE_MESH_CONSISTENCY
#include <boost/test/included/unit_test.hpp>

#include "coal/narrowphase/narrowphase.h"
#include "coal/shape/geometric_shape_to_BVH_model.h"
#include "coal/distance.h"
#include "coal/collision.h"
#include "utility.h"

using namespace coal;

#define BOOST_CHECK_FALSE(p) BOOST_CHECK(!(p))

Scalar extents[6] = {0, 0, 0, 10, 10, 10};

BOOST_AUTO_TEST_CASE(consistency_distance_spheresphere) {
  Sphere s1(20);
  Sphere s2(20);
  BVHModel<RSS> s1_rss;
  BVHModel<RSS> s2_rss;

  generateBVHModel(s1_rss, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_rss, s2, Transform3s(), 16, 16);

  DistanceRequest request;
  DistanceResult res, res1;

  Transform3s pose;

  pose.setTranslation(Vec3s(50, 0, 0));

  res.clear();
  res1.clear();
  distance(&s1, Transform3s(), &s2, pose, request, res);
  distance(&s1_rss, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.05);

  res1.clear();
  distance(&s1, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.05);

  res1.clear();
  distance(&s1_rss, Transform3s(), &s2, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.05);

  for (std::size_t i = 0; i < 10; ++i) {
    Transform3s t;
    generateRandomTransform(extents, t);

    Transform3s pose1(t);
    Transform3s pose2 = t * pose;

    res.clear();
    res1.clear();
    distance(&s1, pose1, &s2, pose2, request, res);
    distance(&s1_rss, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                0.05);

    res1.clear();
    distance(&s1, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                0.05);

    res1.clear();
    distance(&s1_rss, pose1, &s2, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                0.05);
  }

  pose.setTranslation(Vec3s(40.1, 0, 0));

  res.clear(), res1.clear();
  distance(&s1, Transform3s(), &s2, pose, request, res);
  distance(&s1_rss, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  res1.clear();
  distance(&s1, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  res1.clear();
  distance(&s1_rss, Transform3s(), &s2, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  for (std::size_t i = 0; i < 10; ++i) {
    Transform3s t;
    generateRandomTransform(extents, t);

    Transform3s pose1(t);
    Transform3s pose2 = t * pose;

    res.clear();
    res1.clear();
    distance(&s1, pose1, &s2, pose2, request, res);
    distance(&s1_rss, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);

    res1.clear();
    distance(&s1, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);

    res1.clear();
    distance(&s1_rss, pose1, &s2, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);
  }
}

BOOST_AUTO_TEST_CASE(consistency_distance_boxbox) {
  Box s1(20, 40, 50);
  Box s2(10, 10, 10);

  BVHModel<RSS> s1_rss;
  BVHModel<RSS> s2_rss;

  generateBVHModel(s1_rss, s1, Transform3s());
  generateBVHModel(s2_rss, s2, Transform3s());

  DistanceRequest request;
  DistanceResult res, res1;

  Transform3s pose;

  pose.setTranslation(Vec3s(50, 0, 0));

  res.clear();
  res1.clear();
  distance(&s1, Transform3s(), &s2, pose, request, res);
  distance(&s1_rss, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.01);

  res1.clear();
  distance(&s1, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.01);

  res1.clear();
  distance(&s1_rss, Transform3s(), &s2, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.01);

  for (std::size_t i = 0; i < 10; ++i) {
    Transform3s t;
    generateRandomTransform(extents, t);

    Transform3s pose1(t);
    Transform3s pose2 = t * pose;

    res.clear();
    res1.clear();
    distance(&s1, pose1, &s2, pose2, request, res);
    distance(&s1_rss, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                0.01);

    res1.clear();
    distance(&s1, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                0.01);

    res1.clear();
    distance(&s1_rss, pose1, &s2, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                0.01);
  }

  pose.setTranslation(Vec3s(15.1, 0, 0));

  res.clear();
  res1.clear();
  distance(&s1, Transform3s(), &s2, pose, request, res);
  distance(&s1_rss, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  res1.clear();
  distance(&s1, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  res1.clear();
  distance(&s1_rss, Transform3s(), &s2, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  for (std::size_t i = 0; i < 10; ++i) {
    Transform3s t;
    generateRandomTransform(extents, t);

    Transform3s pose1(t);
    Transform3s pose2 = t * pose;

    res.clear();
    res1.clear();
    distance(&s1, pose1, &s2, pose2, request, res);
    distance(&s1_rss, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);

    res1.clear();
    distance(&s1, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);

    res1.clear();
    distance(&s1_rss, pose1, &s2, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);
  }
}

BOOST_AUTO_TEST_CASE(consistency_distance_cylindercylinder) {
  Cylinder s1(5, 10);
  Cylinder s2(5, 10);

  BVHModel<RSS> s1_rss;
  BVHModel<RSS> s2_rss;

  generateBVHModel(s1_rss, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_rss, s2, Transform3s(), 16, 16);

  DistanceRequest request;
  DistanceResult res, res1;

  Transform3s pose;

  pose.setTranslation(Vec3s(20, 0, 0));

  res.clear();
  res1.clear();
  distance(&s1, Transform3s(), &s2, pose, request, res);
  distance(&s1_rss, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.01);

  res1.clear();
  distance(&s1, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.01);

  res1.clear();
  distance(&s1_rss, Transform3s(), &s2, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.01);

  for (std::size_t i = 0; i < 10; ++i) {
    Transform3s t;
    generateRandomTransform(extents, t);

    Transform3s pose1(t);
    Transform3s pose2 = t * pose;

    res.clear();
    res1.clear();
    distance(&s1, pose1, &s2, pose2, request, res);
    distance(&s1_rss, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                0.01);

    res1.clear();
    distance(&s1, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                0.01);

    res1.clear();
    distance(&s1_rss, pose1, &s2, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                0.01);
  }

  pose.setTranslation(
      Vec3s(15, 0, 0));  // libccd cannot use small value here :(

  res.clear();
  res1.clear();
  distance(&s1, Transform3s(), &s2, pose, request, res);
  distance(&s1_rss, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  res1.clear();
  distance(&s1, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  res1.clear();
  distance(&s1_rss, Transform3s(), &s2, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  for (std::size_t i = 0; i < 10; ++i) {
    Transform3s t;
    generateRandomTransform(extents, t);

    Transform3s pose1(t);
    Transform3s pose2 = t * pose;

    res.clear();
    res1.clear();
    distance(&s1, pose1, &s2, pose2, request, res);
    distance(&s1_rss, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);

    res1.clear();
    distance(&s1, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);

    res1.clear();
    distance(&s1_rss, pose1, &s2, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);
  }
}

BOOST_AUTO_TEST_CASE(consistency_distance_conecone) {
  Cone s1(5, 10);
  Cone s2(5, 10);

  BVHModel<RSS> s1_rss;
  BVHModel<RSS> s2_rss;

  generateBVHModel(s1_rss, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_rss, s2, Transform3s(), 16, 16);

  DistanceRequest request;
  DistanceResult res, res1;

  Transform3s pose;

  pose.setTranslation(Vec3s(20, 0, 0));

  res.clear();
  res1.clear();
  distance(&s1, Transform3s(), &s2, pose, request, res);
  distance(&s1_rss, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.05);

  res1.clear();
  distance(&s1, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.05);

  distance(&s1_rss, Transform3s(), &s2, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.05);

  for (std::size_t i = 0; i < 10; ++i) {
    Transform3s t;
    generateRandomTransform(extents, t);

    Transform3s pose1(t);
    Transform3s pose2 = t * pose;

    res.clear();
    res1.clear();
    distance(&s1, pose1, &s2, pose2, request, res);
    distance(&s1_rss, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                0.05);

    res1.clear();
    distance(&s1, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                0.05);

    res1.clear();
    distance(&s1_rss, pose1, &s2, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                0.05);
  }

  pose.setTranslation(
      Vec3s(15, 0, 0));  // libccd cannot use small value here :(

  res.clear();
  res1.clear();
  distance(&s1, Transform3s(), &s2, pose, request, res);
  distance(&s1_rss, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  res1.clear();
  distance(&s1, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  res1.clear();
  distance(&s1_rss, Transform3s(), &s2, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  for (std::size_t i = 0; i < 10; ++i) {
    Transform3s t;
    generateRandomTransform(extents, t);

    Transform3s pose1(t);
    Transform3s pose2 = t * pose;

    res.clear();
    res1.clear();
    distance(&s1, pose1, &s2, pose2, request, res);
    distance(&s1_rss, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);

    res1.clear();
    distance(&s1, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);

    res1.clear();
    distance(&s1_rss, pose1, &s2, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);
  }
}

BOOST_AUTO_TEST_CASE(consistency_distance_spheresphere_GJK) {
  Sphere s1(20);
  Sphere s2(20);
  BVHModel<RSS> s1_rss;
  BVHModel<RSS> s2_rss;

  generateBVHModel(s1_rss, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_rss, s2, Transform3s(), 16, 16);

  DistanceRequest request;
  DistanceResult res, res1;

  Transform3s pose;

  pose.setTranslation(Vec3s(50, 0, 0));

  res.clear();
  res1.clear();
  distance(&s1, Transform3s(), &s2, pose, request, res);
  distance(&s1_rss, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.05);

  res1.clear();
  distance(&s1, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.05);

  res1.clear();
  distance(&s1_rss, Transform3s(), &s2, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.05);

  for (std::size_t i = 0; i < 10; ++i) {
    Transform3s t;
    generateRandomTransform(extents, t);

    Transform3s pose1(t);
    Transform3s pose2 = t * pose;

    res.clear();
    res1.clear();
    distance(&s1, pose1, &s2, pose2, request, res);
    distance(&s1_rss, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                0.05);

    res1.clear();
    distance(&s1, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                0.05);

    res1.clear();
    distance(&s1_rss, pose1, &s2, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                0.05);
  }

  pose.setTranslation(Vec3s(40.1, 0, 0));

  res.clear();
  res1.clear();
  distance(&s1, Transform3s(), &s2, pose, request, res);
  distance(&s1_rss, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  res1.clear();
  distance(&s1, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  res1.clear();
  distance(&s1_rss, Transform3s(), &s2, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              4);

  for (std::size_t i = 0; i < 10; ++i) {
    Transform3s t;
    generateRandomTransform(extents, t);

    Transform3s pose1(t);
    Transform3s pose2 = t * pose;

    res.clear();
    res1.clear();
    distance(&s1, pose1, &s2, pose2, request, res);
    distance(&s1_rss, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);

    res1.clear();
    distance(&s1, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);

    res1.clear();
    distance(&s1_rss, pose1, &s2, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                4);
  }
}

BOOST_AUTO_TEST_CASE(consistency_distance_boxbox_GJK) {
  Box s1(20, 40, 50);
  Box s2(10, 10, 10);

  BVHModel<RSS> s1_rss;
  BVHModel<RSS> s2_rss;

  generateBVHModel(s1_rss, s1, Transform3s());
  generateBVHModel(s2_rss, s2, Transform3s());

  DistanceRequest request;
  DistanceResult res, res1;

  Transform3s pose;

  pose.setTranslation(Vec3s(50, 0, 0));

  res.clear();
  res1.clear();
  distance(&s1, Transform3s(), &s2, pose, request, res);
  distance(&s1_rss, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.01);

  res1.clear();
  distance(&s1, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.01);

  res1.clear();
  distance(&s1_rss, Transform3s(), &s2, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.01);

  for (std::size_t i = 0; i < 10; ++i) {
    Transform3s t;
    generateRandomTransform(extents, t);

    Transform3s pose1(t);
    Transform3s pose2 = t * pose;

    res.clear();
    res1.clear();
    distance(&s1, pose1, &s2, pose2, request, res);
    distance(&s1_rss, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                0.01);

    res1.clear();
    distance(&s1, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                0.01);

    res1.clear();
    distance(&s1_rss, pose1, &s2, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                0.01);
  }

  pose.setTranslation(Vec3s(15.1, 0, 0));

  res.clear();
  res1.clear();
  distance(&s1, Transform3s(), &s2, pose, request, res);
  distance(&s1_rss, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  res1.clear();
  distance(&s1, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  res1.clear();
  distance(&s1_rss, Transform3s(), &s2, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  for (std::size_t i = 0; i < 10; ++i) {
    Transform3s t;
    generateRandomTransform(extents, t);

    Transform3s pose1(t);
    Transform3s pose2 = t * pose;

    res.clear();
    res1.clear();
    distance(&s1, pose1, &s2, pose2, request, res);
    distance(&s1_rss, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);

    res1.clear();
    distance(&s1, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);

    res1.clear();
    distance(&s1_rss, pose1, &s2, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);
  }
}

BOOST_AUTO_TEST_CASE(consistency_distance_cylindercylinder_GJK) {
  Cylinder s1(5, 10);
  Cylinder s2(5, 10);

  BVHModel<RSS> s1_rss;
  BVHModel<RSS> s2_rss;

  generateBVHModel(s1_rss, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_rss, s2, Transform3s(), 16, 16);

  DistanceRequest request;
  DistanceResult res, res1;

  Transform3s pose;

  pose.setTranslation(Vec3s(20, 0, 0));

  res.clear();
  res1.clear();
  distance(&s1, Transform3s(), &s2, pose, request, res);
  distance(&s1_rss, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.05);

  res1.clear();
  distance(&s1, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.05);

  res1.clear();
  distance(&s1_rss, Transform3s(), &s2, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.05);

  for (std::size_t i = 0; i < 10; ++i) {
    Transform3s t;
    generateRandomTransform(extents, t);

    Transform3s pose1(t);
    Transform3s pose2 = t * pose;

    res.clear();
    res1.clear();
    distance(&s1, pose1, &s2, pose2, request, res);
    distance(&s1_rss, pose1, &s2_rss, pose2, request, res1);
    if (fabs(res1.min_distance - res.min_distance) / res.min_distance > 0.05)
      std::cout << "low resolution: " << res1.min_distance << " "
                << res.min_distance << std::endl;
    else
      BOOST_CHECK(
          fabs(res1.min_distance - res.min_distance) / res.min_distance < 0.05);

    res1.clear();
    distance(&s1, pose1, &s2_rss, pose2, request, res1);
    if (fabs(res1.min_distance - res.min_distance) / res.min_distance > 0.05)
      std::cout << "low resolution: " << res1.min_distance << " "
                << res.min_distance << std::endl;
    else
      BOOST_CHECK(
          fabs(res1.min_distance - res.min_distance) / res.min_distance < 0.05);

    res1.clear();
    distance(&s1_rss, pose1, &s2, pose2, request, res1);
    if (fabs(res1.min_distance - res.min_distance) / res.min_distance > 0.05)
      std::cout << "low resolution: " << res1.min_distance << " "
                << res.min_distance << std::endl;
    else
      BOOST_CHECK(
          fabs(res1.min_distance - res.min_distance) / res.min_distance < 0.05);
  }

  pose.setTranslation(Vec3s(10.1, 0, 0));

  res.clear();
  res1.clear();
  distance(&s1, Transform3s(), &s2, pose, request, res);
  distance(&s1_rss, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  res1.clear();
  distance(&s1, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  res1.clear();
  distance(&s1_rss, Transform3s(), &s2, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  for (std::size_t i = 0; i < 10; ++i) {
    Transform3s t;
    generateRandomTransform(extents, t);

    Transform3s pose1(t);
    Transform3s pose2 = t * pose;

    res.clear();
    res1.clear();
    distance(&s1, pose1, &s2, pose2, request, res);
    distance(&s1_rss, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);

    res1.clear();
    distance(&s1, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);

    res1.clear();
    distance(&s1_rss, pose1, &s2, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);
  }
}

BOOST_AUTO_TEST_CASE(consistency_distance_conecone_GJK) {
  Cone s1(5, 10);
  Cone s2(5, 10);

  BVHModel<RSS> s1_rss;
  BVHModel<RSS> s2_rss;

  generateBVHModel(s1_rss, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_rss, s2, Transform3s(), 16, 16);

  DistanceRequest request;
  DistanceResult res, res1;

  Transform3s pose;

  pose.setTranslation(Vec3s(20, 0, 0));

  res.clear();
  res1.clear();
  distance(&s1, Transform3s(), &s2, pose, request, res);
  distance(&s1_rss, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.05);

  res1.clear();
  distance(&s1, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.05);

  res1.clear();
  distance(&s1_rss, Transform3s(), &s2, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              0.05);

  for (std::size_t i = 0; i < 10; ++i) {
    Transform3s t;
    generateRandomTransform(extents, t);

    Transform3s pose1(t);
    Transform3s pose2 = t * pose;

    res.clear();
    res1.clear();
    distance(&s1, pose1, &s2, pose2, request, res);
    distance(&s1_rss, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                0.05);

    res1.clear();
    distance(&s1, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                0.05);

    res1.clear();
    distance(&s1_rss, pose1, &s2, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                0.05);
  }

  pose.setTranslation(Vec3s(10.1, 0, 0));

  res.clear();
  res1.clear();
  distance(&s1, Transform3s(), &s2, pose, request, res);
  distance(&s1_rss, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  res1.clear();
  distance(&s1, Transform3s(), &s2_rss, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  res1.clear();
  distance(&s1_rss, Transform3s(), &s2, pose, request, res1);
  BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
              2);

  for (std::size_t i = 0; i < 10; ++i) {
    Transform3s t;
    generateRandomTransform(extents, t);

    Transform3s pose1(t);
    Transform3s pose2 = t * pose;

    res.clear();
    res1.clear();
    distance(&s1, pose1, &s2, pose2, request, res);
    distance(&s1_rss, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);

    res1.clear();
    distance(&s1, pose1, &s2_rss, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);

    res1.clear();
    distance(&s1_rss, pose1, &s2, pose2, request, res1);
    BOOST_CHECK(fabs(res1.min_distance - res.min_distance) / res.min_distance <
                2);
  }
}

BOOST_AUTO_TEST_CASE(consistency_collision_spheresphere) {
  Sphere s1(20);
  Sphere s2(10);
  BVHModel<AABB> s1_aabb;
  BVHModel<AABB> s2_aabb;
  BVHModel<OBB> s1_obb;
  BVHModel<OBB> s2_obb;

  generateBVHModel(s1_aabb, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_aabb, s2, Transform3s(), 16, 16);
  generateBVHModel(s1_obb, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_obb, s2, Transform3s(), 16, 16);

  CollisionRequest request(false, 1, false);
  CollisionResult result;

  bool res;

  Transform3s pose, pose_aabb, pose_obb;

  // s2 is within s1
  // both are shapes --> collision
  // s1 is shape, s2 is mesh --> in collision
  // s1 is mesh, s2 is shape --> collision free
  // all are reasonable
  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);

  pose.setTranslation(Vec3s(40, 0, 0));
  pose_aabb.setTranslation(Vec3s(40, 0, 0));
  pose_obb.setTranslation(Vec3s(40, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);

  pose.setTranslation(Vec3s(30, 0, 0));
  pose_aabb.setTranslation(Vec3s(30, 0, 0));
  pose_obb.setTranslation(Vec3s(30, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);

  pose.setTranslation(Vec3s(29.9, 0, 0));
  pose_aabb.setTranslation(
      Vec3s(29.8, 0, 0));  // 29.9 fails, result depends on mesh precision
  pose_obb.setTranslation(
      Vec3s(29.8, 0, 0));  // 29.9 fails, result depends on mesh precision

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);

  pose.setTranslation(Vec3s(-29.9, 0, 0));
  pose_aabb.setTranslation(
      Vec3s(-29.8, 0, 0));  // 29.9 fails, result depends on mesh precision
  pose_obb.setTranslation(
      Vec3s(-29.8, 0, 0));  // 29.9 fails, result depends on mesh precision

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);

  pose.setTranslation(Vec3s(-30, 0, 0));
  pose_aabb.setTranslation(Vec3s(-30, 0, 0));
  pose_obb.setTranslation(Vec3s(-30, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
}

BOOST_AUTO_TEST_CASE(consistency_collision_boxbox) {
  Box s1(20, 40, 50);
  Box s2(10, 10, 10);

  BVHModel<AABB> s1_aabb;
  BVHModel<AABB> s2_aabb;
  BVHModel<OBB> s1_obb;
  BVHModel<OBB> s2_obb;

  generateBVHModel(s1_aabb, s1, Transform3s());
  generateBVHModel(s2_aabb, s2, Transform3s());
  generateBVHModel(s1_obb, s1, Transform3s());
  generateBVHModel(s2_obb, s2, Transform3s());

  CollisionRequest request(false, 1, false);
  CollisionResult result;

  bool res;

  Transform3s pose, pose_aabb, pose_obb;

  // s2 is within s1
  // both are shapes --> collision
  // s1 is shape, s2 is mesh --> in collision
  // s1 is mesh, s2 is shape --> collision free
  // all are reasonable
  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);

  pose.setTranslation(Vec3s(15.01, 0, 0));
  pose_aabb.setTranslation(Vec3s(15.01, 0, 0));
  pose_obb.setTranslation(Vec3s(15.01, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);

  pose.setTranslation(Vec3s(14.99, 0, 0));
  pose_aabb.setTranslation(Vec3s(14.99, 0, 0));
  pose_obb.setTranslation(Vec3s(14.99, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
}

BOOST_AUTO_TEST_CASE(consistency_collision_spherebox) {
  Sphere s1(20);
  Box s2(5, 5, 5);

  BVHModel<AABB> s1_aabb;
  BVHModel<AABB> s2_aabb;
  BVHModel<OBB> s1_obb;
  BVHModel<OBB> s2_obb;

  generateBVHModel(s1_aabb, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_aabb, s2, Transform3s());
  generateBVHModel(s1_obb, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_obb, s2, Transform3s());

  CollisionRequest request(false, 1, false);
  CollisionResult result;

  bool res;

  Transform3s pose, pose_aabb, pose_obb;

  // s2 is within s1
  // both are shapes --> collision
  // s1 is shape, s2 is mesh --> in collision
  // s1 is mesh, s2 is shape --> collision free
  // all are reasonable
  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);

  pose.setTranslation(Vec3s(22.4, 0, 0));
  pose_aabb.setTranslation(Vec3s(22.4, 0, 0));
  pose_obb.setTranslation(Vec3s(22.4, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);

  pose.setTranslation(Vec3s(22.51, 0, 0));
  pose_aabb.setTranslation(Vec3s(22.51, 0, 0));
  pose_obb.setTranslation(Vec3s(22.51, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
}

BOOST_AUTO_TEST_CASE(consistency_collision_cylindercylinder) {
  Cylinder s1(5, 10);
  Cylinder s2(5, 10);

  BVHModel<AABB> s1_aabb;
  BVHModel<AABB> s2_aabb;
  BVHModel<OBB> s1_obb;
  BVHModel<OBB> s2_obb;

  generateBVHModel(s1_aabb, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_aabb, s2, Transform3s(), 16, 16);
  generateBVHModel(s1_obb, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_obb, s2, Transform3s(), 16, 16);

  CollisionRequest request(false, 1, false);
  CollisionResult result;

  bool res;

  Transform3s pose, pose_aabb, pose_obb;

  pose.setTranslation(Vec3s(9.99, 0, 0));
  pose_aabb.setTranslation(Vec3s(9.99, 0, 0));
  pose_obb.setTranslation(Vec3s(9.99, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);

  pose.setTranslation(Vec3s(10.01, 0, 0));
  pose_aabb.setTranslation(Vec3s(10.01, 0, 0));
  pose_obb.setTranslation(Vec3s(10.01, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
}

BOOST_AUTO_TEST_CASE(consistency_collision_conecone) {
  Cone s1(5, 10);
  Cone s2(5, 10);

  BVHModel<AABB> s1_aabb;
  BVHModel<AABB> s2_aabb;
  BVHModel<OBB> s1_obb;
  BVHModel<OBB> s2_obb;

  generateBVHModel(s1_aabb, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_aabb, s2, Transform3s(), 16, 16);
  generateBVHModel(s1_obb, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_obb, s2, Transform3s(), 16, 16);

  CollisionRequest request(false, 1, false);
  CollisionResult result;

  bool res;

  Transform3s pose, pose_aabb, pose_obb;

  pose.setTranslation(Vec3s(9.9, 0, 0));
  pose_aabb.setTranslation(Vec3s(9.9, 0, 0));
  pose_obb.setTranslation(Vec3s(9.9, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);

  pose.setTranslation(Vec3s(10.1, 0, 0));
  pose_aabb.setTranslation(Vec3s(10.1, 0, 0));
  pose_obb.setTranslation(Vec3s(10.1, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);

  pose.setTranslation(Vec3s(0, 0, 9.9));
  pose_aabb.setTranslation(Vec3s(0, 0, 9.9));
  pose_obb.setTranslation(Vec3s(0, 0, 9.9));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);

  pose.setTranslation(Vec3s(0, 0, 10.1));
  pose_aabb.setTranslation(Vec3s(0, 0, 10.1));
  pose_obb.setTranslation(Vec3s(0, 0, 10.1));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
}

BOOST_AUTO_TEST_CASE(consistency_collision_spheresphere_GJK) {
  Sphere s1(20);
  Sphere s2(10);
  BVHModel<AABB> s1_aabb;
  BVHModel<AABB> s2_aabb;
  BVHModel<OBB> s1_obb;
  BVHModel<OBB> s2_obb;

  generateBVHModel(s1_aabb, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_aabb, s2, Transform3s(), 16, 16);
  generateBVHModel(s1_obb, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_obb, s2, Transform3s(), 16, 16);

  CollisionRequest request(false, 1, false);

  CollisionResult result;

  bool res;

  Transform3s pose, pose_aabb, pose_obb;

  // s2 is within s1
  // both are shapes --> collision
  // s1 is shape, s2 is mesh --> in collision
  // s1 is mesh, s2 is shape --> collision free
  // all are reasonable
  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);

  pose.setTranslation(Vec3s(40, 0, 0));
  pose_aabb.setTranslation(Vec3s(40, 0, 0));
  pose_obb.setTranslation(Vec3s(40, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);

  pose.setTranslation(Vec3s(30, 0, 0));
  pose_aabb.setTranslation(Vec3s(30, 0, 0));
  pose_obb.setTranslation(Vec3s(30, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);

  pose.setTranslation(Vec3s(29.9, 0, 0));
  pose_aabb.setTranslation(
      Vec3s(29.8, 0, 0));  // 29.9 fails, result depends on mesh precision
  pose_obb.setTranslation(
      Vec3s(29.8, 0, 0));  // 29.9 fails, result depends on mesh precision

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);

  pose.setTranslation(Vec3s(-29.9, 0, 0));
  pose_aabb.setTranslation(
      Vec3s(-29.8, 0, 0));  // 29.9 fails, result depends on mesh precision
  pose_obb.setTranslation(
      Vec3s(-29.8, 0, 0));  // 29.9 fails, result depends on mesh precision

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);

  pose.setTranslation(Vec3s(-30, 0, 0));
  pose_aabb.setTranslation(Vec3s(-30, 0, 0));
  pose_obb.setTranslation(Vec3s(-30, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
}

BOOST_AUTO_TEST_CASE(consistency_collision_boxbox_GJK) {
  Box s1(20, 40, 50);
  Box s2(10, 10, 10);

  BVHModel<AABB> s1_aabb;
  BVHModel<AABB> s2_aabb;
  BVHModel<OBB> s1_obb;
  BVHModel<OBB> s2_obb;

  generateBVHModel(s1_aabb, s1, Transform3s());
  generateBVHModel(s2_aabb, s2, Transform3s());
  generateBVHModel(s1_obb, s1, Transform3s());
  generateBVHModel(s2_obb, s2, Transform3s());

  CollisionRequest request(false, 1, false);

  CollisionResult result;

  bool res;

  Transform3s pose, pose_aabb, pose_obb;

  // s2 is within s1
  // both are shapes --> collision
  // s1 is shape, s2 is mesh --> in collision
  // s1 is mesh, s2 is shape --> collision free
  // all are reasonable
  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);

  pose.setTranslation(Vec3s(15.01, 0, 0));
  pose_aabb.setTranslation(Vec3s(15.01, 0, 0));
  pose_obb.setTranslation(Vec3s(15.01, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);

  pose.setTranslation(Vec3s(14.99, 0, 0));
  pose_aabb.setTranslation(Vec3s(14.99, 0, 0));
  pose_obb.setTranslation(Vec3s(14.99, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
}

BOOST_AUTO_TEST_CASE(consistency_collision_spherebox_GJK) {
  Sphere s1(20);
  Box s2(5, 5, 5);

  BVHModel<AABB> s1_aabb;
  BVHModel<AABB> s2_aabb;
  BVHModel<OBB> s1_obb;
  BVHModel<OBB> s2_obb;

  generateBVHModel(s1_aabb, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_aabb, s2, Transform3s());
  generateBVHModel(s1_obb, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_obb, s2, Transform3s());

  CollisionRequest request(false, 1, false);

  CollisionResult result;

  bool res;

  Transform3s pose, pose_aabb, pose_obb;

  // s2 is within s1
  // both are shapes --> collision
  // s1 is shape, s2 is mesh --> in collision
  // s1 is mesh, s2 is shape --> collision free
  // all are reasonable
  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);

  pose.setTranslation(Vec3s(22.4, 0, 0));
  pose_aabb.setTranslation(Vec3s(22.4, 0, 0));
  pose_obb.setTranslation(Vec3s(22.4, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);

  pose.setTranslation(Vec3s(22.51, 0, 0));
  pose_aabb.setTranslation(Vec3s(22.51, 0, 0));
  pose_obb.setTranslation(Vec3s(22.51, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
}

BOOST_AUTO_TEST_CASE(consistency_collision_cylindercylinder_GJK) {
  Cylinder s1(5, 10);
  Cylinder s2(5, 10);

  BVHModel<AABB> s1_aabb;
  BVHModel<AABB> s2_aabb;
  BVHModel<OBB> s1_obb;
  BVHModel<OBB> s2_obb;

  generateBVHModel(s1_aabb, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_aabb, s2, Transform3s(), 16, 16);
  generateBVHModel(s1_obb, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_obb, s2, Transform3s(), 16, 16);

  CollisionRequest request(false, 1, false);

  CollisionResult result;

  bool res;

  Transform3s pose, pose_aabb, pose_obb;

  pose.setTranslation(Vec3s(9.99, 0, 0));
  pose_aabb.setTranslation(Vec3s(9.99, 0, 0));
  pose_obb.setTranslation(Vec3s(9.99, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);

  pose.setTranslation(Vec3s(10.01, 0, 0));
  pose_aabb.setTranslation(Vec3s(10.01, 0, 0));
  pose_obb.setTranslation(Vec3s(10.01, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
}

BOOST_AUTO_TEST_CASE(consistency_collision_conecone_GJK) {
  Cone s1(5, 10);
  Cone s2(5, 10);

  BVHModel<AABB> s1_aabb;
  BVHModel<AABB> s2_aabb;
  BVHModel<OBB> s1_obb;
  BVHModel<OBB> s2_obb;

  generateBVHModel(s1_aabb, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_aabb, s2, Transform3s(), 16, 16);
  generateBVHModel(s1_obb, s1, Transform3s(), 16, 16);
  generateBVHModel(s2_obb, s2, Transform3s(), 16, 16);

  CollisionRequest request(false, 1, false);

  CollisionResult result;

  bool res;

  Transform3s pose, pose_aabb, pose_obb;

  pose.setTranslation(Vec3s(9.9, 0, 0));
  pose_aabb.setTranslation(Vec3s(9.9, 0, 0));
  pose_obb.setTranslation(Vec3s(9.9, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);

  pose.setTranslation(Vec3s(10.1, 0, 0));
  pose_aabb.setTranslation(Vec3s(10.1, 0, 0));
  pose_obb.setTranslation(Vec3s(10.1, 0, 0));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);

  pose.setTranslation(Vec3s(0, 0, 9.9));
  pose_aabb.setTranslation(Vec3s(0, 0, 9.9));
  pose_obb.setTranslation(Vec3s(0, 0, 9.9));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK(res);

  pose.setTranslation(Vec3s(0, 0, 10.1));
  pose_aabb.setTranslation(Vec3s(0, 0, 10.1));
  pose_obb.setTranslation(Vec3s(0, 0, 10.1));

  result.clear();
  res = (collide(&s1, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_aabb, pose_aabb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2_obb, pose_obb, &s1, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_aabb, pose_aabb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1, Transform3s(), &s2_obb, pose_obb, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_aabb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s2, pose, &s1_obb, Transform3s(), request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
  result.clear();
  res = (collide(&s1_aabb, Transform3s(), &s2, pose, request, result) > 0);
  BOOST_CHECK_FALSE(res);
}
