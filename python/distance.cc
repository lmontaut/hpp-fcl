//
// Software License Agreement (BSD License)
//
//  Copyright (c) 2019 CNRS-LAAS INRIA
//  Author: Joseph Mirabel
//  All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions
//  are met:
//
//   * Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above
//     copyright notice, this list of conditions and the following
//     disclaimer in the documentation and/or other materials provided
//     with the distribution.
//   * Neither the name of CNRS-LAAS. nor the names of its
//     contributors may be used to endorse or promote products derived
//     from this software without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
//  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
//  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
//  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
//  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//  POSSIBILITY OF SUCH DAMAGE.

#include <eigenpy/eigenpy.hpp>

#include "fcl.hh"

#include <hpp/fcl/fwd.hh>
#include <hpp/fcl/distance.h>

#ifdef HPP_FCL_HAS_DOXYGEN_AUTODOC
#include "doxygen_autodoc/functions.h"
#include "doxygen_autodoc/hpp/fcl/collision_data.h"
#endif

using namespace boost::python;
using namespace hpp::fcl;
using hpp::fcl::details::GJK;

namespace dv = doxygen::visitor;

struct DistanceResultWrapper {
  static Vec3f getNearestPoint1(const DistanceResult& res) {
    return res.nearest_points[0];
  }
  static Vec3f getNearestPoint2(const DistanceResult& res) {
    return res.nearest_points[1];
  }

  static std::size_t getNumNeighbors1(const DistanceResult& res) {
    return res.nearest_points_neighbors[0].size();
  }

  static Vec3f getNeighborIndex1(const DistanceResult& res, const int& i) {
    return res.nearest_points_neighbors[0][static_cast<std::size_t>(i)];
  }

  static FCL_REAL getSoftmaxWeight1(const DistanceResult& res, const int& i) {
    return res.softmax_weights[0][i];
  }

  static Vec3f getNeighborIndex2(const DistanceResult& res, const int& i) {
    return res.nearest_points_neighbors[1][static_cast<std::size_t>(i)];
  }

  static std::size_t getNumNeighbors2(const DistanceResult& res) {
    return res.nearest_points_neighbors[1].size();
  }

  static FCL_REAL getSoftmaxWeight2(const DistanceResult& res, const int& i) {
    return res.softmax_weights[1][i];
  }
};

struct SimplexSupportWrapper {
  static int getIndexw0(const GJK::SimplexSupport& simplex_support, const int& i) {
    return simplex_support.index_w[0][static_cast<std::size_t>(i)];
  }
  static int getIndexw1(const GJK::SimplexSupport& simplex_support, const int& i) {
    return simplex_support.index_w[1][static_cast<std::size_t>(i)];
  }
  static int getRank(const GJK::SimplexSupport& simplex_support) {
    return simplex_support.rank;
  }
};

void exposeDistanceAPI() {
  if (!eigenpy::register_symbolic_link_to_registered_type<DistanceRequest>()) {
    class_<DistanceRequest, bases<QueryRequest> >(
        "DistanceRequest", doxygen::class_doc<DistanceRequest>(),
        init<optional<bool, FCL_REAL, FCL_REAL> >(
            (arg("self"), arg("enable_nearest_points"), arg("rel_err"),
             arg("abs_err")),
            "Constructor"))
        .DEF_RW_CLASS_ATTRIB(DistanceRequest, enable_nearest_points)
        .DEF_RW_CLASS_ATTRIB(DistanceRequest, derivative_type)
        .DEF_RW_CLASS_ATTRIB(DistanceRequest, derivative_options)
        .DEF_RW_CLASS_ATTRIB(DistanceRequest, max_neighbors_search_level)
        .DEF_RW_CLASS_ATTRIB(DistanceRequest, rel_err)
        .DEF_RW_CLASS_ATTRIB(DistanceRequest, abs_err);
  }

  if (!eigenpy::register_symbolic_link_to_registered_type<
          std::vector<DistanceRequest> >()) {
    class_<std::vector<DistanceRequest> >("StdVec_DistanceRequest")
        .def(vector_indexing_suite<std::vector<DistanceRequest> >());
  }

  if (!eigenpy::register_symbolic_link_to_registered_type<GJK::SimplexSupport>()) {
    class_<GJK::SimplexSupport>("SimplexSupport", doxygen::class_doc<GJK::SimplexSupport>(),
                          no_init)
        .def(dv::init<GJK::SimplexSupport>())
        .def("getIndexw0", &SimplexSupportWrapper::getIndexw0)
        .def("getIndexw1", &SimplexSupportWrapper::getIndexw1)
        .def("getRank", &SimplexSupportWrapper::getRank);
  }

  if (!eigenpy::register_symbolic_link_to_registered_type<DistanceResult>()) {
    class_<DistanceResult, bases<QueryResult> >(
        "DistanceResult", doxygen::class_doc<DistanceResult>(), no_init)
        .def(dv::init<DistanceResult>())
        .DEF_RW_CLASS_ATTRIB(DistanceResult, min_distance)
        .DEF_RW_CLASS_ATTRIB(DistanceResult, normal)
        //.def_readwrite ("nearest_points", &DistanceResult::nearest_points)
        .def("getNearestPoint1", &DistanceResultWrapper::getNearestPoint1,
             doxygen::class_attrib_doc<DistanceResult>("nearest_points"))
        .def("getNearestPoint2", &DistanceResultWrapper::getNearestPoint2,
             doxygen::class_attrib_doc<DistanceResult>("nearest_points"))
        .def("getNumNeighbors1", &DistanceResultWrapper::getNumNeighbors1)
        .def("getNumNeighbors2", &DistanceResultWrapper::getNumNeighbors2)
        .def("getNeighborIndex1", &DistanceResultWrapper::getNeighborIndex1)
        .def("getNeighborIndex2", &DistanceResultWrapper::getNeighborIndex2)
        .def("getSoftmaxWeight1", &DistanceResultWrapper::getSoftmaxWeight1)
        .def("getSoftmaxWeight2", &DistanceResultWrapper::getSoftmaxWeight2)
        .DEF_RO_CLASS_ATTRIB(DistanceResult, nearest_points)
        .DEF_RO_CLASS_ATTRIB(DistanceResult, o1)
        .DEF_RO_CLASS_ATTRIB(DistanceResult, o2)
        .DEF_RO_CLASS_ATTRIB(DistanceResult, nearest_points_neighbors)
        .DEF_RO_CLASS_ATTRIB(DistanceResult, softmax_weights)
        .DEF_RO_CLASS_ATTRIB(DistanceResult, simplex_support)
        .DEF_RO_CLASS_ATTRIB(DistanceResult, time_distance_derivatives)
        .DEF_RO_CLASS_ATTRIB(DistanceResult, ds1_ddir1)
        .DEF_RO_CLASS_ATTRIB(DistanceResult, ds2_ddir2)
        .DEF_RW_CLASS_ATTRIB(DistanceResult, b1)
        .DEF_RW_CLASS_ATTRIB(DistanceResult, b2)
        .DEF_RW_CLASS_ATTRIB(DistanceResult, w)
        .DEF_RW_CLASS_ATTRIB(DistanceResult, w1)
        .DEF_RW_CLASS_ATTRIB(DistanceResult, w2)
        .DEF_RW_CLASS_ATTRIB(DistanceResult, dw_dq)
        .DEF_RW_CLASS_ATTRIB(DistanceResult, dw1_dq)
        .DEF_RW_CLASS_ATTRIB(DistanceResult, dw2_dq)
        .DEF_RW_CLASS_ATTRIB(DistanceResult, dw_dq1)
        .DEF_RW_CLASS_ATTRIB(DistanceResult, dw_dq2)
        .DEF_RW_CLASS_ATTRIB(DistanceResult, dw1_dq1)
        .DEF_RW_CLASS_ATTRIB(DistanceResult, dw1_dq2)
        .DEF_RW_CLASS_ATTRIB(DistanceResult, dw2_dq1)
        .DEF_RW_CLASS_ATTRIB(DistanceResult, dw2_dq2)

        .def("clear", &DistanceResult::clear,
             doxygen::member_func_doc(&DistanceResult::clear));
  }

  if (!eigenpy::register_symbolic_link_to_registered_type<DerivativeType>()) {
    enum_<DerivativeType>("DerivativeType")
        .value("FiniteDifference", DerivativeType::FiniteDifference)
        .value("ZeroOrderRS", DerivativeType::ZeroOrderRS)
        .value("FirstOrderRS", DerivativeType::FirstOrderRS)
        .value("FirstOrderGumbel", DerivativeType::FirstOrderGumbel)
        .export_values();
  }

  if (!eigenpy::register_symbolic_link_to_registered_type<DerivativeOptions>()) {
    class_<DerivativeOptions>(
        "DerivativeOptions", doxygen::class_doc<DerivativeOptions>(), no_init)
        .def(dv::init<DerivativeOptions, FCL_REAL, int, Vec3f, support_func_guess_t, bool>())
        .DEF_RW_CLASS_ATTRIB(DerivativeOptions, noise)
        .DEF_RW_CLASS_ATTRIB(DerivativeOptions, num_samples)
        .DEF_RW_CLASS_ATTRIB(DerivativeOptions, warm_start)
        .DEF_RW_CLASS_ATTRIB(DerivativeOptions, use_analytic_hessians)
        .DEF_RW_CLASS_ATTRIB(DerivativeOptions, hint);
  }

  if (!eigenpy::register_symbolic_link_to_registered_type<
          std::vector<DistanceResult> >()) {
    class_<std::vector<DistanceResult> >("StdVec_DistanceResult")
        .def(vector_indexing_suite<std::vector<DistanceResult> >());
  }

  doxygen::def(
      "distance",
      static_cast<FCL_REAL (*)(const CollisionObject*, const CollisionObject*,
                               const DistanceRequest&, DistanceResult&)>(
          &distance));
  doxygen::def(
      "distance",
      static_cast<FCL_REAL (*)(const CollisionGeometry*, const Transform3f&,
                               const CollisionGeometry*, const Transform3f&,
                               DistanceRequest&, DistanceResult&)>(&distance));

  class_<ComputeDistance>("ComputeDistance",
                          doxygen::class_doc<ComputeDistance>(), no_init)
      .def(dv::init<ComputeDistance, const CollisionGeometry*,
                    const CollisionGeometry*>())
      .def("__call__",
           static_cast<FCL_REAL (ComputeDistance::*)(
               const Transform3f&, const Transform3f&, DistanceRequest&,
               DistanceResult&) const>(&ComputeDistance::operator()));
}
