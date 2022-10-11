//
// Software License Agreement (BSD License)
//
//  Copyright (c) 2020 CNRS-LAAS
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
#include <hpp/fcl/narrowphase/gjk.h>
#include <hpp/fcl/internal/intersect.h>

#ifdef HPP_FCL_HAS_DOXYGEN_AUTODOC
#include "doxygen_autodoc/functions.h"
#include "doxygen_autodoc/hpp/fcl/narrowphase/gjk.h"
#endif

using namespace boost::python;
using namespace hpp::fcl;
using hpp::fcl::details::EPA;
using hpp::fcl::details::GJK;
using hpp::fcl::details::MinkowskiDiff;

struct SupportData
{
  public:
    Vec3f support;
    int hint;
    SupportData(){
      support.setZero();
      hint = 0;
    }
};

void getSupport(const ShapeBase* shape, const Vec3f& dir, bool is_normalized, SupportData& support_data)
{
  int hint = support_data.hint;
  Vec3f support = hpp::fcl::details::getSupport(shape, dir, is_normalized, hint);
  support_data.hint = hint;
  support_data.support = support;
}

struct SimplexWrapper 
{
  hpp::fcl::details::GJK::SimplexV vertices[4];
  hpp::fcl::details::GJK::vertex_id_t rank;
  inline void copy(const SimplexWrapper& other) {
    rank = other.rank;
    for (int i = 0; i < rank; i++) {
      vertices[i].index_w0 = other.vertices[i].index_w0;
      vertices[i].index_w1 = other.vertices[i].index_w1;
      vertices[i].w = other.vertices[i].w;
      vertices[i].w0 = other.vertices[i].w0;
      vertices[i].w1 = other.vertices[i].w1;
    }
  }

  hpp::fcl::details::GJK::SimplexV getVertex(int i){
    return vertices[i];
  }

  void setVertex(const hpp::fcl::details::GJK::SimplexV& vertex, int i){
    vertices[i] = vertex;
  }
};

Vec3f projectOriginOntoSimplex(const SimplexWrapper& curr_simplex, SimplexWrapper& next_simplex)
{
  Project::ProjectResult projection_result;
  // Project simplex
  switch (curr_simplex.rank) {
    case 1:
      next_simplex.copy(curr_simplex);  
      break;
    case 2:
      projection_result = Project::projectLineOrigin(curr_simplex.vertices[0].w,
                                                     curr_simplex.vertices[1].w);
      break;
    case 3:
      projection_result = Project::projectTriangleOrigin(curr_simplex.vertices[0].w,
                                                         curr_simplex.vertices[1].w,
                                                         curr_simplex.vertices[2].w);
      break;
    case 4:
      projection_result = Project::projectTetrahedraOrigin(curr_simplex.vertices[0].w,
                                                           curr_simplex.vertices[1].w,
                                                           curr_simplex.vertices[2].w,
                                                           curr_simplex.vertices[3].w);
      break;
    default:
      throw std::logic_error("Invalid simplex rank!");
  }
  if (curr_simplex.rank == 1){
    return curr_simplex.vertices[0].w;
  }

  // Construct next_simplex and projection
  Vec3f projection = Vec3f::Zero();
  next_simplex.rank = 0;
  for (int i = 0; i < curr_simplex.rank; i++) {
    projection += projection_result.parameterization[i] * curr_simplex.vertices[i].w;
    if (projection_result.parameterization[i] > 1e-12){
      next_simplex.vertices[next_simplex.rank] = curr_simplex.vertices[i];
      next_simplex.rank++;
    }
  }
  return projection;
}

void exposeGJK() {
  if (!eigenpy::register_symbolic_link_to_registered_type<GJK::Status>()) {
    enum_<GJK::Status>("GJKStatus")
        .value("Valid", GJK::Valid)
        .value("Inside", GJK::Inside)
        .value("Failed", GJK::Failed)
        .export_values();
  }

  if (!eigenpy::register_symbolic_link_to_registered_type<GJKVariant>()) {
    enum_<GJKVariant>("GJKVariant")
        .value("DefaultGJK", GJKVariant::DefaultGJK)
        .value("NesterovAcceleration", GJKVariant::NesterovAcceleration)
        .value("PolyakAcceleration", GJKVariant::PolyakAcceleration)
        .export_values();
  }

  if (!eigenpy::register_symbolic_link_to_registered_type<GJKInitialGuess>()) {
    enum_<GJKInitialGuess>("GJKInitialGuess")
        .value("DefaultGuess", GJKInitialGuess::DefaultGuess)
        .value("CachedGuess", GJKInitialGuess::CachedGuess)
        .value("BoundingVolumeGuess", GJKInitialGuess::BoundingVolumeGuess)
        .export_values();
  }

  if (!eigenpy::register_symbolic_link_to_registered_type<
          GJKConvergenceCriterion>()) {
    enum_<GJKConvergenceCriterion>("GJKConvergenceCriterion")
        .value("VDB", GJKConvergenceCriterion::VDB)
        .value("DualityGap", GJKConvergenceCriterion::DualityGap)
        .value("Hybrid", GJKConvergenceCriterion::Hybrid)
        .export_values();
  }

  if (!eigenpy::register_symbolic_link_to_registered_type<
          GJKConvergenceCriterionType>()) {
    enum_<GJKConvergenceCriterionType>("GJKConvergenceCriterionType")
        .value("Absolute", GJKConvergenceCriterionType::Absolute)
        .value("Relative", GJKConvergenceCriterionType::Relative)
        .export_values();
  }

  if (!eigenpy::register_symbolic_link_to_registered_type<SupportData>()) {
    class_<SupportData>("SupportData", doxygen::class_doc<SupportData>(),
                          no_init)
      .def(doxygen::visitor::init<SupportData>())
        .DEF_RW_CLASS_ATTRIB(SupportData, support)
        .DEF_RW_CLASS_ATTRIB(SupportData, hint);
  }

  def("getSupport", &getSupport);

  if (!eigenpy::register_symbolic_link_to_registered_type<GJK::SimplexV>()) {
    class_<GJK::SimplexV>("SimplexV", doxygen::class_doc<GJK::SimplexV>(),
                          no_init)
        .def(doxygen::visitor::init<GJK::SimplexV>())
        .DEF_RW_CLASS_ATTRIB(GJK::SimplexV, w0)
        .DEF_RW_CLASS_ATTRIB(GJK::SimplexV, w1)
        .DEF_RW_CLASS_ATTRIB(GJK::SimplexV, w);
  }

  if (!eigenpy::register_symbolic_link_to_registered_type<SimplexWrapper>()) {
    class_<SimplexWrapper>("Simplex", doxygen::class_doc<SimplexWrapper>(),
                          no_init)
        .def(doxygen::visitor::init<SimplexWrapper>())
        .def("getVertex", &SimplexWrapper::getVertex)
        .def("setVertex", &SimplexWrapper::setVertex)
        .def("copy", &SimplexWrapper::copy)
        .DEF_RW_CLASS_ATTRIB(SimplexWrapper, rank);
  }

  def("projectOriginOntoSimplex", &projectOriginOntoSimplex);

  if (!eigenpy::register_symbolic_link_to_registered_type<MinkowskiDiff>()) {
    class_<MinkowskiDiff>("MinkowskiDiff", doxygen::class_doc<MinkowskiDiff>(),
                          no_init)
        .def(doxygen::visitor::init<MinkowskiDiff>())
        .def("set",
             static_cast<void (MinkowskiDiff::*)(
                 const ShapeBase*, const ShapeBase*)>(&MinkowskiDiff::set),
             doxygen::member_func_doc(
                 static_cast<void (MinkowskiDiff::*)(
                     const ShapeBase*, const ShapeBase*)>(&MinkowskiDiff::set)))
        .def("set",
             static_cast<void (MinkowskiDiff::*)(
                 const ShapeBase*, const ShapeBase*, const Transform3f& tf0,
                 const Transform3f& tf1)>(&MinkowskiDiff::set),
             doxygen::member_func_doc(
                 static_cast<void (MinkowskiDiff::*)(
                     const ShapeBase*, const ShapeBase*, const Transform3f& tf0,
                     const Transform3f& tf1)>(&MinkowskiDiff::set)))
        .DEF_RW_CLASS_ATTRIB(MinkowskiDiff, inflation)
        .DEF_RW_CLASS_ATTRIB(MinkowskiDiff, normalize_support_direction);
  }

  if (!eigenpy::register_symbolic_link_to_registered_type<GJK>()) {
    class_<GJK>("GJK", doxygen::class_doc<GJK>(), no_init)
        .def(doxygen::visitor::init<GJK, unsigned int, FCL_REAL>())
        .DEF_RW_CLASS_ATTRIB(GJK, distance)
        .DEF_RW_CLASS_ATTRIB(GJK, ray)
        .DEF_RW_CLASS_ATTRIB(GJK, x0)
        .DEF_RW_CLASS_ATTRIB(GJK, x1)
        .DEF_RW_CLASS_ATTRIB(GJK, support_hint)
        .DEF_RW_CLASS_ATTRIB(GJK, nfree)
        .DEF_RW_CLASS_ATTRIB(GJK, gjk_variant)
        .DEF_RW_CLASS_ATTRIB(GJK, convergence_criterion)
        .DEF_RW_CLASS_ATTRIB(GJK, convergence_criterion_type)
        .DEF_CLASS_FUNC(GJK, evaluate)
        .DEF_CLASS_FUNC(GJK, hasClosestPoints)
        .DEF_CLASS_FUNC(GJK, hasPenetrationInformation)
        .DEF_CLASS_FUNC(GJK, getClosestPoints)
        .DEF_CLASS_FUNC(GJK, computeClosestPoints)
        .DEF_CLASS_FUNC(GJK, setDistanceEarlyBreak)
        .DEF_CLASS_FUNC(GJK, getGuessFromSimplex)
        .DEF_CLASS_FUNC(GJK, getIterations)
        .DEF_CLASS_FUNC(GJK, getIterationsEarly)
        .DEF_CLASS_FUNC(GJK, getNumCallSupport)
        .DEF_CLASS_FUNC(GJK, getNumCallSupportEarly)
        .DEF_CLASS_FUNC(GJK, getCumulativeSupportDotprods)
        .DEF_CLASS_FUNC(GJK, getCumulativeSupportDotprodsEarly)
        .DEF_CLASS_FUNC(GJK, getNumCallProjection)
        .DEF_CLASS_FUNC(GJK, getNumCallProjectionEarly)
        .DEF_CLASS_FUNC(GJK, measureRunTime)
        .DEF_CLASS_FUNC(GJK, getGJKRunTime)
        .DEF_CLASS_FUNC(GJK, getGJKRunTimeEarly)
        .DEF_CLASS_FUNC(GJK, computeGJKAverageRunTime)
        .DEF_CLASS_FUNC(GJK, getAverageGJKRunTime)
        .DEF_CLASS_FUNC(GJK, getAverageGJKRunTimeEarly);
  }

  def("projectLineOrigin", &Project::projectLineOrigin);
  def("projectTriangleOrigin", &Project::projectTriangleOrigin);
  def("projectTetrahedraOrigin", &Project::projectTetrahedraOrigin);

  if (!eigenpy::register_symbolic_link_to_registered_type<
          Project::ProjectResult>()) {
    class_<Project::ProjectResult>(
        "ProjectResult", doxygen::class_doc<Project::ProjectResult>(), no_init)
        .def(doxygen::visitor::init<Project::ProjectResult>())
        .DEF_RW_CLASS_ATTRIB(Project::ProjectResult, parameterization_eigen)
        .DEF_CLASS_FUNC(Project::ProjectResult, updateParameterization);
  }
}
