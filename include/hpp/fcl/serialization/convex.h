//
// Copyright (c) 2022 INRIA
//

#ifndef HPP_FCL_SERIALIZATION_CONVEX_H
#define HPP_FCL_SERIALIZATION_CONVEX_H

#include "hpp/fcl/shape/convex.h"

#include "hpp/fcl/serialization/fwd.h"
#include "hpp/fcl/serialization/geometric_shapes.h"
#include "hpp/fcl/serialization/memory.h"
#include "hpp/fcl/serialization/triangle.h"
#include "hpp/fcl/serialization/quadrilateral.h"

namespace boost {
namespace serialization {

namespace internal {
struct ConvexBaseAccessor : hpp::fcl::ConvexBase {
  typedef hpp::fcl::ConvexBase Base;
  using Base::nneighbors_;
  using Base::own_storage_;
};

}  // namespace internal

template <class Archive>
void serialize(Archive &ar, hpp::fcl::ConvexBase &convex_base,
               const unsigned int /*version*/) {
  using namespace hpp::fcl;

  typedef internal::ConvexBaseAccessor Accessor;
  Accessor &accessor = reinterpret_cast<Accessor &>(convex_base);

  ar &make_nvp("base", boost::serialization::base_object<hpp::fcl::ShapeBase>(
                           convex_base));
  const unsigned int num_points_previous = convex_base.num_points;
  ar &make_nvp("num_points", convex_base.num_points);

  if (Archive::is_loading::value) {
    if (num_points_previous != convex_base.num_points ||
        !accessor.own_storage_) {
      delete[] convex_base.points;
      convex_base.points = new hpp::fcl::Vec3f[convex_base.num_points];
      accessor.own_storage_ = true;
    }
  }

  {
    typedef Eigen::Matrix<FCL_REAL, 3, Eigen::Dynamic> MatrixPoints;
    Eigen::Map<MatrixPoints> points_map(
        reinterpret_cast<double *>(convex_base.points), 3,
        convex_base.num_points);
    ar &make_nvp("points", points_map);
  }

  ar &make_nvp("center", convex_base.center);
  // We don't save neighbors as they will be computed directly by calling
  // fillNeighbors.
  // Update: actually, it's important to save neighbors if we save a ConvexBase
  // and not a Convex<PolygonT>.
  // For example the GJK algorithm requires neighbors for the support function
  // of large ConvexBase objects.
  {
    std::vector<unsigned int> neighbors_indexes;
    std::vector<unsigned char> num_neighbors;
    if (Archive::is_loading::value) {
      ar &make_nvp("neighbors_indexes", neighbors_indexes);
      ar &make_nvp("num_neighbors", num_neighbors);
      assert(num_neighbors.size() == convex_base.num_points);
      convex_base.neighbors = new ConvexBase::Neighbors[convex_base.num_points];
      accessor.nneighbors_ = new unsigned int[neighbors_indexes.size()];
      std::copy(neighbors_indexes.begin(), neighbors_indexes.end(),
                accessor.nneighbors_);

      unsigned int *p_nneighbors = accessor.nneighbors_;
      for (size_t i = 0; i < static_cast<size_t>(convex_base.num_points); ++i) {
        ConvexBase::Neighbors &n = convex_base.neighbors[i];
        n.count_ = num_neighbors[i];
        n.n_ = p_nneighbors;
        p_nneighbors += num_neighbors[i];
      }
    } else {
      neighbors_indexes.clear();
      num_neighbors.clear();
      num_neighbors.reserve(convex_base.num_points);
      for (size_t i = 0; i < static_cast<size_t>(convex_base.num_points); ++i) {
        ConvexBase::Neighbors &n = convex_base.neighbors[i];
        num_neighbors.push_back(n.count_);
        for (size_t j = 0; j < static_cast<size_t>(n.count_); ++j) {
          neighbors_indexes.push_back(n.n_[j]);
        }
      }
      ar &make_nvp("neighbors_indexes", neighbors_indexes);
      ar &make_nvp("num_neighbors", num_neighbors);
    }
  }
}

namespace internal {
template <typename PolygonT>
struct ConvexAccessor : hpp::fcl::Convex<PolygonT> {
  typedef hpp::fcl::Convex<PolygonT> Base;
  using Base::fillNeighbors;
};

}  // namespace internal

template <class Archive, typename PolygonT>
void serialize(Archive &ar, hpp::fcl::Convex<PolygonT> &convex_,
               const unsigned int /*version*/) {
  using namespace hpp::fcl;
  typedef internal::ConvexAccessor<PolygonT> Accessor;

  Accessor &convex = reinterpret_cast<Accessor &>(convex_);
  ar &make_nvp("base", boost::serialization::base_object<ConvexBase>(convex));

  const unsigned int num_polygons_previous = convex.num_polygons;
  ar &make_nvp("num_polygons", convex.num_polygons);

  if (Archive::is_loading::value) {
    if (num_polygons_previous != convex.num_polygons) {
      delete[] convex.polygons;
      convex.polygons = new PolygonT[convex.num_polygons];
    }
  }

  ar &make_array<PolygonT>(convex.polygons, convex.num_polygons);

  if (Archive::is_loading::value) convex.fillNeighbors();
}

}  // namespace serialization
}  // namespace boost

namespace hpp {
namespace fcl {

// namespace internal {
// template <typename BV>
// struct memory_footprint_evaluator< ::hpp::fcl::BVHModel<BV> > {
//   static size_t run(const ::hpp::fcl::BVHModel<BV> &bvh_model) {
//     return static_cast<size_t>(bvh_model.memUsage(false));
//   }
// };
// }  // namespace internal

}  // namespace fcl
}  // namespace hpp

#endif  // ifndef HPP_FCL_SERIALIZATION_CONVEX_H
