#define BOOST_TEST_MODULE FCL_EPA_DEBUG
#include <boost/test/included/unit_test.hpp>
#include <hpp/fcl/shape/geometric_shapes.h>
#include "hpp/fcl/narrowphase/gjk.h"
#include "hpp/fcl/collision_data.h"
#include "hpp/fcl/distance.h"

using hpp::fcl::ShapeBase;
using hpp::fcl::Ellipsoid;
using hpp::fcl::Vec3f;
using hpp::fcl::Transform3f;
using hpp::fcl::details::MinkowskiDiff;
using hpp::fcl::details::GJK;
using hpp::fcl::details::EPA;
// using hpp::fcl::DistanceRequest;
// using hpp::fcl::DistanceResult;
using hpp::fcl::FCL_REAL;

void run_gjk_epa(const ShapeBase* shape1, const Transform3f& tf1,
                 const ShapeBase* shape2, const Transform3f& tf2)
{

  MinkowskiDiff mink;
  mink.set(shape1, shape2, tf1, tf2);

  GJK gjk(100, 1e-8);
  gjk.convergence_criterion = hpp::fcl::GJKConvergenceCriterion::DualityGap;
  gjk.convergence_criterion_type = hpp::fcl::GJKConvergenceCriterionType::Absolute;
  Vec3f guess(1, 1, 1);
  hpp::fcl::support_func_guess_t hint(0, 0);
  GJK::Status gjk_status = gjk.evaluate(mink, guess, hint);
  std::cout << "GJK status: " << gjk_status << std::endl;
  std::cout << "GJK distance: " << (gjk.ray).norm() << std::endl;
  std::cout << "GJK numit: " << gjk.getIterations() << std::endl;
  if (gjk.getSimplex())
    std::cout << "GJK final simplex rank: " << (int)(gjk.getSimplex()->rank) << std::endl;
  if (gjk.getSimplex()->rank >= 2){
    const Vec3f& a = gjk.getSimplex()->vertex[0]->w;
    const Vec3f& b = gjk.getSimplex()->vertex[1]->w;
    // std::cout << "a: " << a.transpose() << std::endl;
    // std::cout << "b: " << b.transpose() << std::endl;
    // std::cout << "(a - b).norm(): " << (a-b).norm() << std::endl;
    BOOST_CHECK((a - b).squaredNorm() >= gjk.getTolerance());
  }

  EPA epa(128, 64, 100, gjk.getTolerance());
  EPA::Status epa_status;
  double gjk_time = 0.;
  double epa_time = 0.;
  hpp::fcl::Timer timer;
  for (int i = 0; i < 100; i++) {
    timer.stop(); timer.start();
    gjk_status = gjk.evaluate(mink, guess, hint);
    timer.stop();
    gjk_time += timer.elapsed().user;
    timer.stop(); timer.start();
    epa_status = epa.evaluate(gjk, guess);
    timer.stop();
    epa_time += timer.elapsed().user;
  }
  std::cout << "GJK time: " << gjk_time / 100 << " us"<< std::endl;
  std::cout << "EPA status: " << epa_status << std::endl;
  std::cout << "EPA depth: " << epa.depth << std::endl;
  std::cout << "EPA numit: " << epa.getIterations() << std::endl;
  std::cout << "EPA normal: " << epa.normal.transpose() << std::endl;
  std::cout << "EPA time: " << epa_time / 100 << " us" << std::endl;
  std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(epa_debug_sphere) {
  double r = 0.5;
  Vec3f r1(r, r, r);
  Ellipsoid shape1(r1);
  Vec3f r2(0.5, 0.5, 0.5);
  Ellipsoid shape2(r2);

  Transform3f tf1 = Transform3f::Identity();
  Transform3f tf2 = Transform3f::Identity();
  tf2.setTranslation(Vec3f(r + 0.4, 0., 0.));

  std::cout << "--- SPHERES ---" << std::endl;
  run_gjk_epa(&shape1, tf1, &shape2, tf2);

  // std::cout << std::endl;

  // DistanceRequest req;
  // req.gjk_convergence_criterion = hpp::fcl::GJKConvergenceCriterion::DualityGap;
  // req.gjk_convergence_criterion_type = hpp::fcl::GJKConvergenceCriterionType::Absolute;
  // req.gjk_tolerance = 1e-8;
  // DistanceResult res;
  // FCL_REAL dist = hpp::fcl::distance(&shape1, tf1, &shape2, tf2, req, res);
  // std::cout << "DISTANCE: " << dist << std::endl;
  // std::cout << "EPA NUMIT: " << res.epa_numit << std::endl;
}

#include "utility.h"
using hpp::fcl::Triangle;
using hpp::fcl::Convex;

BOOST_AUTO_TEST_CASE(epa_debug_polyhedral_ellipsoids) {
  Ellipsoid ellipsoid1(0.5, 0.5, 0.5);
  Ellipsoid ellipsoid2(0.5, 0.5, 0.5);
  Convex<Triangle> shape1 = hpp::fcl::constructPolytopeFromEllipsoid(ellipsoid1);
  Convex<Triangle> shape2 = hpp::fcl::constructPolytopeFromEllipsoid(ellipsoid2);

  Transform3f tf1 = Transform3f::Identity();
  Transform3f tf2 = Transform3f::Identity();
  tf2.setTranslation(Vec3f(0.5, 0.5, 0.));

  std::cout << "--- POLYHEDRAL ELLIPSOIDS ---" << std::endl;
  run_gjk_epa(&shape1, tf1, &shape2, tf2);
}

Convex<Triangle> buildDiamond(){
  Vec3f* pts = new Vec3f[6];
  pts[0] = Vec3f(1, 0, 0);
  pts[1] = Vec3f(0, 1, 0);
  pts[2] = Vec3f(-1, 0, 0);
  pts[3] = Vec3f(0, -1, 0);
  pts[4] = Vec3f(0, 0, 0.25);
  pts[5] = Vec3f(0, 0, -0.25);

  Triangle* tris = new Triangle[8];
  tris[0].set(0, 1, 4);
  tris[1].set(1, 2, 4);
  tris[2].set(3, 0, 4);
  tris[3].set(0, 5, 1);
  tris[4].set(1, 5, 2);
  tris[5].set(2, 5, 3);
  tris[6].set(2, 5, 3);
  tris[7].set(3, 5, 0);

  return Convex<Triangle>(true,
                           pts,  // points
                           6,    // num points
                           tris,
                           8  // number of polygons
  );
}

BOOST_AUTO_TEST_CASE(epa_debug_diamond) {
  Convex<Triangle> shape1 = buildDiamond();
  Convex<Triangle> shape2 = buildDiamond();

  Transform3f tf1 = Transform3f::Identity();
  Transform3f tf2 = Transform3f::Identity();
  tf2.setTranslation(Vec3f(1.5, 0., 0.));

  std::cout << "--- DIAMONDS ---" << std::endl;
  run_gjk_epa(&shape1, tf1, &shape2, tf2);
}

