//===-- Test driver for OMEGA tendency terms ---------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA tendency term functors
///
/// This driver tests the functors used to calculate the tendencies used to
/// update OMEGA state variables. The tests are designed to be run with the
/// planar and spherical meshes described in the OMEGA Quick Start. For each
/// functor, input arrays are initialized based on arbitrary periodic functions
/// defined in the structs for the planar and spherical configurations. The
/// difference between analytical solutions and the output of each function
/// are used to calculate L2 and L-Infinity error norms, which are compared to
/// expected values for the given mesh.
///
//
//===-----------------------------------------------------------------------===/
#include "TendencyTerms.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Error.h"
#include "Field.h"
#include "GlobalConstants.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "TimeStepper.h"
#include "Tracers.h"
#include "VertCoord.h"
#include "VertMix.h"
#include "mpi.h"

#include <cmath>
#include <limits>
#include <string>
#include <vector>

using namespace OMEGA;

struct TestSetupPlane {

   Real Lx = 1;
   Real Ly = SqrtThree / 2;

   ErrorMeasures ExpectedDivErrors              = {0.00124886886594453264,
                                                   0.00124886886590977139};
   ErrorMeasures ExpectedPVErrors               = {0.00807347170900282914,
                                                   0.00794755105765788429};
   ErrorMeasures ExpectedGradErrors             = {0.00125026071878537952,
                                                   0.00134354611117262161};
   ErrorMeasures ExpectedLaplaceErrors          = {0.00113090174765822192,
                                                   0.00134324628763667899};
   ErrorMeasures ExpectedTrHAdvErrors           = {0.0029211089892916243,
                                                   0.0024583038518548855};
   ErrorMeasures ExpectedFCTErrors              = {0.00, 0.00};
   ErrorMeasures ExpectedFCTHProv               = {0.00, 0.00};
   ErrorMeasures ExpectedFCTHInv                = {0.00, 0.00};
   ErrorMeasures ExpectedFCTHNew                = {0.00, 0.00};
   ErrorMeasures ExpectedFCT_High               = {1.0408340855860843e-17,
                                                   1.2836244918860014e-17};
   ErrorMeasures ExpectedFCT_Low                = {0, 0};
   ErrorMeasures ExpectedFCTToNonFCT            = {9.524542988823467e-16,
                                                   5.336197203674528e-16};
   ErrorMeasures ExpectedTrDel2Errors           = {0.00334357193650093847,
                                                   0.00290978146207349032};
   ErrorMeasures ExpectedTrDel4Errors           = {0.00508833446725232875,
                                                   0.00523080740758275625};
   ErrorMeasures ExpectedSurfTrRestErrors       = {0, 0};
   ErrorMeasures ExpectedSfcStressForcingErrors = {0, 0};
   ErrorMeasures ExpectedBottomDragErrors       = {0.033848740052302935,
                                                   0.01000133508329411};
   ErrorMeasures ExpectedCoriolis2DErrors       = {0.01322503349256595,
                                                   0.011816356596905431};
   ErrorMeasures ExpectedCoriolis1DErrors       = {0.012883757868753419,
                                                   0.011805932389752003};

   KOKKOS_FUNCTION Real velEdgeX(Real X, Real Y) const { return 2 / Lx; }

   KOKKOS_FUNCTION Real velEdgeY(Real X, Real Y) const { return 3 / Ly; }

   KOKKOS_FUNCTION Real vectorX(Real X, Real Y) const {
      return std::sin(TwoPi * X / Lx) * std::cos(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real vectorY(Real X, Real Y) const {
      return std::cos(TwoPi * X / Lx) * std::sin(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real divergence(Real X, Real Y) const {
      return TwoPi * (1. / Lx + 1. / Ly) * std::cos(TwoPi * X / Lx) *
             std::cos(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real scalar(Real X, Real Y) const {
      return std::sin(TwoPi * X / Lx) * std::sin(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real gradX(Real X, Real Y) const {
      return TwoPi / Lx * std::cos(TwoPi * X / Lx) * std::sin(TwoPi * Y / Ly);
   }
   KOKKOS_FUNCTION Real gradY(Real X, Real Y) const {
      return TwoPi / Ly * std::sin(TwoPi * X / Lx) * std::cos(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real curl(Real X, Real Y) const {
      return TwoPi * (-1. / Lx + 1. / Ly) * std::sin(TwoPi * X / Lx) *
             std::sin(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real laplaceVecX(Real X, Real Y) const {
      return -TwoPi * TwoPi * (1. / Lx / Lx + 1. / Ly / Ly) *
             std::sin(TwoPi * X / Lx) * std::cos(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real laplaceVecY(Real X, Real Y) const {
      return -TwoPi * TwoPi * (1. / Lx / Lx + 1. / Ly / Ly) *
             std::cos(TwoPi * X / Lx) * std::sin(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real pseudoThick(Real X, Real Y) const {
      return 2. + std::sin(TwoPi * X / Lx) * std::cos(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real planetaryVort(Real X, Real Y) const {
      return std::cos(TwoPi * X / Lx) * std::cos(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real normRelVort(Real X, Real Y) const {
      return curl(X, Y) / pseudoThick(X, Y);
   }

   KOKKOS_FUNCTION Real normPlanetVort(Real X, Real Y) const {
      return planetaryVort(X, Y) / pseudoThick(X, Y);
   }

   KOKKOS_FUNCTION Real tracerFluxDiv(Real X, Real Y) const {
      return (TwoPi / (Lx * Ly)) *
             (std::cos(TwoPi * X / Lx) *
              (2 * (Lx + Ly) * std::cos(TwoPi * Y / Ly) +
               (Lx + 2 * Ly) * std::sin(TwoPi * X / Lx) *
                   std::pow(std::cos(TwoPi * Y / Ly), 2) -
               Lx * std::sin(TwoPi * X / Lx) *
                   std::pow(std::sin(TwoPi * Y / Ly), 2)));
   }

   KOKKOS_FUNCTION Real scalarA(Real X, Real Y) const {
      return std::cos(TwoPi * X / Lx) * std::sin(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real scalarB(Real X, Real Y) const {
      return 2. + std::cos(TwoPi * X / Lx) * std::cos(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real tracerDiff(Real X, Real Y, Real EddyDiff2) const {
      return -EddyDiff2 * TwoPi * TwoPi * std::sin(TwoPi * Y / Ly) *
             (2 * (1 / Lx / Lx + 1 / Ly / Ly) * std::cos(TwoPi * X / Lx) +
              (1 / Ly / Ly +
               (1 / Lx / Lx + 1 / Ly / Ly) * std::cos(2 * TwoPi * X / Lx)) *
                  std::cos(TwoPi * Y / Ly));
   }

   KOKKOS_FUNCTION Real scalarC(Real X, Real Y) const {
      return std::pow(std::cos(TwoPi * X / Lx), 2) -
             std::pow(std::sin(TwoPi * Y / Ly), 2);
   }

   KOKKOS_FUNCTION Real tracerHyperDiff(Real X, Real Y, Real EddyDiff4) const {
      return -EddyDiff4 * 2 * TwoPi * TwoPi *
             (std::cos(2 * TwoPi * X / Lx) / Lx / Lx +
              std::cos(2 * TwoPi * Y / Ly) / Ly / Ly);
   }

   KOKKOS_FUNCTION Real sfcStressForcingX(Real X, Real Y) const {
      const Real StressU = vectorX(X, Y);
      const Real Thick   = scalarB(X, Y);
      return StressU / (Thick * RhoSw);
   }

   KOKKOS_FUNCTION Real sfcStressForcingY(Real X, Real Y) const {
      const Real StressV = vectorY(X, Y);
      const Real Thick   = scalarB(X, Y);
      return StressV / (Thick * RhoSw);
   }

   KOKKOS_FUNCTION Real bottomDragX(Real X, Real Y, Real Coeff) const {
      const Real UVel = vectorX(X, Y);
      return -Coeff * std::abs(scalarA(X, Y)) / scalarB(X, Y) * UVel;
   }

   KOKKOS_FUNCTION Real bottomDragY(Real X, Real Y, Real Coeff) const {
      const Real VVel = vectorY(X, Y);
      return -Coeff * std::abs(scalarA(X, Y)) / scalarB(X, Y) * VVel;
   }

}; // end TestSetupPlane

struct TestSetupSphere {
   // radius of spherical mesh
   Real Radius = REarth;

   ErrorMeasures ExpectedDivErrors        = {0.013659556526126423,
                                             0.0036698023569596279};
   ErrorMeasures ExpectedPVErrors         = {0.021937128290005589,
                                             0.012253081309866708};
   ErrorMeasures ExpectedGradErrors       = {0.0018790366180557886,
                                             0.0014984647360648147};
   ErrorMeasures ExpectedLaplaceErrors    = {0.2819223346929331,
                                             0.27053032782204006};
   ErrorMeasures ExpectedTrHAdvErrors     = {0.013259084476801913,
                                             0.004075236217689366};
   ErrorMeasures ExpectedFCTErrors        = {0.00, 0.00};
   ErrorMeasures ExpectedFCTHProv         = {3.054304092320326e-05,
                                             1.077950920692842e-06};
   ErrorMeasures ExpectedFCTHInv          = {3.054397382906693e-05,
                                             1.0779509220589763e-06};
   ErrorMeasures ExpectedFCTHNew          = {3.0543973828844884e-05,
                                             1.0779233406323438e-06};
   ErrorMeasures ExpectedFCT_High         = {0.00146484375, 16478.526025524854};
   ErrorMeasures ExpectedFCT_Low          = {0, 0};
   ErrorMeasures ExpectedFCTToNonFCT      = {1.5120472343676013e-15,
                                             6.215666716489436e-16};
   ErrorMeasures ExpectedTrDel2Errors     = {0.015620987792697782,
                                             0.0032365878896264888};
   ErrorMeasures ExpectedTrDel4Errors     = {0.00081985237645818541,
                                             0.00064699727012728107};
   ErrorMeasures ExpectedSurfTrRestErrors = {0, 0};
   ErrorMeasures ExpectedSfcStressForcingErrors = {0, 0};
   ErrorMeasures ExpectedBottomDragErrors       = {0.0015343843060931499,
                                                   0.0014861805971789941};
   ErrorMeasures ExpectedCoriolis2DErrors       = {0.017808748900211516,
                                                   0.0095867629560083354};
   ErrorMeasures ExpectedCoriolis1DErrors       = {0.017532349880532885,
                                                   0.011526365807152578};

   KOKKOS_FUNCTION Real velEdgeX(Real X, Real Y) const { return 2; }

   KOKKOS_FUNCTION Real velEdgeY(Real X, Real Y) const { return 3; }

   KOKKOS_FUNCTION Real vectorX(Real Lon, Real Lat) const {
      return -Radius * std::pow(std::sin(Lon), 2) * std::pow(std::cos(Lat), 3);
   }

   KOKKOS_FUNCTION Real vectorY(Real Lon, Real Lat) const {
      return -4 * Radius * std::sin(Lon) * std::cos(Lon) *
             std::pow(std::cos(Lat), 3) * std::sin(Lat);
   }

   KOKKOS_FUNCTION Real divergence(Real Lon, Real Lat) const {
      return std::sin(Lon) * std::cos(Lon) * std::pow(std::cos(Lat), 2) *
             (20 * std::pow(std::sin(Lat), 2) - 6);
   }

   KOKKOS_FUNCTION Real scalar(Real Lon, Real Lat) const {
      return Radius * std::cos(Lon) * std::pow(std::cos(Lat), 4);
   }

   KOKKOS_FUNCTION Real gradX(Real Lon, Real Lat) const {
      return -std::sin(Lon) * std::pow(std::cos(Lat), 3);
   }

   KOKKOS_FUNCTION Real gradY(Real Lon, Real Lat) const {
      return -4 * std::cos(Lon) * std::pow(std::cos(Lat), 3) * std::sin(Lat);
   }

   KOKKOS_FUNCTION Real curl(Real Lon, Real Lat) const {
      return -4 * std::pow(std::cos(Lon), 2) * std::pow(std::cos(Lat), 2) *
             std::sin(Lat);
   }

   KOKKOS_FUNCTION Real laplaceVecX(Real Lon, Real Lat) const {
      return std::cos(Lat) *
             (std::pow(std::sin(Lat), 2) *
                  (17 - 37 * std::pow(std::sin(Lon), 2)) +
              11 * std::pow(std::sin(Lon), 2) - 5) /
             Radius;
   }

   KOKKOS_FUNCTION Real laplaceVecY(Real Lon, Real Lat) const {
      return std::sin(Lon) * std::cos(Lon) * std::sin(Lat) * std::cos(Lat) *
             (96 * std::pow(std::cos(Lat), 2) - 22) / Radius;
   }

   KOKKOS_FUNCTION Real pseudoThick(Real Lon, Real Lat) const {
      return (2 + std::cos(Lon) * std::pow(std::cos(Lat), 4));
   }

   KOKKOS_FUNCTION Real planetaryVort(Real Lon, Real Lat) const {
      return std::sin(Lat);
   }

   KOKKOS_FUNCTION Real normRelVort(Real Lon, Real Lat) const {
      return curl(Lon, Lat) / pseudoThick(Lon, Lat);
   }

   KOKKOS_FUNCTION Real normPlanetVort(Real Lon, Real Lat) const {
      return planetaryVort(Lon, Lat) / pseudoThick(Lon, Lat);
   }

   KOKKOS_FUNCTION Real tracerFluxDiv(Real Lon, Real Lat) const {
      return std::sin(Lon) * std::pow(std::cos(Lat), 2) *
             (std::cos(Lon) * (8 - 20 * std::cos(2 * Lat)) -
              6 * std::pow(std::cos(Lon), 2) * std::pow(std::cos(Lat), 4) *
                  (-2 + 3 * std::cos(2 * Lat)) +
              std::pow(std::cos(Lat), 4) * std::pow(std::sin(Lon), 2));
   }

   KOKKOS_FUNCTION Real scalarA(Real Lon, Real Lat) const {
      return Radius * std::pow(std::sin(Lon), 2) * std::pow(std::cos(Lat), 2);
   }

   // scalarB is differentiated by the tracer diffusion test, so it must be
   // smooth on the sphere: cos(Lon)*cos(Lat)*sin(Lat) = x*z/R^2 vanishes at
   // the poles, where Lon is undefined.
   KOKKOS_FUNCTION Real scalarB(Real Lon, Real Lat) const {
      return 2. + std::cos(Lon) * std::cos(Lat) * std::sin(Lat);
   }

   // div(scalarB grad(scalarA)) on the sphere
   KOKKOS_FUNCTION Real tracerDiff(Real Lon, Real Lat, Real EddyDiff2) const {
      const Real CosLon = std::cos(Lon);
      const Real SinLon = std::sin(Lon);
      const Real CosLat = std::cos(Lat);
      const Real SinLat = std::sin(Lat);
      return 2 * EddyDiff2 *
             (2 * std::cos(2 * Lon) +
              CosLat * SinLat * CosLon *
                  (CosLon * CosLon - 2 * SinLon * SinLon) -
              SinLon * SinLon *
                  (2 * (CosLat * CosLat - 2 * SinLat * SinLat) +
                   CosLon * CosLat * SinLat *
                       (2 * CosLat * CosLat - 3 * SinLat * SinLat))) /
             Radius;
   }

   KOKKOS_FUNCTION Real scalarC(Real Lon, Real Lat) const {
      return -(Radius / 2) * std::sqrt(3. / 2. / Pi) * std::cos(Lat) *
             std::cos(Lon);
   }

   KOKKOS_FUNCTION Real tracerHyperDiff(Real Lon, Real Lat,
                                        Real EddyDiff4) const {
      return EddyDiff4 * std::sqrt(3. / 2. / Pi) * std::cos(Lat) *
             std::cos(Lon) / Radius;
   }

   KOKKOS_FUNCTION Real sfcStressForcingX(Real Lon, Real Lat) const {
      const Real StressU = vectorX(Lon, Lat);
      const Real Thick   = scalarB(Lon, Lat);
      return StressU / (Thick * RhoSw);
   }

   KOKKOS_FUNCTION Real sfcStressForcingY(Real Lon, Real Lat) const {
      const Real StressV = vectorY(Lon, Lat);
      const Real Thick   = scalarB(Lon, Lat);
      return StressV / (Thick * RhoSw);
   }

   KOKKOS_FUNCTION Real bottomDragX(Real Lon, Real Lat, Real Coeff) const {
      const Real UVel = vectorX(Lon, Lat);
      return -Coeff * std::abs(scalarA(Lon, Lat)) / scalarB(Lon, Lat) * UVel;
   }

   KOKKOS_FUNCTION Real bottomDragY(Real Lon, Real Lat, Real Coeff) const {
      const Real VVel = vectorY(Lon, Lat);
      return -Coeff * std::abs(scalarA(Lon, Lat)) / scalarB(Lon, Lat) * VVel;
   }

}; // end TestSetupSphere

#ifdef TENDENCYTERMS_TEST_PLANE
constexpr Geometry Geom          = Geometry::Planar;
constexpr char DefaultMeshFile[] = "OmegaPlanarMesh.nc";
using TestSetup                  = TestSetupPlane;
#else
constexpr Geometry Geom          = Geometry::Spherical;
constexpr char DefaultMeshFile[] = "OmegaSphereMesh.nc";
using TestSetup                  = TestSetupSphere;
#endif

int setupBottomDragTestFields(const int NVertLayers, const Real Coeff,
                              const Array2DReal &ExactBottomDrag,
                              const Array2DReal &NormalVelEdge,
                              const Array2DReal &KECell,
                              const Array2DReal &PseudoThickEdge) {
   int Err = 0;
   TestSetup Setup;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   // Note: this computes bottom drag at every layer.
   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.bottomDragX(X, Y, Coeff);
          VecField[1] = Setup.bottomDragY(X, Y, Coeff);
       },
       ExactBottomDrag, EdgeComponent::Normal, Geom, Mesh, ExchangeHalos::No);

   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   // Reset bottom drag to zero above the lowest layer.
   parallelFor(
       "ResetExactBottomDrag", {Mesh->NEdgesOwned, NVertLayers},
       KOKKOS_LAMBDA(int IEdge, int K) {
          if (K != MaxLayerEdgeTop(IEdge)) {
             ExactBottomDrag(IEdge, K) = 0.0_Real;
          }
       });

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.vectorX(X, Y);
          VecField[1] = Setup.vectorY(X, Y);
       },
       NormalVelEdge, EdgeComponent::Normal, Geom, Mesh);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) {
          return Setup.scalarA(X, Y) * Setup.scalarA(X, Y) / 2;
       },
       KECell, Geom, Mesh, OnCell);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.scalarB(X, Y); },
       PseudoThickEdge, Geom, Mesh, OnEdge);

   return Err;
}

int testThickFluxDiv(int NVertLayers, Real RTol) {

   int Err = 0;
   TestSetup Setup;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   // Compute exact result
   Array2DReal ExactThickFluxDiv("ExactThickFluxDiv", Mesh->NCellsOwned,
                                 NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return -Setup.divergence(X, Y); },
       ExactThickFluxDiv, Geom, Mesh, OnCell, ExchangeHalos::No);

   // Set input array
   Array2DReal ThickFluxEdge("ThickFluxEdge", Mesh->NEdgesSize, NVertLayers);

   // TODO(mwarusz) temporary fix for this test
   Array2DReal OnesEdge("OnesEdge", Mesh->NEdgesSize, NVertLayers);
   deepCopy(OnesEdge, 1);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.vectorX(X, Y);
          VecField[1] = Setup.vectorY(X, Y);
       },
       ThickFluxEdge, EdgeComponent::Normal, Geom, Mesh);

   // Compute numerical result
   Array2DReal NumThickFluxDiv("NumThickFluxDiv", Mesh->NCellsOwned,
                               NVertLayers);
   PseudoThicknessFluxDivOnCell ThickFluxDivOnC(Mesh, VCoord);
   parallelForOuter(
       LaunchConfig({Mesh->NCellsOwned}, TeamScratch<Real>(NVertLayers)),
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          ThickFluxDivOnC(Team, NumThickFluxDiv, ICell, OnesEdge,
                          ThickFluxEdge);
       });

   // Compute errors
   ErrorMeasures TFDivErrors;
   Err += computeErrors(TFDivErrors, NumThickFluxDiv, ExactThickFluxDiv, Mesh,
                        OnCell);

   // Check error values
   Err += checkErrors("TendencyTermsTest", "ThickFluxDiv", TFDivErrors,
                      Setup.ExpectedDivErrors, RTol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: ThickFluxDiv PASS");
   }

   return Err;
} // end testThickFluxDiv

int testPotVortHAdv(int NVertLayers, Real RTol) {

   int Err = 0;
   TestSetup Setup;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   // Compute exact result
   Array2DReal ExactPotVortHAdv("ExactPotVortHAdv", Mesh->NEdgesOwned,
                                NVertLayers);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = (Setup.normRelVort(X, Y) + Setup.normPlanetVort(X, Y)) *
                        Setup.pseudoThick(X, Y) * Setup.vectorX(X, Y);
          VecField[1] = (Setup.normRelVort(X, Y) + Setup.normPlanetVort(X, Y)) *
                        Setup.pseudoThick(X, Y) * Setup.vectorY(X, Y);
       },
       ExactPotVortHAdv, EdgeComponent::Tangential, Geom, Mesh,
       ExchangeHalos::No);

   // Set input arrays
   Array2DReal NormRelVortEdge("NormRelVortEdge", Mesh->NEdgesSize,
                               NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.normRelVort(X, Y); },
       NormRelVortEdge, Geom, Mesh, OnEdge);

   Array2DReal NormPlanetVortEdge("NormPlanetVortEdge", Mesh->NEdgesSize,
                                  NVertLayers);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.normPlanetVort(X, Y); },
       NormPlanetVortEdge, Geom, Mesh, OnEdge);

   Array2DReal PseudoThickEdge("PseudoThickEdge", Mesh->NEdgesSize,
                               NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.pseudoThick(X, Y); },
       PseudoThickEdge, Geom, Mesh, OnEdge);

   Array2DReal NormVelEdge("NormVelEdge", Mesh->NEdgesSize, NVertLayers);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.vectorX(X, Y);
          VecField[1] = Setup.vectorY(X, Y);
       },
       NormVelEdge, EdgeComponent::Normal, Geom, Mesh);

   // Compute numerical result
   Array2DReal NumPotVortHAdv("NumPotVortHAdv", Mesh->NEdgesOwned, NVertLayers);

   PotentialVortHAdvOnEdge PotVortHAdvOnE(Mesh, VCoord);
   parallelForOuter(
       LaunchConfig({Mesh->NEdgesOwned}, TeamScratch<Real>(NVertLayers)),
       KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          PotVortHAdvOnE(Team, NumPotVortHAdv, IEdge, NormRelVortEdge,
                         NormPlanetVortEdge, PseudoThickEdge, NormVelEdge);
       });

   // Compute errors
   ErrorMeasures PotVortHAdvErrors;
   Err += computeErrors(PotVortHAdvErrors, NumPotVortHAdv, ExactPotVortHAdv,
                        Mesh, OnEdge);

   // Check error values
   Err += checkErrors("TendencyTermsTest", "PotVortHAdv", PotVortHAdvErrors,
                      Setup.ExpectedPVErrors, RTol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: PotVortHAdv PASS");
   }

   return Err;
} // end testPotVortHAdv

int testCoriolisAccelerationOnEdge(int NVertLayers, Real RTol) {

   int Err = 0;
   TestSetup Setup;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   // Set input arrays

   Array2DReal NormalVelEdge("NormalVelEdge", Mesh->NEdgesSize, NVertLayers);
   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.vectorX(X, Y);
          VecField[1] = Setup.vectorY(X, Y);
       },
       NormalVelEdge, EdgeComponent::Normal, Geom, Mesh);

   Array1DReal NormalBarotropicVelEdge("NormalBarotropicVelEdge",
                                       Mesh->NEdgesSize);
   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.vectorY(X, Y);
          VecField[1] = Setup.vectorX(X, Y);
       },
       NormalBarotropicVelEdge, EdgeComponent::Normal, Geom, Mesh);

   Array1DReal FEdge("FEdge", Mesh->NEdgesSize);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.planetaryVort(X, Y); },
       FEdge, Geom, Mesh, OnEdge);

   // Compute exact results

   Array2DReal ExactCoriolis2D("ExactCoriolis2D", Mesh->NEdgesOwned,
                               NVertLayers);
   Array1DReal ExactCoriolis1D("ExactCoriolis1D", Mesh->NEdgesOwned);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          const Real PVort = Setup.planetaryVort(X, Y);
          VecField[0]      = PVort * Setup.vectorX(X, Y);
          VecField[1]      = PVort * Setup.vectorY(X, Y);
       },
       ExactCoriolis2D, EdgeComponent::Tangential, Geom, Mesh,
       ExchangeHalos::No);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          const Real PVort = Setup.planetaryVort(X, Y);
          VecField[0]      = PVort * Setup.vectorY(X, Y);
          VecField[1]      = PVort * Setup.vectorX(X, Y);
       },
       ExactCoriolis1D, EdgeComponent::Tangential, Geom, Mesh,
       ExchangeHalos::No);

   // Compute numerical results

   Array2DReal NumCoriolis2D("NumCoriolis2D", Mesh->NEdgesOwned, NVertLayers);
   Array1DReal NumCoriolis1D("NumCoriolis1D", Mesh->NEdgesOwned);

   CoriolisAccelerationOnEdge CoriolisAccelOnE(Mesh, VCoord);

   parallelForOuter(
       LaunchConfig({Mesh->NEdgesOwned},
                    TeamScratch<Real>(VCoord->NVertLayers)),
       KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          CoriolisAccelOnE(Team, NumCoriolis2D, NumCoriolis2D, IEdge,
                           NormalVelEdge, FEdge);
       });

   parallelFor(
       {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge) {
          CoriolisAccelOnE(NumCoriolis1D, IEdge, NormalBarotropicVelEdge,
                           FEdge);
       });

   // Compute errors and check error values

   ErrorMeasures Coriolis2DErrors;
   Err += computeErrors(Coriolis2DErrors, NumCoriolis2D, ExactCoriolis2D, Mesh,
                        OnEdge);
   Err += checkErrors("TendencyTermsTest", "CoriolisAcceleration2D",
                      Coriolis2DErrors, Setup.ExpectedCoriolis2DErrors, RTol);

   ErrorMeasures Coriolis1DErrors;
   Err += computeErrors(Coriolis1DErrors, NumCoriolis1D, ExactCoriolis1D, Mesh,
                        OnEdge);
   Err += checkErrors("TendencyTermsTest", "CoriolisAcceleration1D",
                      Coriolis1DErrors, Setup.ExpectedCoriolis1DErrors, RTol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: CoriolisAccelerationOnEdge PASS");
   }

   return Err;
} // end testCoriolisAccelerationOnEdge

int testKEGrad(int NVertLayers, Real RTol) {

   int Err = 0;
   TestSetup Setup;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   // Compute exact result
   Array2DReal ExactKEGrad("ExactKEGrad", Mesh->NEdgesOwned, NVertLayers);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = -Setup.gradX(X, Y);
          VecField[1] = -Setup.gradY(X, Y);
       },
       ExactKEGrad, EdgeComponent::Normal, Geom, Mesh, ExchangeHalos::No);

   // Set input array
   Array2DReal KECell("KECell", Mesh->NCellsSize, NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.scalar(X, Y); }, KECell,
       Geom, Mesh, OnCell);

   // Compute numerical result
   Array2DReal NumKEGrad("NumKEGrad", Mesh->NEdgesOwned, NVertLayers);

   KEGradOnEdge KEGradOnE(Mesh, VCoord);
   parallelForOuter(
       {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          KEGradOnE(Team, NumKEGrad, IEdge, KECell);
       });

   // Compute errors
   ErrorMeasures KEGradErrors;
   Err += computeErrors(KEGradErrors, NumKEGrad, ExactKEGrad, Mesh, OnEdge);

   // Check error values
   Err += checkErrors("TendencyTermsTest", "KEGrad", KEGradErrors,
                      Setup.ExpectedGradErrors, RTol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: KEGrad PASS");
   }

   return Err;
} // end testKEGrad

int testSSHGrad(int NVertLayers, Real RTol) {

   int Err = 0;
   TestSetup Setup;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   // Compute exact result
   Array2DReal ExactSSHGrad("ExactSSHGrad", Mesh->NEdgesOwned, NVertLayers);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = -Gravity * Setup.gradX(X, Y);
          VecField[1] = -Gravity * Setup.gradY(X, Y);
       },
       ExactSSHGrad, EdgeComponent::Normal, Geom, Mesh, ExchangeHalos::No);

   Array1DReal SSHCell("SSHCell", Mesh->NCellsSize);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.scalar(X, Y); }, SSHCell,
       Geom, Mesh, OnCell);

   // Compute numerical result
   Array2DReal NumSSHGrad("NumSSHGrad", Mesh->NEdgesOwned, NVertLayers);

   SSHGradOnEdge SSHGradOnE(Mesh, VCoord);
   parallelForOuter(
       {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          SSHGradOnE(Team, NumSSHGrad, IEdge, SSHCell);
       });

   // Compute errors
   ErrorMeasures SSHGradErrors;
   Err += computeErrors(SSHGradErrors, NumSSHGrad, ExactSSHGrad, Mesh, OnEdge);

   // Check error values
   Err += checkErrors("TendencyTermsTest", "SSHGrad", SSHGradErrors,
                      Setup.ExpectedGradErrors, RTol);

   return Err;
} // end testSSHGrad

int testVelDiff(int NVertLayers, Real RTol) {

   int Err = 0;
   Error Err1;
   TestSetup Setup;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   Config *OmegaConfig = Config::getOmegaConfig();
   Config TendConfig("Tendencies");
   Err1 += OmegaConfig->get(TendConfig);
   CHECK_ERROR_ABORT(Err1, "Tendencies: Tendencies group not found in Config");

   VelocityDiffusionOnEdge VelDiffOnE(Mesh, VCoord);
   Err1 += TendConfig.get("ViscDel2", VelDiffOnE.ViscDel2);
   CHECK_ERROR_ABORT(Err1, "Tendencies: ViscDel2 not found in TendConfig");

   const Real ViscDel2 = VelDiffOnE.ViscDel2;

   // Compute exact result
   Array2DReal ExactVelDiff("ExactVelDiff", Mesh->NEdgesOwned, NVertLayers);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = ViscDel2 * Setup.laplaceVecX(X, Y);
          VecField[1] = ViscDel2 * Setup.laplaceVecY(X, Y);
       },
       ExactVelDiff, EdgeComponent::Normal, Geom, Mesh, ExchangeHalos::No);

   // Set input arrays
   Array2DReal DivCell("DivCell", Mesh->NCellsSize, NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.divergence(X, Y); },
       DivCell, Geom, Mesh, OnCell);

   Array2DReal RVortVertex("RVortVertex", Mesh->NVerticesSize, NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.curl(X, Y); }, RVortVertex,
       Geom, Mesh, OnVertex);

   // Compute numerical result
   Array2DReal NumVelDiff("NumVelDiff", Mesh->NEdgesOwned, NVertLayers);

   parallelForOuter(
       {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          VelDiffOnE(Team, NumVelDiff, IEdge, DivCell, RVortVertex);
       });

   // Compute errors
   ErrorMeasures VelDiffErrors;
   Err += computeErrors(VelDiffErrors, NumVelDiff, ExactVelDiff, Mesh, OnEdge);

   // Check error values
   Err += checkErrors("TendencyTermsTest", "VelDiff", VelDiffErrors,
                      Setup.ExpectedLaplaceErrors, RTol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: VelDiff PASS");
   }

   return Err;
} // end testVelDiff

int testVelHyperDiff(int NVertLayers, Real RTol) {

   int Err = 0;
   Error Err1;
   TestSetup Setup;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   Config *OmegaConfig = Config::getOmegaConfig();
   Config TendConfig("Tendencies");
   Err1 += OmegaConfig->get(TendConfig);
   CHECK_ERROR_ABORT(Err1, "Tendencies: Tendencies group not found in Config");

   VelocityHyperDiffOnEdge VelHyperDiffOnE(Mesh, VCoord);
   Err1 += TendConfig.get("ViscDel4", VelHyperDiffOnE.ViscDel4);
   CHECK_ERROR_ABORT(Err1, "Tendencies: ViscDel4 not found in TendConfig");

   Err1 += TendConfig.get("DivFactor", VelHyperDiffOnE.DivFactor);
   CHECK_ERROR_ABORT(Err1, "Tendencies: DivFactor not found in TendConfig");

   const Real ViscDel4 = VelHyperDiffOnE.ViscDel4;

   // Compute exact result
   Array2DReal ExactVelHyperDiff("ExactVelHyperDiff", Mesh->NEdgesOwned,
                                 NVertLayers);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = -ViscDel4 * Setup.laplaceVecX(X, Y);
          VecField[1] = -ViscDel4 * Setup.laplaceVecY(X, Y);
       },
       ExactVelHyperDiff, EdgeComponent::Normal, Geom, Mesh, ExchangeHalos::No);

   // Set input arrays
   Array2DReal DivCell("DivCell", Mesh->NCellsSize, NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.divergence(X, Y); },
       DivCell, Geom, Mesh, OnCell);

   Array2DReal RVortVertex("RVortVertex", Mesh->NVerticesSize, NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.curl(X, Y); }, RVortVertex,
       Geom, Mesh, OnVertex);

   // Compute numerical result
   Array2DReal NumVelHyperDiff("NumVelHyperDiff", Mesh->NEdgesOwned,
                               NVertLayers);

   parallelForOuter(
       {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          VelHyperDiffOnE(Team, NumVelHyperDiff, IEdge, DivCell, RVortVertex);
       });

   // Compute errors
   ErrorMeasures VelHyperDiffErrors;
   Err += computeErrors(VelHyperDiffErrors, NumVelHyperDiff, ExactVelHyperDiff,
                        Mesh, OnEdge);

   // Check error values
   Err += checkErrors("TendencyTermsTest", "VelHyperDiff", VelHyperDiffErrors,
                      Setup.ExpectedLaplaceErrors, RTol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: VelHyperDiff PASS");
   }

   return Err;
} // end testVelHyperDiff

int testSfcStressForcing(int NVertLayers) {

   int Err = 0;
   TestSetup Setup;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   // Compute exact result
   Array2DReal ExactSfcStressForcing("ExactSfcStressForcing", Mesh->NEdgesOwned,
                                     NVertLayers);

   // Note: this computes surface stress forcing at every layer
   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.sfcStressForcingX(X, Y);
          VecField[1] = Setup.sfcStressForcingY(X, Y);
       },
       ExactSfcStressForcing, EdgeComponent::Normal, Geom, Mesh,
       ExchangeHalos::No);

   // Reset surface stress forcing to zero below the surface
   deepCopy(Kokkos::subview(ExactSfcStressForcing, Kokkos::ALL,
                            Kokkos::make_pair(1, NVertLayers)),
            0);

   // Set input arrays
   Array1DReal NormalStressEdge("NormalStressEdge", Mesh->NEdgesSize);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.vectorX(X, Y);
          VecField[1] = Setup.vectorY(X, Y);
       },
       NormalStressEdge, EdgeComponent::Normal, Geom, Mesh);

   Array2DReal PseudoThickEdge("PseudoThickEdge", Mesh->NEdgesSize,
                               NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.scalarB(X, Y); },
       PseudoThickEdge, Geom, Mesh, OnEdge);

   // Compute numerical result
   Array2DReal NumSfcStressForcing("NumSfcStressForcing", Mesh->NEdgesOwned,
                                   NVertLayers);

   SfcStressForcingOnEdge SfcStressForcingOnE(Mesh, VCoord);

   parallelForOuter(
       {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          SfcStressForcingOnE(Team, NumSfcStressForcing, IEdge,
                              NormalStressEdge, PseudoThickEdge);
       });

   // Compute errors
   ErrorMeasures SfcStressForcingErrors;
   Err += computeErrors(SfcStressForcingErrors, NumSfcStressForcing,
                        ExactSfcStressForcing, Mesh, OnEdge);

   // Check error values
   const Real RTol = 0;
   const Real ATol = 100 * std::numeric_limits<Real>::epsilon();
   Err += checkErrors("TendencyTermsTest", "SfcStressForcing",
                      SfcStressForcingErrors,
                      Setup.ExpectedSfcStressForcingErrors, RTol, ATol);

   return Err;
} // end testSfcStressForcing

int testExplicitBottomDrag(int NVertLayers, Real RTol) {

   int Err = 0;
   TestSetup Setup;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   const Real Coeff = 1.123456789;

   // Compute exact result
   Array2DReal ExactBottomDrag("ExactBottomDrag", Mesh->NEdgesOwned,
                               NVertLayers);

   // Set input arrays
   Array2DReal NormalVelEdge("NormalVelEdge", Mesh->NEdgesSize, NVertLayers);

   Array2DReal KECell("KECell", Mesh->NCellsSize, NVertLayers);

   Array2DReal PseudoThickEdge("PseudoThickEdge", Mesh->NEdgesSize,
                               NVertLayers);

   Err += setupBottomDragTestFields(NVertLayers, Coeff, ExactBottomDrag,
                                    NormalVelEdge, KECell, PseudoThickEdge);

   // Compute numerical result
   Array2DReal NumBottomDrag("NumBottomDrag", Mesh->NEdgesOwned, NVertLayers);

   BottomDragOnEdge BottomDragOnE(Mesh, VCoord);
   BottomDragOnE.Coeff = Coeff;

   parallelFor(
       {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge) {
          BottomDragOnE(NumBottomDrag, IEdge, NormalVelEdge, KECell,
                        PseudoThickEdge);
       });

   // Compute errors
   ErrorMeasures BottomDragErrors;
   Err += computeErrors(BottomDragErrors, NumBottomDrag, ExactBottomDrag, Mesh,
                        OnEdge);

   // Check error values
   Err += checkErrors("TendencyTermsTest", "ExplicitBottomDrag",
                      BottomDragErrors, Setup.ExpectedBottomDragErrors, RTol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: ExplicitBottomDrag PASS");
   }

   return Err;
} // end testExplicitBottomDrag

int testImplicitBottomDrag(int NVertLayers, Real RTol) {

   int Err = 0;
   TestSetup Setup;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   const Real Coeff = 1.123456789_Real;
   const Real DT    = 5.0_Real;

   Array2DReal SpecVol("SpecVol", Mesh->NCellsSize, NVertLayers);
   Array2DReal KineticEnergyCell("KineticEnergyCell", Mesh->NCellsSize,
                                 NVertLayers);
   Array2DReal PseudoThickCell("PseudoThickCell", Mesh->NCellsSize,
                               NVertLayers);
   Array2DReal PseudoThickEdge("PseudoThickEdge", Mesh->NEdgesSize,
                               NVertLayers);
   Array2DReal VertVisc("VertVisc", Mesh->NCellsSize, NVertLayers + 1);
   Array2DReal NormalVelEdge("NormalVelEdge", Mesh->NEdgesSize, NVertLayers);
   Array2DReal NumImplicitBottomDrag("NumImplicitBottomDrag", Mesh->NEdgesOwned,
                                     NVertLayers);
   Array2DReal ExactBottomDrag("ExactBottomDrag", Mesh->NEdgesOwned,
                               NVertLayers);

   VelVertMixSetupOnEdge ImplicitBottomDragOnE(Mesh, VCoord);
   ImplicitBottomDragOnE.ImplicitBottomDragEnabled = true;
   ImplicitBottomDragOnE.BottomDragCoeff           = Coeff;

   Err += setupBottomDragTestFields(NVertLayers, Coeff, ExactBottomDrag,
                                    NormalVelEdge, KineticEnergyCell,
                                    PseudoThickEdge);

   deepCopy(SpecVol, 1.0_Real / ImplicitBottomDragOnE.LocRhoSw);
   deepCopy(VertVisc, 0.0_Real);
   deepCopy(NumImplicitBottomDrag, 0.0_Real);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.scalarB(X, Y); },
       PseudoThickCell, Geom, Mesh, OnCell);

   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);
   OMEGA_SCOPE(CellsOnEdge, Mesh->CellsOnEdge);

   // Compute numerical result
   parallelFor(
       "NumImplicitBottomDrag", {Mesh->NEdgesOwned, NVertLayers},
       KOKKOS_LAMBDA(int IEdge, int K) {
          const int KMin = MinLayerEdgeBot(IEdge);
          const int KMax = MaxLayerEdgeTop(IEdge);
          if (K < KMin || K > KMax)
             return;

          Real G, H, X;
          ImplicitBottomDragOnE(IEdge, K, KMin, KMax, DT, SpecVol,
                                KineticEnergyCell, PseudoThickCell, VertVisc,
                                NormalVelEdge, G, H, X);

          const int JCell0 = CellsOnEdge(IEdge, 0);
          const int JCell1 = CellsOnEdge(IEdge, 1);
          const Real PseudoThickCellEdge =
              0.5_Real *
              (PseudoThickCell(JCell0, K) + PseudoThickCell(JCell1, K));
          const Real DragDiagContribution = H - PseudoThickCellEdge;

          NumImplicitBottomDrag(IEdge, K) = -DragDiagContribution *
                                            NormalVelEdge(IEdge, K) /
                                            (DT * PseudoThickEdge(IEdge, K));
       });

   // Compute errors
   ErrorMeasures ImplicitBottomDragErrors;
   Err += computeErrors(ImplicitBottomDragErrors, NumImplicitBottomDrag,
                        ExactBottomDrag, Mesh, OnEdge);

   // Check error values
   Err += checkErrors("TendencyTermsTest", "ImplicitBottomDrag",
                      ImplicitBottomDragErrors, Setup.ExpectedBottomDragErrors,
                      RTol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: ImplicitBottomDrag PASS");
   }

   return Err;
} // end testImplicitBottomDrag

/// Bottom drag must leave edges with no active layer on both sides alone.
///
/// VertCoord sets MaxLayerEdgeTop to -1 on land edges and on the outermost
/// edges of the halo, and Tendencies applies bottom drag with a flat loop over
/// all edges rather than over the chunked vertical range that the other
/// tendency terms use. Without an explicit guard the functor indexes the
/// bottom layer at -1, reading and writing outside the intended row: for a
/// LayoutRight array Tend(IEdge, -1) aliases Tend(IEdge - 1, NVertLayers - 1),
/// so an inactive edge silently corrupts its neighbor's deepest layer.
///
/// Mark a subset of edges inactive, then compare a pass over all edges against
/// a reference pass over only the active ones. The two must agree exactly.
int testBottomDragInactiveEdges(int NVertLayers) {

   int Err = 0;
   TestSetup Setup;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   const Real Coeff = 1.123456789;

   // Set input arrays
   Array2DReal NormalVelEdge("NormalVelEdge", Mesh->NEdgesSize, NVertLayers);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.vectorX(X, Y);
          VecField[1] = Setup.vectorY(X, Y);
       },
       NormalVelEdge, EdgeComponent::Normal, Geom, Mesh);

   Array2DReal KECell("KECell", Mesh->NCellsSize, NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) {
          return Setup.scalarA(X, Y) * Setup.scalarA(X, Y) / 2;
       },
       KECell, Geom, Mesh, OnCell);

   Array2DReal PseudoThickEdge("PseudoThickEdge", Mesh->NEdgesSize,
                               NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.scalarB(X, Y); },
       PseudoThickEdge, Geom, Mesh, OnEdge);

   // Save the layer indices so they can be restored for the later tests
   Array1DI4 SavedMaxLayerEdgeTop("SavedMaxLayerEdgeTop", Mesh->NEdgesSize);
   deepCopy(SavedMaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   // Mark every third edge inactive, the way VertCoord marks land edges and
   // the outermost edges of the halo
   OMEGA_SCOPE(LocMaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);
   parallelFor(
       {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge) {
          if (IEdge % 3 == 0)
             LocMaxLayerEdgeTop(IEdge) = -1;
       });

   // The functor keeps a shallow copy of MaxLayerEdgeTop, so construct it
   // after the array has been modified to make the intent clear
   BottomDragOnEdge BottomDragOnE(Mesh, VCoord);
   BottomDragOnE.Coeff = Coeff;

   // Reference: only the active edges contribute
   Array2DReal RefBottomDrag("RefBottomDrag", Mesh->NEdgesSize, NVertLayers);
   deepCopy(RefBottomDrag, 0);
   parallelFor(
       {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge) {
          if (LocMaxLayerEdgeTop(IEdge) >= 0)
             BottomDragOnE(RefBottomDrag, IEdge, NormalVelEdge, KECell,
                           PseudoThickEdge);
       });

   // What Tendencies actually does: loop over every edge unconditionally
   Array2DReal NumBottomDrag("NumBottomDrag", Mesh->NEdgesSize, NVertLayers);
   deepCopy(NumBottomDrag, 0);
   parallelFor(
       {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge) {
          BottomDragOnE(NumBottomDrag, IEdge, NormalVelEdge, KECell,
                        PseudoThickEdge);
       });

   // Restore the layer indices before checking, so that returning early on
   // failure cannot leave the vertical coordinate modified for later tests
   deepCopy(VCoord->MaxLayerEdgeTop, SavedMaxLayerEdgeTop);

   I4 NumBad = 0;
   parallelReduce(
       {Mesh->NEdgesAll, NVertLayers},
       KOKKOS_LAMBDA(int IEdge, int K, I4 &Accum) {
          const Real Num = NumBottomDrag(IEdge, K);
          const Real Ref = RefBottomDrag(IEdge, K);
          if (Num != Ref or Kokkos::isnan(Num) or Kokkos::isinf(Num))
             Accum += 1;
       },
       NumBad);

   if (NumBad > 0) {
      LOG_ERROR("TendencyTermsTest: BottomDragInactiveEdges FAIL, {} entries "
                "differ from the active-edge-only result",
                NumBad);
      Err += 1;
   } else {
      LOG_INFO("TendencyTermsTest: BottomDragInactiveEdges PASS");
   }

   return Err;
} // end testBottomDragInactiveEdges

int testTracerHorzAdvOnCell(int NVertLayers, int NTracers, Real RTol) {

   I4 Err = 0;
   TestSetup Setup;
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   // Compute exact result
   Array3DReal ExactTrFluxDiv("ExactTrFluxDiv", NTracers, Mesh->NCellsOwned,
                              NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.tracerFluxDiv(X, Y); },
       ExactTrFluxDiv, Geom, Mesh, OnCell, ExchangeHalos::No);

   // Set input arrays
   Array2DReal NormalVelocity("NormalVelocity", Mesh->NEdgesSize, NVertLayers);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.vectorX(X, Y);
          VecField[1] = Setup.vectorY(X, Y);
       },
       NormalVelocity, EdgeComponent::Normal, Geom, Mesh);

   Array3DReal TrCell("TrCell", NTracers, Mesh->NCellsSize, NVertLayers);
   Array2DReal ThickEdge("ThickEdh", Mesh->NEdgesSize, NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return -Setup.pseudoThick(X, Y); },
       TrCell, Geom, Mesh, OnCell);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return 1; }, ThickEdge, Geom, Mesh,
       OnEdge);

   // Compute numerical result
   Array3DReal NumTrFluxDiv("NumTrFluxDiv", NTracers, Mesh->NCellsOwned,
                            NVertLayers);

   VertAdv::init();
   const auto VAdv = VertAdv::getDefault();

   TracerHorzAdvOnCell TrHorzAdvOnC(Mesh, VCoord, VAdv);
   TrHorzAdvOnC.ForceLowOrder = true;
   TrHorzAdvOnC.init();

   parallelForOuter(
       LaunchConfig({NTracers, Mesh->NEdgesAll},
                    TeamScratch<Real>(NVertLayers)),
       KOKKOS_LAMBDA(int L, int IEdge, const TeamMember &Team) {
          TrHorzAdvOnC(Team, L, IEdge, TrCell, ThickEdge, NormalVelocity);
       });

   parallelForOuter(
       LaunchConfig({NTracers, Mesh->NCellsOwned},
                    TeamScratch<Real>(NVertLayers)),
       KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
          TrHorzAdvOnC(Team, NumTrFluxDiv, L, ICell);
       });

   ErrorMeasures TrHAdvErrors;
   Err +=
       computeErrors(TrHAdvErrors, NumTrFluxDiv, ExactTrFluxDiv, Mesh, OnCell);

   Err += checkErrors("TendencyTermsTest", "TracerHorzAdv", TrHAdvErrors,
                      Setup.ExpectedTrHAdvErrors, RTol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: TracerHorzAdv PASS");
   }

   VertAdv::clear();
   return Err;
} // end testTracerHorzAdvOnCell

class TracerHorzAdvOnCellTest : public TracerHorzAdvOnCell {
 public:
   TracerHorzAdvOnCellTest(const HorzMesh *Mesh, const VertCoord *VCoord,
                           const VertAdv *VAdv)
       : TracerHorzAdvOnCell(Mesh, VCoord, VAdv) {};
   Array2DReal GetHProvInv() const { return HProvInv; };
   Array2DReal GetHProv() const { return HProv; };
   Array2DReal GetHNewInv() const { return HNewInv; };
   Array2DReal GetTracerCur() const { return TracerCur; };
   Array2DReal GetTracerMin() const { return TracerMin; };
   Array2DReal GetTracerMax() const { return TracerMax; };
   Array3DReal GetHighOrderFlxHorz() const { return HighOrderFlxHorz; };
   Array2DReal GetLowOrderFlx() const { return LowOrderFlx; };
   Array2DReal GetHighOrderFlx() const { return HighOrderFlx; };
   Array2DReal GetFlxIn() const { return FlxIn; };
   Array2DReal GetFlxOut() const { return FlxOut; };
};

int testFCTTracerHorzAdvOnCell(int NVertLayers, int NTracers, Real RTol) {

   I4 Err = 0, TotErr = 0;
   TestSetup Setup;
   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   const int NCellsFCT = Mesh->NCellsSize;
   const int NEdgesFCT = Mesh->NEdgesSize;
   const int NVertsFCT = NVertLayers + 1;

   // Compute exact result
   Array3DReal ExactTrFluxDiv("ExactTrFluxDiv", NTracers, Mesh->NCellsOwned,
                              NVertsFCT);

   setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.tracerFluxDiv(X, Y); },
       ExactTrFluxDiv, Geom, Mesh, OnCell, ExchangeHalos::No);

   // Set input arrays
   Array2DReal NormalVelocity("NormalVelocity", NEdgesFCT, NVertsFCT);

   Array3DReal TrCell("TrCell", NTracers, NCellsFCT, NVertsFCT);
   Array2DReal Tend("Tend", NCellsFCT, NVertsFCT);
   Array3DReal TendNoFCT("TendNoFCT", NTracers, NCellsFCT, NVertsFCT);
   Array2DReal FluxPseudoThickEdge("FluxPseudoThickEdge", NEdgesFCT, NVertsFCT);
   Array2DReal LayerThickness("LayerThickness", NCellsFCT, NVertsFCT);
   Array2DReal NormVelEdge("NormVelEdge", NEdgesFCT, NVertsFCT);
   Array2DReal HProvInvExact("HProvInvExact", NCellsFCT, NVertsFCT);
   Array2DReal HProvExact("HProvExact", NCellsFCT, NVertsFCT);
   Array2DReal HNewInvExact("HProvExact", NCellsFCT, NVertsFCT);
   Array2DReal TracerSubView("TracerSubView", NCellsFCT, NVertsFCT);
   Array2DReal CellSubView("CellSubView", NCellsFCT, NVertsFCT);
   Array2DReal FluxSubView("FluxSubView", NEdgesFCT, NVertsFCT);
   Array2DReal ReferenceSolution("ReferenceSolution", NEdgesFCT, NVertsFCT);
   deepCopy(TrCell, 1._Real);
   deepCopy(Tend, 0._Real);
   deepCopy(TendNoFCT, 0._Real);
   deepCopy(Kokkos::subview(TrCell, Kokkos::ALL, Kokkos::ALL, NVertLayers),
            0._Real);
   deepCopy(FluxPseudoThickEdge, 1._Real);
   deepCopy(LayerThickness, 1._Real);
   deepCopy(HProvExact, 1._Real);
   deepCopy(Kokkos::subview(HProvExact, Kokkos::ALL, NVertLayers), 0._Real);
   deepCopy(HProvInvExact, 1._Real);
   deepCopy(Kokkos::subview(HProvInvExact, Kokkos::ALL, NVertLayers), 0._Real);
   deepCopy(HNewInvExact, 1._Real);
   deepCopy(Kokkos::subview(HNewInvExact, Kokkos::ALL, NVertLayers), 0._Real);
   deepCopy(NormVelEdge, 1._Real);
   deepCopy(ReferenceSolution, 0._Real);
   deepCopy(CellSubView, 0._Real);

   // Compute numerical result
   Array3DReal NumTrFluxDiv("NumTrFluxDiv", NTracers, Mesh->NCellsOwned,
                            NVertsFCT);

   VertAdv::init();
   const auto VAdv = VertAdv::getDefault();
   deepCopy(VAdv->VerticalPseudoVelocity, 10._Real);
   deepCopy(VAdv->TotalVerticalPseudoVelocity, 10._Real);
   deepCopy(VAdv->TotalVerticalTransportPseudoVelocity, 10._Real);

   TracerHorzAdvOnCellTest TrHorzAdvOnC(Mesh, VCoord, VAdv);
   TrHorzAdvOnC.ForceLowOrder = false;
   TrHorzAdvOnC.FCT           = true;
   TrHorzAdvOnC.init();

   const Real Dt = 1;
   ErrorMeasures FCTErrors;

   {
      setVectorEdge(
          KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
             VecField[0] = Setup.velEdgeX(X, Y);
             VecField[1] = Setup.velEdgeY(X, Y);
          },
          NormalVelocity, EdgeComponent::Normal, Geom, Mesh);

      parallelFor(
          {Mesh->NCellsOwned, NVertLayers},
          KOKKOS_LAMBDA(const int ICell, const int K) {
             TrHorzAdvOnC.FCTProvisionaLayerThicknesses(
                 ICell, K, Dt, FluxPseudoThickEdge, LayerThickness,
                 NormalVelocity);
          });

      const auto HProvInv = TrHorzAdvOnC.GetHProvInv();
      const auto HProv    = TrHorzAdvOnC.GetHProv();
      const auto HNewInv  = TrHorzAdvOnC.GetHNewInv();

      const Real ATol = sizeof(Real) == 4 ? 1e-4 : 1e-10;
      Err += computeErrors(FCTErrors, HProv, HProvExact, Mesh, OnCell);
      Err += checkErrors("TendencyTermsTest", "FCTHProv", FCTErrors,
                         Setup.ExpectedFCTHProv, RTol, ATol);

      if (Err == 0)
         LOG_INFO("TendencyTermsTest: FCTHProv PASS");
      TotErr += Err;
      Err = 0;

      Err += computeErrors(FCTErrors, HProvInv, HProvInvExact, Mesh, OnCell);
      Err += checkErrors("TendencyTermsTest", "FCTHProvInv", FCTErrors,
                         Setup.ExpectedFCTHInv, RTol, ATol);
      if (Err == 0)
         LOG_INFO("TendencyTermsTest: FCTHProvInv PASS");
      TotErr += Err;
      Err = 0;

      Err += computeErrors(FCTErrors, HNewInv, HNewInvExact, Mesh, OnCell);
      Err += checkErrors("TendencyTermsTest", "FCTHNewInv", FCTErrors,
                         Setup.ExpectedFCTHNew, RTol, ATol);
      if (Err == 0)
         LOG_INFO("TendencyTermsTest: FCTHNewInv PASS");
      TotErr += Err;
      Err = 0;
   }
   for (int L = 0; L < NTracers; ++L) {
      parallelFor(
          {Mesh->NCellsAll, NVertLayers},
          KOKKOS_LAMBDA(const int ICell, const int K) {
             TrCell(L, ICell, K) = L + ICell + K;
          });
      parallelFor(
          {Mesh->NCellsAll, NVertLayers},
          KOKKOS_LAMBDA(const int ICell, const int K) {
             TrHorzAdvOnC.FCTTracerCurFill(L, ICell, K, TrCell);
          });
      const Array2DReal TracerCur = TrHorzAdvOnC.GetTracerCur();
      const Array2DReal Tracers =
          Kokkos::subview(TrCell, L, Kokkos::ALL, Kokkos::ALL);
      deepCopy(TracerSubView, Tracers);
      Err += computeErrors(FCTErrors, TracerCur, TracerSubView, Mesh, OnCell);
      Err += checkErrors("TendencyTermsTest", "FCTTracerCurFill", FCTErrors,
                         Setup.ExpectedFCTErrors, RTol);
      if (Err == 0)
         LOG_INFO("TendencyTermsTest: FCTTracerCurFill_" + std::to_string(L) +
                  " PASS");
      TotErr += Err;
      Err = 0;
   }
   {
      setVectorEdge(
          KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
             VecField[0] = Setup.velEdgeX(X, Y);
             VecField[1] = Setup.velEdgeY(X, Y);
          },
          NormalVelocity, EdgeComponent::Normal, Geom, Mesh);
      deepCopy(TrCell, 1._Real);
      deepCopy(Kokkos::subview(TrCell, Kokkos::ALL, Kokkos::ALL, NVertLayers),
               0._Real);
      const int L = 0;
      parallelFor(
          {Mesh->NCellsAll, NVertLayers},
          KOKKOS_LAMBDA(const int ICell, const int K) {
             TrHorzAdvOnC.FCTTracerCurFill(L, ICell, K, TrCell);
          });
      parallelFor(
          {Mesh->NCellsOwned, NVertLayers},
          KOKKOS_LAMBDA(const int ICell, const int K) {
             TrHorzAdvOnC.FCTProvisionaLayerThicknesses(
                 ICell, K, Dt, FluxPseudoThickEdge, LayerThickness,
                 NormalVelocity);
          });
      parallelForOuter(
          {Mesh->NEdgesHaloH(1)},
          KOKKOS_LAMBDA(const int IEdge, const TeamMember &Team) {
             TrHorzAdvOnC.FCTHighAndLowOrderFlux(
                 Team, IEdge, FluxPseudoThickEdge, NormalVelocity);
          });

      Array2DReal LowOrderFlx = TrHorzAdvOnC.GetLowOrderFlx();
      const auto DvEdge       = Mesh->DvEdge;

      parallelFor(
          {Mesh->NEdgesHaloH(1), NVertLayers},
          KOKKOS_LAMBDA(const int IEdge, const int K) {
             FluxSubView(IEdge, K) = DvEdge(IEdge) * NormalVelocity(IEdge, K);
          });

      const Real ATol = 1.0e-10;
      Err += computeErrors(FCTErrors, LowOrderFlx, FluxSubView, Mesh, OnEdge);
      Err += checkErrors("TendencyTermsTest", "FCTHighAndLowOrderFlux_Low",
                         FCTErrors, Setup.ExpectedFCT_Low, RTol, ATol);
      if (Err == 0)
         LOG_INFO("TendencyTermsTest: FCTHighAndLowOrderFlux_Low PASS");
      TotErr += Err;
      Err = 0;
   }
   if (Geom == Geometry::Planar) {
      setVectorEdge(
          KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
             VecField[0] = Setup.vectorX(X, Y);
             VecField[1] = Setup.vectorY(X, Y);
          },
          NormalVelocity, EdgeComponent::Normal, Geom, Mesh);

      deepCopy(TrCell, 1._Real);
      deepCopy(Kokkos::subview(TrCell, Kokkos::ALL, Kokkos::ALL, NVertLayers),
               0._Real);
      const int L = 0;
      parallelFor(
          {Mesh->NCellsAll, NVertLayers},
          KOKKOS_LAMBDA(const int ICell, const int K) {
             TrHorzAdvOnC.FCTTracerCurFill(L, ICell, K, TrCell);
          });
      parallelFor(
          {Mesh->NCellsOwned, NVertLayers},
          KOKKOS_LAMBDA(const int ICell, const int K) {
             TrHorzAdvOnC.FCTProvisionaLayerThicknesses(
                 ICell, K, Dt, FluxPseudoThickEdge, LayerThickness,
                 NormalVelocity);
          });
      parallelForOuter(
          {Mesh->NEdgesHaloH(1)},
          KOKKOS_LAMBDA(const int IEdge, const TeamMember &Team) {
             TrHorzAdvOnC.FCTHighAndLowOrderFlux(
                 Team, IEdge, FluxPseudoThickEdge, NormalVelocity);
          });

      parallelFor(
          {Mesh->NEdgesHaloH(1), NVertLayers},
          KOKKOS_LAMBDA(const int IEdge, const int K) {
             FluxSubView(IEdge, K) = 0;
          });

      const Real ATol          = sizeof(Real) == 4 ? 1e-8 : 1e-10;
      Array2DReal HighOrderFlx = TrHorzAdvOnC.GetHighOrderFlx();
      Err += computeErrors(FCTErrors, HighOrderFlx, FluxSubView, Mesh, OnEdge);
      Err += checkErrors("TendencyTermsTest", "FCTHighAndLowOrderFlux_High",
                         FCTErrors, Setup.ExpectedFCT_High, RTol, ATol);
      if (Err == 0)
         LOG_INFO("TendencyTermsTest: FCTHighAndLowOrderFlux_High PASS");
      TotErr += Err;
      Err = 0;
   }
   if (Geom == Geometry::Planar) {
      setVectorEdge(
          KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
             VecField[0] = Setup.velEdgeX(X, Y);
             VecField[1] = Setup.velEdgeY(X, Y);
          },
          NormalVelocity, EdgeComponent::Normal, Geom, Mesh);

      deepCopy(TrCell, 1._Real);
      deepCopy(Kokkos::subview(TrCell, Kokkos::ALL, Kokkos::ALL, NVertLayers),
               0._Real);
      const int L = 0;
      parallelFor(
          {Mesh->NCellsAll, NVertLayers},
          KOKKOS_LAMBDA(const int ICell, const int K) {
             TrHorzAdvOnC.FCTTracerCurFill(L, ICell, K, TrCell);
          });
      parallelFor(
          {Mesh->NCellsOwned, NVertLayers},
          KOKKOS_LAMBDA(const int ICell, const int K) {
             TrHorzAdvOnC.FCTProvisionaLayerThicknesses(
                 ICell, K, Dt, FluxPseudoThickEdge, LayerThickness,
                 NormalVelocity);
          });
      parallelForOuter(
          {Mesh->NEdgesHaloH(1)},
          KOKKOS_LAMBDA(const int IEdge, const TeamMember &Team) {
             TrHorzAdvOnC.FCTHighAndLowOrderFlux(
                 Team, IEdge, FluxPseudoThickEdge, NormalVelocity);
          });
      parallelForOuter(
          {Mesh->NCellsHaloH(0)},
          KOKKOS_LAMBDA(const int ICell, const TeamMember &Team) {
             TrHorzAdvOnC.FCTFluxInOut(Team, ICell, Dt, LayerThickness);
          });

      Array2DReal FlxOut = TrHorzAdvOnC.GetFlxOut();
      Err += computeErrors(FCTErrors, FlxOut, CellSubView, Mesh, OnCell);
      Err += checkErrors("TendencyTermsTest", "FCTFluxOut", FCTErrors,
                         Setup.ExpectedFCTErrors, RTol);
      if (Err == 0)
         LOG_INFO("TendencyTermsTest: FCTFluxOut PASS");
      TotErr += Err;
      Err = 0;
   }
   {
      parallelFor(
          {Mesh->NEdgesHaloH(0), NVertLayers},
          KOKKOS_LAMBDA(const int IEdge, const int K) {
             TrHorzAdvOnC.FCTRescaleHighOrderFlux(IEdge, K);
          });
      const Array2DReal HighOrderFlx = TrHorzAdvOnC.GetHighOrderFlx();
      Err += computeErrors(FCTErrors, HighOrderFlx, ReferenceSolution, Mesh,
                           OnEdge);
      Err += checkErrors("TendencyTermsTest", "FCTRescaleHighOrderFlux",
                         FCTErrors, Setup.ExpectedFCTErrors, RTol);
      if (Err == 0)
         LOG_INFO("TendencyTermsTest: FCTRescaleHighOrderFlux PASS");
      TotErr += Err;
      Err = 0;
   }
   {
      setVectorEdge(
          KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
             VecField[0] = Setup.vectorX(X, Y);
             VecField[1] = Setup.vectorY(X, Y);
          },
          NormalVelocity, EdgeComponent::Normal, Geom, Mesh);

      deepCopy(TrCell, 1._Real);
      deepCopy(Kokkos::subview(TrCell, Kokkos::ALL, Kokkos::ALL, NVertLayers),
               0._Real);
      const int L = 0;
      parallelFor(
          {Mesh->NCellsAll, NVertLayers},
          KOKKOS_LAMBDA(const int ICell, const int K) {
             TrHorzAdvOnC.FCTTracerCurFill(L, ICell, K, TrCell);
          });
      parallelFor(
          {Mesh->NCellsOwned, NVertLayers},
          KOKKOS_LAMBDA(const int ICell, const int K) {
             TrHorzAdvOnC.FCTProvisionaLayerThicknesses(
                 ICell, K, Dt, FluxPseudoThickEdge, LayerThickness,
                 NormalVelocity);
          });
      parallelForOuter(
          {Mesh->NEdgesHaloH(1)},
          KOKKOS_LAMBDA(const int IEdge, const TeamMember &Team) {
             TrHorzAdvOnC.FCTHighAndLowOrderFlux(
                 Team, IEdge, FluxPseudoThickEdge, NormalVelocity);
          });
      parallelForOuter(
          {Mesh->NCellsHaloH(0)},
          KOKKOS_LAMBDA(const int ICell, const TeamMember &Team) {
             TrHorzAdvOnC.FCTFluxInOut(Team, ICell, Dt, LayerThickness);
          });
      parallelForOuter(
          {Mesh->NCellsOwned},
          KOKKOS_LAMBDA(const int ICell, const TeamMember &Team) {
             TrHorzAdvOnC.FCTAccumulateHighOrderFlux(Team, ICell, Dt, Tend,
                                                     LayerThickness);
          });

      const auto &MinLayerEdgeBot = VCoord->MinLayerEdgeBot;
      const auto &MaxLayerEdgeTop = VCoord->MaxLayerEdgeTop;
      const auto &MinLayerCell    = VCoord->MinLayerCell;
      const auto &MaxLayerCell    = VCoord->MaxLayerCell;
      parallelForOuter(
          "Tend:TracerHorzAdv.init",
          LaunchConfig({Mesh->NEdgesAll}, TeamScratch<Real>(NVertLayers)),
          KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             TrHorzAdvOnC(Team, L, IEdge, TrCell, FluxPseudoThickEdge,
                          NormalVelocity);
          });
      parallelForOuter(
          "Tend:TracerHorzAdv.exec",
          LaunchConfig({Mesh->NCellsAll}, TeamScratch<Real>(NVertLayers)),
          KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
             TrHorzAdvOnC(Team, TendNoFCT, L, ICell);
          });
      deepCopy(TracerSubView,
               Kokkos::subview(TendNoFCT, L, Kokkos::ALL, Kokkos::ALL));

      const Real ATol = sizeof(Real) == 4 ? 1e-6 : 1e-10;
      Err += computeErrors(FCTErrors, Tend, TracerSubView, Mesh, OnCell);
      Err += checkErrors("TendencyTermsTest",
                         "FCTAccumulateHighOrderFlux_" + std::to_string(L),
                         FCTErrors, Setup.ExpectedFCTToNonFCT, RTol, ATol);
      if (Err == 0)
         LOG_INFO("TendencyTermsTest: FCTAccumulateHighOrderFlux_" +
                  std::to_string(L) + " PASS");
      TotErr += Err;
      Err = 0;
   }
   VertAdv::clear();
   return TotErr;
} // end testTracerHorzAdvOnCell

int testTracerDiffOnCell(int NVertLayers, int NTracers, Real RTol) {

   I4 Err = 0;
   TestSetup Setup;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   const Real EddyDiff2 = 0.84_Real;

   // Compute exact result
   Array3DReal ExactTracerDiff("ExactTracerDiff", NTracers, Mesh->NCellsOwned,
                               NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) {
          return Setup.tracerDiff(X, Y, EddyDiff2);
       },
       ExactTracerDiff, Geom, Mesh, OnCell, ExchangeHalos::No);

   // Set input arrays
   Array3DReal TracerCell("TracerCell", NTracers, Mesh->NCellsSize,
                          NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.scalarA(X, Y); },
       TracerCell, Geom, Mesh, OnCell);

   Array2DReal PseudoThickEdge("PseudoThickEdge", Mesh->NEdgesSize,
                               NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.scalarB(X, Y); },
       PseudoThickEdge, Geom, Mesh, OnEdge);

   // Compute numerical result
   Array3DReal NumTracerDiff("NumTracerDiff", NTracers, Mesh->NCellsOwned,
                             NVertLayers);
   TracerDiffOnCell TrDiffOnC(Mesh, VCoord);
   TrDiffOnC.EddyDiff2 = EddyDiff2;

   parallelForOuter(
       LaunchConfig({NTracers, Mesh->NCellsOwned},
                    TeamScratch<Real>(NVertLayers)),
       KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
          TrDiffOnC(Team, NumTracerDiff, L, ICell, TracerCell, PseudoThickEdge);
       });

   ErrorMeasures TrDiffErrors;
   Err += computeErrors(TrDiffErrors, NumTracerDiff, ExactTracerDiff, Mesh,
                        OnCell);

   Err += checkErrors("TendencyTermsTest", "TracerDiff", TrDiffErrors,
                      Setup.ExpectedTrDel2Errors, RTol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: TracerDiff PASS");
   }

   return Err;
} // end testTracerDiffOnCell

int testKPPNonLocalTracerFluxOnCell(int NVertLayers, int NTracers, Real RTol) {

   I4 Err = 0;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   constexpr Real InitialTendency = 0.375_Real;
   constexpr Real FluxScale       = 0.25_Real;
   constexpr Real ProfileScale    = 0.01_Real;

   Array2DReal SurfaceTracerFlux("SurfaceTracerFlux", NTracers,
                                 Mesh->NCellsSize);
   Array2DReal VertNonLocalFlux("VertNonLocalFlux", Mesh->NCellsSize,
                                NVertLayers + 1);
   Array3DReal ExactTendency("ExactKPPNonLocalTracerTendency", NTracers,
                             Mesh->NCellsOwned, NVertLayers);
   Array3DReal NumTendency("NumKPPNonLocalTracerTendency", NTracers,
                           Mesh->NCellsOwned, NVertLayers);

   deepCopy(ExactTendency, InitialTendency);
   deepCopy(NumTendency, InitialTendency);

   parallelFor(
       {NTracers, Mesh->NCellsSize}, KOKKOS_LAMBDA(int L, int ICell) {
          SurfaceTracerFlux(L, ICell) = FluxScale * Real(L + 1);
       });
   parallelFor(
       {Mesh->NCellsSize, NVertLayers + 1}, KOKKOS_LAMBDA(int ICell, int K) {
          const Real DistanceFromBottom = Real(NVertLayers - K);
          VertNonLocalFlux(ICell, K)    = ProfileScale * Real(ICell + 1) *
                                       DistanceFromBottom * DistanceFromBottom;
       });

   const auto MinLayerCell = VCoord->MinLayerCell;
   const auto MaxLayerCell = VCoord->MaxLayerCell;
   parallelForOuter(
       {NTracers, Mesh->NCellsOwned},
       KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
          const int KMin = MinLayerCell(ICell);
          const int KMax = MaxLayerCell(ICell);
          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                 const Real DistanceFromBottom = Real(NVertLayers - K);
                 ExactTendency(L, ICell, K) +=
                     FluxScale * Real(L + 1) * ProfileScale * Real(ICell + 1) *
                     (2.0_Real * DistanceFromBottom - 1.0_Real);
              });
       });

   KPPNonLocalTracerFluxOnCell KPPNonLocalFlux(Mesh, VCoord);
   parallelForOuter(
       {NTracers, Mesh->NCellsOwned},
       KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
          KPPNonLocalFlux(Team, NumTendency, L, ICell, SurfaceTracerFlux,
                          VertNonLocalFlux);
       });

   ErrorMeasures KPPNonLocalErrors;
   Err += computeErrors(KPPNonLocalErrors, NumTendency, ExactTendency, Mesh,
                        OnCell);

   const ErrorMeasures ExpectedErrors{0.0_Real, 0.0_Real};
   const Real ATol = 100 * std::numeric_limits<Real>::epsilon();
   Err += checkErrors("TendencyTermsTest", "KPPNonLocalTracerFlux",
                      KPPNonLocalErrors, ExpectedErrors, RTol, ATol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: KPPNonLocalTracerFluxOnCell PASS");
   }

   return Err;
} // end testKPPNonLocalTracerFluxOnCell

int testTracerHyperDiffOnCell(int NVertLayers, int NTracers, Real RTol) {

   I4 Err = 0;
   TestSetup Setup;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   const Real EddyDiff4 = 1.23_Real;

   // Compute exact result
   Array3DReal ExactTracerHyperDiff("ExactTracerHyperDiff", NTracers,
                                    Mesh->NCellsOwned, NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) {
          return -Setup.tracerHyperDiff(X, Y, EddyDiff4);
       },
       ExactTracerHyperDiff, Geom, Mesh, OnCell, ExchangeHalos::No);

   // Set input arrays
   Array3DReal TrDel2Cell("TracerCell", NTracers, Mesh->NCellsSize,
                          NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.scalarC(X, Y); },
       TrDel2Cell, Geom, Mesh, OnCell);

   // Compute numerical result
   Array3DReal NumTracerHyperDiff("NumTracerHyperDiff", NTracers,
                                  Mesh->NCellsOwned, NVertLayers);
   TracerHyperDiffOnCell TrHypDiffOnC(Mesh, VCoord);
   TrHypDiffOnC.EddyDiff4 = EddyDiff4;
   parallelForOuter(
       LaunchConfig({NTracers, Mesh->NCellsOwned},
                    TeamScratch<Real>(NVertLayers)),
       KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
          TrHypDiffOnC(Team, NumTracerHyperDiff, L, ICell, TrDel2Cell);
       });

   ErrorMeasures TrHyperDiffErrors;
   Err += computeErrors(TrHyperDiffErrors, NumTracerHyperDiff,
                        ExactTracerHyperDiff, Mesh, OnCell);

   Err += checkErrors("TendencyTermsTest", "TracerHyperDiff", TrHyperDiffErrors,
                      Setup.ExpectedTrDel4Errors, RTol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: TracerHyperDiff PASS");
   }

   return Err;
} // end testTracerHyperDiffOnCell

int testSurfaceTracerRestoringOnCell(int NVertLayers, int NTracers, Real RTol) {

   I4 Err = 0;
   TestSetup Setup;

   const auto Mesh = HorzMesh::getDefault();

   if (NTracers < 2) {
      LOG_ERROR("TendencyTermsTest: SurfaceTracerRestoring requires at least 2 "
                "tracers, found {}",
                NTracers);
      return 1;
   }

   // Test multiple cases with different combinations of tracers being restored
   const char *CaseLabels[3] = {"SurfaceTracerRestoringSalinityOnly",
                                "SurfaceTracerRestoringTemperatureOnly",
                                "SurfaceTracerRestoringTempSalt"};

   const std::vector<std::vector<I4>> CaseTracerIds = {
       {1},    // Salinity Only
       {0},    // Temperature Only
       {0, 1}, // Temperature and Salinity
   };

   // Loop over cases, computing exact result, numerical result, and
   // errors for each.
   for (int Case = 0; Case < 3; ++Case) {

      Array3DReal ExactSurfRest("ExactSurfRest", NTracers, Mesh->NCellsOwned,
                                NVertLayers);
      Array2DReal InputField("InputField", Mesh->NCellsSize, 1);
      Array2DReal TracersMonthlySurfClimoCell("TracersMonthlySurfClimoCell",
                                              NTracers, Mesh->NCellsSize);
      Array3DReal TracersOnCell("TracersOnCell", NTracers, Mesh->NCellsSize,
                                NVertLayers);
      Array3DReal NumSurfRest("NumSurfRest", NTracers, Mesh->NCellsOwned,
                              NVertLayers);

      deepCopy(ExactSurfRest, 0);
      deepCopy(InputField, 0);
      deepCopy(TracersMonthlySurfClimoCell, 0);
      deepCopy(TracersOnCell, 0);
      deepCopy(NumSurfRest, 0);

      // Set Input Field values. Use a combination of scalarB and vectorX to
      // ensure the full surface tracer restoring logic is exercised.
      Err += setScalar(
          KOKKOS_LAMBDA(Real X, Real Y) {
             return Setup.scalarB(X, Y) + Setup.vectorX(X, Y);
          },
          InputField, Geom, Mesh, OnCell);

      // Set TracersOnCell values (use scalarB for simplicity, but could be
      // any field).
      Err += setScalar(
          KOKKOS_LAMBDA(Real X, Real Y) { return Setup.scalarB(X, Y); },
          TracersOnCell, Geom, Mesh, OnCell);
      parallelFor(
          {NTracers, Mesh->NCellsSize, NVertLayers},
          KOKKOS_LAMBDA(int L, int ICell, int K) {
             TracersOnCell(L, ICell, K) += 0.04_Real * K;
          });

      SurfaceTracerRestoringOnCell SurfRestOnC(Mesh);
      SurfRestOnC.PistonVelocity = 1.585e-5;

      // Build host-selected tracer IDs for restoring and copy to device.
      const I4 NTracersToRestore = static_cast<I4>(CaseTracerIds[Case].size());
      Array1DI4 TracerIdsToRestore("TracerIdsToRestore", NTracersToRestore);
      HostArray1DI4 TracerIdsToRestoreH("TracerIdsToRestoreH",
                                        NTracersToRestore);
      for (I4 R = 0; R < NTracersToRestore; ++R) {
         TracerIdsToRestoreH(R) = CaseTracerIds[Case][R];
      }
      deepCopy(TracerIdsToRestore, TracerIdsToRestoreH);

      // Compute exact result using the same logic as
      // SurfaceTracerRestoringOnCell, but iterating
      // selected tracer IDs to match restoring implementation.
      parallelFor(
          {NTracersToRestore, Mesh->NCellsOwned},
          KOKKOS_LAMBDA(int R, int ICell) {
             const int L                           = TracerIdsToRestore(R);
             TracersMonthlySurfClimoCell(L, ICell) = InputField(ICell, 0);
             const Real Diff = TracersMonthlySurfClimoCell(L, ICell) -
                               TracersOnCell(L, ICell, 0);

             ExactSurfRest(L, ICell, 0) = SurfRestOnC.PistonVelocity * Diff;
          });

      // Compute numerical result
      parallelFor(
          {NTracersToRestore, Mesh->NCellsOwned},
          KOKKOS_LAMBDA(int R, int ICell) {
             const int L = TracerIdsToRestore(R);
             SurfRestOnC(NumSurfRest, L, ICell, 0, TracersMonthlySurfClimoCell,
                         TracersOnCell);
          });

      ErrorMeasures SurfRestErrors;
      Err += computeErrors(SurfRestErrors, NumSurfRest, ExactSurfRest, Mesh,
                           OnCell);

      Err += checkErrors("TendencyTermsTest", CaseLabels[Case], SurfRestErrors,
                         Setup.ExpectedSurfTrRestErrors, RTol);
   }

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: SurfaceTracerRestoring PASS");
   }

   return Err;
} // end testSurfaceTracerRestoringOnCell

int testPenetratingShortwaveOnCell(int NVertLayers, Real RTol) {

   I4 Err = 0;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   // Test values for a shallow water column where shortwave radiation
   // penetrates all the way to the seafloor with significant leftover energy. A
   // layer thickness of 1.0 m (16 m total depth) with realistic blue and red
   // extinction coefficients ensures that substantial flux reaches the bottom
   // interface. The test verifies that this leftover flux is fully absorbed
   // into the bottom layer so that column-integrated heating is conserved
   // (integrated fraction is 1.0).
   const Real SurfaceFluxVal    = 235.0_Real; // W/m^2
   const Real ExtinctionRedVal  = 0.2_Real;   // 1/m
   const Real ExtinctionBlueVal = 0.05_Real;  // 1/m
   const Real NearIrFractionVal = 0.58_Real;
   const Real NearIrCoeffVal    = 2.86_Real; // 1/m
   const Real RedFractionVal    = 0.21_Real;
   const Real BlueFractionVal   = 0.21_Real;
   const Real LayerThickness    = 1.0_Real; // m

   Array1DReal ShortWaveHeatFlux("ShortWaveHeatFlux", Mesh->NCellsSize);
   Array1DReal ExtinctionCoeffRed("ExtinctionCoeffRed", Mesh->NCellsSize);
   Array1DReal ExtinctionCoeffBlue("ExtinctionCoeffBlue", Mesh->NCellsSize);
   deepCopy(ShortWaveHeatFlux, SurfaceFluxVal);
   deepCopy(ExtinctionCoeffRed, ExtinctionRedVal);
   deepCopy(ExtinctionCoeffBlue, ExtinctionBlueVal);

   // An independent, exactly-known column geometry (uniform layer
   // thickness), rather than relying on the mesh's default vertical
   // coordinate, so the expected result is known exactly.
   Array2DReal GeomZInterface("GeomZInterface", Mesh->NCellsSize,
                              NVertLayers + 1);
   parallelFor(
       {Mesh->NCellsSize, NVertLayers + 1}, KOKKOS_LAMBDA(int ICell, int K) {
          GeomZInterface(ICell, K) = -LayerThickness * K;
       });

   Array3DReal Tend("Tend", 1, Mesh->NCellsSize, NVertLayers);
   deepCopy(Tend, 0);

   PenetratingShortwaveOnCell PenSwOnC(Mesh, VCoord, /*TempTracerIndex=*/0);
   PenSwOnC.Enabled        = true;
   PenSwOnC.NearIrFraction = NearIrFractionVal;
   PenSwOnC.NearIrCoeff    = NearIrCoeffVal;
   PenSwOnC.RedFraction    = RedFractionVal;
   PenSwOnC.BlueFraction   = BlueFractionVal;
   const auto MinLayerCell = VCoord->MinLayerCell;
   const auto MaxLayerCell = VCoord->MaxLayerCell;

   parallelForOuter(
       LaunchConfig({Mesh->NCellsOwned}, TeamScratch<Real>(NVertLayers + 1)),
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          PenSwOnC(Team, Tend, ICell, GeomZInterface, ShortWaveHeatFlux,
                   ExtinctionCoeffRed, ExtinctionCoeffBlue);
       });

   I4 NumBad                        = 0;
   const Real ExpectedColumnHeating = SurfaceFluxVal * HFluxFac;
   parallelReduce(
       {Mesh->NCellsOwned},
       KOKKOS_LAMBDA(int ICell, I4 &Accum) {
          const I4 KTop = MinLayerCell(ICell);
          const I4 KBot = MaxLayerCell(ICell);
          if (KTop > KBot) {
             return;
          }

          // Analytic flux that would reach the seafloor interface without
          // bottom-layer absorption
          const Real SeafloorDepth = Kokkos::abs(
              GeomZInterface(ICell, KBot + 1) - GeomZInterface(ICell, KTop));
          const Real ResidualSeafloorFlux =
              SurfaceFluxVal *
              (NearIrFractionVal *
                   Kokkos::exp(-NearIrCoeffVal * SeafloorDepth) +
               RedFractionVal * Kokkos::exp(-ExtinctionRedVal * SeafloorDepth) +
               BlueFractionVal *
                   Kokkos::exp(-ExtinctionBlueVal * SeafloorDepth));

          // Ensure the test setup genuinely has significant leftover flux at
          // the bottom interface (at least 5% of incident flux)
          if (ResidualSeafloorFlux < 0.05_Real * SurfaceFluxVal) {
             Accum += 1;
          }

          Real FluxAtLayerTop = SurfaceFluxVal;
          Real ColumnHeating  = 0.0_Real;
          for (I4 K = KTop; K <= KBot; ++K) {
             Real FluxAtLayerBottom = 0.0_Real;
             if (K < KBot) {
                const Real Depth = Kokkos::abs(GeomZInterface(ICell, K + 1) -
                                               GeomZInterface(ICell, KTop));
                FluxAtLayerBottom =
                    SurfaceFluxVal *
                    (NearIrFractionVal * Kokkos::exp(-NearIrCoeffVal * Depth) +
                     RedFractionVal * Kokkos::exp(-ExtinctionRedVal * Depth) +
                     BlueFractionVal * Kokkos::exp(-ExtinctionBlueVal * Depth));
             }
             const Real ExpectedLayerHeating =
                 (FluxAtLayerTop - FluxAtLayerBottom) * HFluxFac;
             ColumnHeating += Tend(0, ICell, K);
             const Real RelErr =
                 Kokkos::abs(Tend(0, ICell, K) - ExpectedLayerHeating) /
                 Kokkos::abs(ExpectedLayerHeating);
             if (RelErr > RTol || Kokkos::isnan(Tend(0, ICell, K)) ||
                 Kokkos::isinf(Tend(0, ICell, K))) {
                Accum += 1;
             }

             // In the bottom layer, verify that the absorbed heating is
             // strictly greater than unadjusted exponential divergence,
             // confirming the leftover seafloor flux was absorbed in the bottom
             // layer
             if (K == KBot) {
                const Real UnadjustedBottomHeating =
                    (FluxAtLayerTop - ResidualSeafloorFlux) * HFluxFac;
                if (Tend(0, ICell, K) <= UnadjustedBottomHeating) {
                   Accum += 1;
                }
             }

             FluxAtLayerTop = FluxAtLayerBottom;
          }
          const Real ColumnRelErr =
              Kokkos::abs(ColumnHeating - ExpectedColumnHeating) /
              Kokkos::abs(ExpectedColumnHeating);
          if (ColumnRelErr > RTol || Kokkos::isnan(ColumnHeating) ||
              Kokkos::isinf(ColumnHeating)) {
             Accum += 1;
          }
       },
       NumBad);

   if (NumBad > 0) {
      LOG_ERROR("TendencyTermsTest: PenetratingShortwave FAIL, {} layer or "
                "column heating values do not match the expected results",
                NumBad);
      Err += 1;
   } else {
      LOG_INFO("TendencyTermsTest: PenetratingShortwave PASS");
   }

   return Err;
} // end testPenetratingShortwaveOnCell

void initTendTest(const std::string &MeshFile, int NVertLayers) {

   Error Err;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // Initialize logging
   initLogging(DefEnv);

   // Open config file
   Config::Initialize();
   Config::readAll("omega.yml");

   // The analytic solutions assume constant viscosity and diffusivity, so
   // turn off the mesh scaling
   Config *OmegaConfig = Config::getOmegaConfig();
   Config HorzMeshConfig("HorzMesh");
   Err = OmegaConfig->get(HorzMeshConfig);
   CHECK_ERROR_ABORT(Err,
                     "TendencyTermsTest: HorzMesh group not found in Config");
   HorzMeshConfig.set("MeshScalingEnable", false);

   // Initialize time stepping and get model clock
   TimeStepper::init1();
   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   IO::init(DefComm);

   Decomp::init(MeshFile);

   int HaloErr = Halo::init();
   if (HaloErr != 0) {
      ABORT_ERROR("TendencyTermsTest: error initializing default halo");
   }

   Field::init(ModelClock);    // Fields are used in streams
   IOStream::init(ModelClock); // initialize streams for mesh reading
   HorzMesh::init(ModelClock);

   // initialize vertical coordinate, do not read stream and use local
   // NVertLayers value
   VertCoord::init(false, NVertLayers);
   Tracers::init();

} // end initTendTest

void finalizeTendTest() {
   Tracers::clear();
   HorzMesh::clear();
   VertCoord::clear();
   Field::clear();
   Dimension::clear();
   TimeStepper::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();

} // end finalizeTendTest

int tendencyTermsTest(const std::string &MeshFile = DefaultMeshFile) {
   int Err         = 0;
   int NVertLayers = 16;

   initTendTest(MeshFile, NVertLayers);

   int NTracers = Tracers::getNumTracers();

   const Real RTol = sizeof(Real) == 4 ? 2e-2 : 1e-5;

   Err += testThickFluxDiv(NVertLayers, RTol);

   Err += testPotVortHAdv(NVertLayers, RTol);

   Err += testCoriolisAccelerationOnEdge(NVertLayers, RTol);

   Err += testKEGrad(NVertLayers, RTol);

   Err += testSSHGrad(NVertLayers, RTol);

   Err += testVelDiff(NVertLayers, RTol);

   Err += testVelHyperDiff(NVertLayers, RTol);

   Err += testSfcStressForcing(NVertLayers);

   Err += testExplicitBottomDrag(NVertLayers, RTol);

   Err += testImplicitBottomDrag(NVertLayers, RTol);

   Err += testBottomDragInactiveEdges(NVertLayers);

   Err += testTracerHorzAdvOnCell(NVertLayers, NTracers, RTol);

   Err += testFCTTracerHorzAdvOnCell(NVertLayers, NTracers, RTol);

   Err += testTracerDiffOnCell(NVertLayers, NTracers, RTol);

   Err += testKPPNonLocalTracerFluxOnCell(NVertLayers, NTracers, RTol);

   Err += testTracerHyperDiffOnCell(NVertLayers, NTracers, RTol);

   Err += testSurfaceTracerRestoringOnCell(NVertLayers, NTracers, RTol);

   Err += testPenetratingShortwaveOnCell(NVertLayers, RTol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: Successful completion");
   }

   finalizeTendTest();

   return Err;

} // end tendencyTermsTest

int main(int argc, char *argv[]) {

   int RetErr = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   RetErr = tendencyTermsTest();

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return RetErr;

} // end of main
//===-----------------------------------------------------------------------===/
