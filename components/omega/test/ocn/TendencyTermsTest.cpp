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
#include "mpi.h"

#include <cmath>
#include <limits>
#include <vector>

using namespace OMEGA;

struct TestSetupPlane {

   Real Lx = 1;
   Real Ly = SqrtThree / 2;

   ErrorMeasures ExpectedDivErrors         = {0.00124886886594453264,
                                              0.00124886886590977139};
   ErrorMeasures ExpectedPVErrors          = {0.00807347170900282914,
                                              0.00794755105765788429};
   ErrorMeasures ExpectedGradErrors        = {0.00125026071878537952,
                                              0.00134354611117262161};
   ErrorMeasures ExpectedLaplaceErrors     = {0.00113090174765822192,
                                              0.00134324628763667899};
   ErrorMeasures ExpectedTrHAdvErrors      = {0.0029211089892916243,
                                              0.0024583038518548855};
   ErrorMeasures ExpectedTrDel2Errors      = {0.00334357193650093847,
                                              0.00290978146207349032};
   ErrorMeasures ExpectedTrDel4Errors      = {0.00508833446725232875,
                                              0.00523080740758275625};
   ErrorMeasures ExpectedSurfTrRestErrors  = {0, 0};
   ErrorMeasures ExpectedWindForcingErrors = {0, 0};
   ErrorMeasures ExpectedBottomDragErrors  = {0.033848740052302935,
                                              0.01000133508329411};

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

   KOKKOS_FUNCTION Real tracerDiff(Real X, Real Y) const {
      return -TwoPi * TwoPi * std::sin(TwoPi * Y / Ly) *
             (2 * (1 / Lx / Lx + 1 / Ly / Ly) * std::cos(TwoPi * X / Lx) +
              (1 / Ly / Ly +
               (1 / Lx / Lx + 1 / Ly / Ly) * std::cos(2 * TwoPi * X / Lx)) *
                  std::cos(TwoPi * Y / Ly));
   }

   KOKKOS_FUNCTION Real scalarC(Real X, Real Y) const {
      return std::pow(std::cos(TwoPi * X / Lx), 2) -
             std::pow(std::sin(TwoPi * Y / Ly), 2);
   }

   KOKKOS_FUNCTION Real tracerHyperDiff(Real X, Real Y) const {
      return -2 * TwoPi * TwoPi *
             (std::cos(2 * TwoPi * X / Lx) / Lx / Lx +
              std::cos(2 * TwoPi * Y / Ly) / Ly / Ly);
   }

   KOKKOS_FUNCTION Real windForcingX(Real X, Real Y) const {
      const Real StressU = vectorX(X, Y);
      const Real Thick   = scalarB(X, Y);
      return StressU / (Thick * RhoSw);
   }

   KOKKOS_FUNCTION Real windForcingY(Real X, Real Y) const {
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
   // TODO: get this from the mesh
   Real Radius = REarth;

   ErrorMeasures ExpectedDivErrors         = {0.013652414501664885,
                                              0.0036904315983599676};
   ErrorMeasures ExpectedPVErrors          = {0.0219217796608757037,
                                              0.0122537418367830303};
   ErrorMeasures ExpectedGradErrors        = {0.0019094381714837498,
                                              0.0015218320661105702};
   ErrorMeasures ExpectedLaplaceErrors     = {0.28193638497826856,
                                              0.270546491554748};
   ErrorMeasures ExpectedTrHAdvErrors      = {0.013259410329645643,
                                              0.004094907022292395};
   ErrorMeasures ExpectedTrDel2Errors      = {0.04865718541236144,
                                              0.005105510870642706};
   ErrorMeasures ExpectedTrDel4Errors      = {0.0008646345116716073,
                                              0.0007118574326665881};
   ErrorMeasures ExpectedSurfTrRestErrors  = {0, 0};
   ErrorMeasures ExpectedWindForcingErrors = {0, 0};
   ErrorMeasures ExpectedBottomDragErrors  = {0.0015333449035655053,
                                              0.0014897009917655022};

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

   KOKKOS_FUNCTION Real scalarB(Real Lon, Real Lat) const {
      return 2. + std::cos(Lon) * std::sin(Lat);
   }

   KOKKOS_FUNCTION Real tracerDiff(Real Lon, Real Lat) const {
      return (4 * std::pow(std::cos(Lon), 2) -
              2 * (1. + 3 * std::cos(2 * Lat)) * std::pow(std::sin(Lon), 2) +
              2 * std::pow(std::cos(Lon), 3) * std::sin(Lat) -
              8 * std::cos(Lon) * std::pow(std::cos(Lat), 2) *
                  std::pow(std::sin(Lon), 2) * std::sin(Lat)) /
             Radius;
   }

   KOKKOS_FUNCTION Real scalarC(Real Lon, Real Lat) const {
      return -(Radius / 2) * std::sqrt(3. / 2. / Pi) * std::cos(Lat) *
             std::cos(Lon);
   }

   KOKKOS_FUNCTION Real tracerHyperDiff(Real Lon, Real Lat) const {
      return std::sqrt(3. / 2. / Pi) * std::cos(Lat) * std::cos(Lon) / Radius;
   }

   KOKKOS_FUNCTION Real windForcingX(Real Lon, Real Lat) const {
      const Real StressU = vectorX(Lon, Lat);
      const Real Thick   = scalarB(Lon, Lat);
      return StressU / (Thick * RhoSw);
   }

   KOKKOS_FUNCTION Real windForcingY(Real Lon, Real Lat) const {
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
   parallelFor(
       {Mesh->NCellsOwned, NVertLayers}, KOKKOS_LAMBDA(int ICell, int KLayer) {
          ThickFluxDivOnC(NumThickFluxDiv, ICell, KLayer, OnesEdge,
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
   parallelFor(
       {Mesh->NEdgesOwned, NVertLayers}, KOKKOS_LAMBDA(int IEdge, int KLayer) {
          PotVortHAdvOnE(NumPotVortHAdv, IEdge, KLayer, NormRelVortEdge,
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
   parallelFor(
       {Mesh->NEdgesOwned, NVertLayers}, KOKKOS_LAMBDA(int IEdge, int KLayer) {
          KEGradOnE(NumKEGrad, IEdge, KLayer, KECell);
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
   parallelFor(
       {Mesh->NEdgesOwned, NVertLayers}, KOKKOS_LAMBDA(int IEdge, int KLayer) {
          SSHGradOnE(NumSSHGrad, IEdge, KLayer, SSHCell);
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

   parallelFor(
       {Mesh->NEdgesOwned, NVertLayers}, KOKKOS_LAMBDA(int IEdge, int KLayer) {
          VelDiffOnE(NumVelDiff, IEdge, KLayer, DivCell, RVortVertex);
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

   parallelFor(
       {Mesh->NEdgesOwned, NVertLayers}, KOKKOS_LAMBDA(int IEdge, int KLayer) {
          VelHyperDiffOnE(NumVelHyperDiff, IEdge, KLayer, DivCell, RVortVertex);
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

int testWindForcing(int NVertLayers) {

   int Err = 0;
   TestSetup Setup;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   // Compute exact result
   Array2DReal ExactWindForcing("ExactWindForcing", Mesh->NEdgesOwned,
                                NVertLayers);

   // Note: this computes wind forcing at every layer
   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.windForcingX(X, Y);
          VecField[1] = Setup.windForcingY(X, Y);
       },
       ExactWindForcing, EdgeComponent::Normal, Geom, Mesh, ExchangeHalos::No);

   // Reset wind forcing to zero below the surface
   deepCopy(Kokkos::subview(ExactWindForcing, Kokkos::ALL,
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
   Array2DReal NumWindForcing("NumWindForcing", Mesh->NEdgesOwned, NVertLayers);

   WindForcingOnEdge WindForcingOnE(Mesh, VCoord);

   parallelFor(
       {Mesh->NEdgesOwned, NVertLayers}, KOKKOS_LAMBDA(int IEdge, int KLayer) {
          WindForcingOnE(NumWindForcing, IEdge, KLayer, NormalStressEdge,
                         PseudoThickEdge);
       });

   // Compute errors
   ErrorMeasures WindForcingErrors;
   Err += computeErrors(WindForcingErrors, NumWindForcing, ExactWindForcing,
                        Mesh, OnEdge);

   // Check error values
   const Real RTol = 0;
   const Real ATol = 100 * std::numeric_limits<Real>::epsilon();
   Err += checkErrors("TendencyTermsTest", "WindForcing", WindForcingErrors,
                      Setup.ExpectedWindForcingErrors, RTol, ATol);

   return Err;
} // end testWindForcing

int testBottomDrag(int NVertLayers, Real RTol) {

   int Err = 0;
   TestSetup Setup;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   const Real Coeff = 1.123456789;

   // Compute exact result
   Array2DReal ExactBottomDrag("ExactBottomDrag", Mesh->NEdgesOwned,
                               NVertLayers);

   // Note: this computes bottom drag at every layer
   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.bottomDragX(X, Y, Coeff);
          VecField[1] = Setup.bottomDragY(X, Y, Coeff);
       },
       ExactBottomDrag, EdgeComponent::Normal, Geom, Mesh, ExchangeHalos::No);

   // Reset bottom drag to zero above the lowest layer
   deepCopy(Kokkos::subview(ExactBottomDrag, Kokkos::ALL,
                            Kokkos::make_pair(0, NVertLayers - 1)),
            0);

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
   Err += checkErrors("TendencyTermsTest", "BottomDrag", BottomDragErrors,
                      Setup.ExpectedBottomDragErrors, RTol);

   return Err;
} // end testBottomDrag

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
   TracerHorzAdvOnCell TrHorzAdvOnC(Mesh, VCoord);
   TrHorzAdvOnC.init();
   TrHorzAdvOnC.ForceLowOrder = true;

   parallelFor(
       {NTracers, Mesh->NEdgesAll, NVertLayers},
       KOKKOS_LAMBDA(int L, int IEdge, int KLayer) {
          TrHorzAdvOnC(L, IEdge, KLayer, TrCell, ThickEdge, NormalVelocity);
       });

   parallelFor(
       {NTracers, Mesh->NCellsOwned, NVertLayers},
       KOKKOS_LAMBDA(int L, int ICell, int KLayer) {
          TrHorzAdvOnC(NumTrFluxDiv, L, ICell, KLayer);
       });

   ErrorMeasures TrHAdvErrors;
   Err +=
       computeErrors(TrHAdvErrors, NumTrFluxDiv, ExactTrFluxDiv, Mesh, OnCell);

   Err += checkErrors("TendencyTermsTest", "TracerHorzAdv", TrHAdvErrors,
                      Setup.ExpectedTrHAdvErrors, RTol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: TracerHorzAdv PASS");
   }

   return Err;
} // end testTracerHorzAdvOnCell

int testTracerDiffOnCell(int NVertLayers, int NTracers, Real RTol) {

   I4 Err = 0;
   TestSetup Setup;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   // Compute exact result
   Array3DReal ExactTracerDiff("ExactTracerDiff", NTracers, Mesh->NCellsOwned,
                               NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.tracerDiff(X, Y); },
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
   TrDiffOnC.EddyDiff2 = 1._Real;

   parallelFor(
       {NTracers, Mesh->NCellsOwned, NVertLayers},
       KOKKOS_LAMBDA(int L, int ICell, int KLayer) {
          TrDiffOnC(NumTracerDiff, L, ICell, KLayer, TracerCell,
                    PseudoThickEdge);
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

int testTracerHyperDiffOnCell(int NVertLayers, int NTracers, Real RTol) {

   I4 Err = 0;
   TestSetup Setup;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   // Compute exact result
   Array3DReal ExactTracerHyperDiff("ExactTracerHyperDiff", NTracers,
                                    Mesh->NCellsOwned, NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return -Setup.tracerHyperDiff(X, Y); },
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
   TrHypDiffOnC.EddyDiff4 = 1._Real;
   parallelFor(
       {NTracers, Mesh->NCellsOwned, NVertLayers},
       KOKKOS_LAMBDA(int L, int ICell, int KLayer) {
          TrHypDiffOnC(NumTracerHyperDiff, L, ICell, KLayer, TrDel2Cell);
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

   if (NTracers < 3) {
      LOG_ERROR("TendencyTermsTest: SurfaceTracerRestoring requires at least 3 "
                "tracers, found {}",
                NTracers);
      return 1;
   }

   // Test multiple cases with different combinations of tracers being restored
   const char *CaseLabels[3] = {"SurfaceTracerRestoringSalinityOnly",
                                "SurfaceTracerRestoringTemperatureOnly",
                                "SurfaceTracerRestoringTempSaltDebug"};

   const std::vector<std::vector<I4>> CaseTracerIds = {
       {1},       // Salinity Only
       {0},       // Temperature Only
       {0, 1, 2}, // Temperature, Salinity, and a Debug Tracer
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

   Err += testKEGrad(NVertLayers, RTol);

   Err += testSSHGrad(NVertLayers, RTol);

   Err += testVelDiff(NVertLayers, RTol);

   Err += testVelHyperDiff(NVertLayers, RTol);

   Err += testWindForcing(NVertLayers);

   Err += testBottomDrag(NVertLayers, RTol);

   Err += testTracerHorzAdvOnCell(NVertLayers, NTracers, RTol);

   Err += testTracerDiffOnCell(NVertLayers, NTracers, RTol);

   Err += testTracerHyperDiffOnCell(NVertLayers, NTracers, RTol);

   Err += testSurfaceTracerRestoringOnCell(NVertLayers, NTracers, RTol);

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
   MPI_Finalize();

   return RetErr;

} // end of main
//===-----------------------------------------------------------------------===/
