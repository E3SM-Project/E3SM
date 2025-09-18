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
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Error.h"
#include "GlobalConstants.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "Logging.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "VertCoord.h"
#include "mpi.h"

#include <cmath>
#include <limits>

using namespace OMEGA;

struct TestSetupPlane {

   Real Lx = 1;
   Real Ly = std::sqrt(3) / 2;

   ErrorMeasures ExpectedDivErrors         = {0.00124886886594453264,
                                              0.00124886886590977139};
   ErrorMeasures ExpectedPVErrors          = {0.00807347170900282914,
                                              0.00794755105765788429};
   ErrorMeasures ExpectedGradErrors        = {0.00125026071878537952,
                                              0.00134354611117262161};
   ErrorMeasures ExpectedLaplaceErrors     = {0.00113090174765822192,
                                              0.00134324628763667899};
   ErrorMeasures ExpectedTrHAdvErrors      = {0.00205864372747571571,
                                              0.00172418025417940784};
   ErrorMeasures ExpectedTrDel2Errors      = {0.00334357193650093847,
                                              0.00290978146207349032};
   ErrorMeasures ExpectedTrDel4Errors      = {0.00508833446725232875,
                                              0.00523080740758275625};
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

   KOKKOS_FUNCTION Real layerThick(Real X, Real Y) const {
      return 2. + std::sin(TwoPi * X / Lx) * std::cos(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real planetaryVort(Real X, Real Y) const {
      return std::cos(TwoPi * X / Lx) * std::cos(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real normRelVort(Real X, Real Y) const {
      return curl(X, Y) / layerThick(X, Y);
   }

   KOKKOS_FUNCTION Real normPlanetVort(Real X, Real Y) const {
      return planetaryVort(X, Y) / layerThick(X, Y);
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

   KOKKOS_FUNCTION Real windForcingX(Real X, Real Y,
                                     Real SaltWaterDensity) const {
      const Real StressU = vectorX(X, Y);
      const Real Thick   = scalarB(X, Y);
      return StressU / (Thick * SaltWaterDensity);
   }

   KOKKOS_FUNCTION Real windForcingY(Real X, Real Y,
                                     Real SaltWaterDensity) const {
      const Real StressV = vectorY(X, Y);
      const Real Thick   = scalarB(X, Y);
      return StressV / (Thick * SaltWaterDensity);
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

   ErrorMeasures ExpectedDivErrors         = {0.0136595773989796766,
                                              0.00367052484586384131};
   ErrorMeasures ExpectedPVErrors          = {0.0219217796608757037,
                                              0.0122537418367830303};
   ErrorMeasures ExpectedGradErrors        = {0.00187912292540623471,
                                              0.00149841802817334935};
   ErrorMeasures ExpectedLaplaceErrors     = {0.281930203304510130,
                                              0.270530313560271740};
   ErrorMeasures ExpectedTrHAdvErrors      = {0.0132310202299444034,
                                              0.0038523368564029538};
   ErrorMeasures ExpectedTrDel2Errors      = {0.0486107109846934185,
                                              0.00507514214194892694};
   ErrorMeasures ExpectedTrDel4Errors      = {0.000819552466009620408,
                                              0.00064700084412871962};
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

   KOKKOS_FUNCTION Real layerThick(Real Lon, Real Lat) const {
      return (2 + std::cos(Lon) * std::pow(std::cos(Lat), 4));
   }

   KOKKOS_FUNCTION Real planetaryVort(Real Lon, Real Lat) const {
      return std::sin(Lat);
   }

   KOKKOS_FUNCTION Real normRelVort(Real Lon, Real Lat) const {
      return curl(Lon, Lat) / layerThick(Lon, Lat);
   }

   KOKKOS_FUNCTION Real normPlanetVort(Real Lon, Real Lat) const {
      return planetaryVort(Lon, Lat) / layerThick(Lon, Lat);
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
      return -(Radius / 2) * std::sqrt(3 / 2 / Pi) * std::cos(Lat) *
             std::cos(Lon);
   }

   KOKKOS_FUNCTION Real tracerHyperDiff(Real Lon, Real Lat) const {
      return std::sqrt(3 / 2 / Pi) * std::cos(Lat) * std::cos(Lon) / Radius;
   }

   KOKKOS_FUNCTION Real windForcingX(Real Lon, Real Lat,
                                     Real SaltWaterDensity) const {
      const Real StressU = vectorX(Lon, Lat);
      const Real Thick   = scalarB(Lon, Lat);
      return StressU / (Thick * SaltWaterDensity);
   }

   KOKKOS_FUNCTION Real windForcingY(Real Lon, Real Lat,
                                     Real SaltWaterDensity) const {
      const Real StressV = vectorY(Lon, Lat);
      const Real Thick   = scalarB(Lon, Lat);
      return StressV / (Thick * SaltWaterDensity);
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

   const auto Mesh = HorzMesh::getDefault();

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
   ThicknessFluxDivOnCell ThickFluxDivOnC(Mesh);
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

   const auto Mesh = HorzMesh::getDefault();

   // Compute exact result
   Array2DReal ExactPotVortHAdv("ExactPotVortHAdv", Mesh->NEdgesOwned,
                                NVertLayers);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = (Setup.normRelVort(X, Y) + Setup.normPlanetVort(X, Y)) *
                        Setup.layerThick(X, Y) * Setup.vectorX(X, Y);
          VecField[1] = (Setup.normRelVort(X, Y) + Setup.normPlanetVort(X, Y)) *
                        Setup.layerThick(X, Y) * Setup.vectorY(X, Y);
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

   Array2DReal LayerThickEdge("LayerThickEdge", Mesh->NEdgesSize, NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.layerThick(X, Y); },
       LayerThickEdge, Geom, Mesh, OnEdge);

   Array2DReal NormVelEdge("NormVelEdge", Mesh->NEdgesSize, NVertLayers);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.vectorX(X, Y);
          VecField[1] = Setup.vectorY(X, Y);
       },
       NormVelEdge, EdgeComponent::Normal, Geom, Mesh);

   // Compute numerical result
   Array2DReal NumPotVortHAdv("NumPotVortHAdv", Mesh->NEdgesOwned, NVertLayers);

   PotentialVortHAdvOnEdge PotVortHAdvOnE(Mesh);
   parallelFor(
       {Mesh->NEdgesOwned, NVertLayers}, KOKKOS_LAMBDA(int IEdge, int KLayer) {
          PotVortHAdvOnE(NumPotVortHAdv, IEdge, KLayer, NormRelVortEdge,
                         NormPlanetVortEdge, LayerThickEdge, NormVelEdge);
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

   const auto Mesh = HorzMesh::getDefault();

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

   KEGradOnEdge KEGradOnE(Mesh);
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

   const auto Mesh = HorzMesh::getDefault();

   // Compute exact result
   Array2DReal ExactSSHGrad("ExactSSHGrad", Mesh->NEdgesOwned, NVertLayers);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = -Gravity * Setup.gradX(X, Y);
          VecField[1] = -Gravity * Setup.gradY(X, Y);
       },
       ExactSSHGrad, EdgeComponent::Normal, Geom, Mesh, ExchangeHalos::No);

   // Set input array
   Array2DReal SSHCell("SSHCell", Mesh->NCellsSize, NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.scalar(X, Y); }, SSHCell,
       Geom, Mesh, OnCell);

   // Compute numerical result
   Array2DReal NumSSHGrad("NumSSHGrad", Mesh->NEdgesOwned, NVertLayers);

   SSHGradOnEdge SSHGradOnE(Mesh);
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

   const auto Mesh = HorzMesh::getDefault();

   Config *OmegaConfig = Config::getOmegaConfig();
   Config TendConfig("Tendencies");
   Err1 += OmegaConfig->get(TendConfig);
   CHECK_ERROR_ABORT(Err1, "Tendencies: Tendencies group not found in Config");

   VelocityDiffusionOnEdge VelDiffOnE(Mesh);
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

   const auto Mesh = HorzMesh::getDefault();

   Config *OmegaConfig = Config::getOmegaConfig();
   Config TendConfig("Tendencies");
   Err1 += OmegaConfig->get(TendConfig);
   CHECK_ERROR_ABORT(Err1, "Tendencies: Tendencies group not found in Config");

   VelocityHyperDiffOnEdge VelHyperDiffOnE(Mesh);
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

   const auto Mesh = HorzMesh::getDefault();

   const Real SaltWaterDensity = 0.987654321;

   // Compute exact result
   Array2DReal ExactWindForcing("ExactWindForcing", Mesh->NEdgesOwned,
                                NVertLayers);

   // Note: this computes wind forcing at every layer
   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.windForcingX(X, Y, SaltWaterDensity);
          VecField[1] = Setup.windForcingY(X, Y, SaltWaterDensity);
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

   Array2DReal LayerThickEdge("LayerThickEdge", Mesh->NEdgesSize, NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.scalarB(X, Y); },
       LayerThickEdge, Geom, Mesh, OnEdge);

   // Compute numerical result
   Array2DReal NumWindForcing("NumWindForcing", Mesh->NEdgesOwned, NVertLayers);

   WindForcingOnEdge WindForcingOnE(Mesh);
   WindForcingOnE.SaltWaterDensity = SaltWaterDensity;

   parallelFor(
       {Mesh->NEdgesOwned, NVertLayers}, KOKKOS_LAMBDA(int IEdge, int KLayer) {
          WindForcingOnE(NumWindForcing, IEdge, KLayer, NormalStressEdge,
                         LayerThickEdge);
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

   Array2DReal LayerThickEdge("LayerThickEdge", Mesh->NEdgesSize, NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.scalarB(X, Y); },
       LayerThickEdge, Geom, Mesh, OnEdge);

   // Compute numerical result
   Array2DReal NumBottomDrag("NumBottomDrag", Mesh->NEdgesOwned, NVertLayers);

   BottomDragOnEdge BottomDragOnE(Mesh, VCoord);
   BottomDragOnE.Coeff = Coeff;

   parallelFor(
       {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge) {
          BottomDragOnE(NumBottomDrag, IEdge, NormalVelEdge, KECell,
                        LayerThickEdge);
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

   const auto Mesh = HorzMesh::getDefault();

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

   Array3DReal HTrOnEdge("HTrOnEdge", NTracers, Mesh->NEdgesSize, NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return -Setup.layerThick(X, Y); },
       HTrOnEdge, Geom, Mesh, OnEdge);

   // Compute numerical result
   Array3DReal NumTrFluxDiv("NumTrFluxDiv", NTracers, Mesh->NCellsOwned,
                            NVertLayers);
   TracerHorzAdvOnCell TrHorzAdvOnC(Mesh);
   parallelFor(
       {NTracers, Mesh->NCellsOwned, NVertLayers},
       KOKKOS_LAMBDA(int L, int ICell, int KLayer) {
          TrHorzAdvOnC(NumTrFluxDiv, L, ICell, KLayer, NormalVelocity,
                       HTrOnEdge);
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

   const auto Mesh = HorzMesh::getDefault();

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

   Array2DReal LayerThickEdge("LayerThickEdge", Mesh->NEdgesSize, NVertLayers);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.scalarB(X, Y); },
       LayerThickEdge, Geom, Mesh, OnEdge);

   // Compute numerical result
   Array3DReal NumTracerDiff("NumTracerDiff", NTracers, Mesh->NCellsOwned,
                             NVertLayers);
   TracerDiffOnCell TrDiffOnC(Mesh);
   TrDiffOnC.EddyDiff2 = 1._Real;

   parallelFor(
       {NTracers, Mesh->NCellsOwned, NVertLayers},
       KOKKOS_LAMBDA(int L, int ICell, int KLayer) {
          TrDiffOnC(NumTracerDiff, L, ICell, KLayer, TracerCell,
                    LayerThickEdge);
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

   const auto Mesh = HorzMesh::getDefault();

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
   TracerHyperDiffOnCell TrHypDiffOnC(Mesh);
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

void initTendTest(const std::string &MeshFile, int NVertLayers) {

   Error Err;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // Initialize logging
   initLogging(DefEnv);

   // Open config file
   Config("Omega");
   Config::readAll("omega.yml");

   IO::init(DefComm);

   Decomp::init(MeshFile);

   int HaloErr = Halo::init();
   if (HaloErr != 0) {
      ABORT_ERROR("TendencyTermsTest: error initializing default halo");
   }

   VertCoord::init1();

   // Reset NVertLayers to the test value
   auto *DefVertCoord        = VertCoord::getDefault();
   DefVertCoord->NVertLayers = NVertLayers;
   Dimension::destroy("NVertLayers");
   std::shared_ptr<Dimension> VertDim =
       Dimension::create("NVertLayers", NVertLayers);

   HorzMesh::init();

} // end initTendTest

void finalizeTendTest() {
   HorzMesh::clear();
   VertCoord::clear();
   Dimension::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();

} // end finalizeTendTest

int tendencyTermsTest(const std::string &MeshFile = DefaultMeshFile) {
   int Err         = 0;
   int NVertLayers = 16;

   initTendTest(MeshFile, NVertLayers);

   const auto &Mesh = HorzMesh::getDefault();
   int NTracers     = 3;

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
