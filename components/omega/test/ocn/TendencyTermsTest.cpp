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
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "Logging.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "mpi.h"

#include <cmath>

using namespace OMEGA;

struct TestSetupPlane {
   Real Pi = M_PI;

   Real Lx = 1;
   Real Ly = std::sqrt(3) / 2;

   ErrorMeasures ExpectedDivErrors     = {0.00124886886594453264,
                                          0.00124886886590977139};
   ErrorMeasures ExpectedPVErrors      = {0.00807347170900282914,
                                          0.00794755105765788429};
   ErrorMeasures ExpectedGradErrors    = {0.00125026071878537952,
                                          0.00134354611117262161};
   ErrorMeasures ExpectedLaplaceErrors = {0.00113090174765822192,
                                          0.00134324628763667899};
   ErrorMeasures ExpectedTrHAdvErrors  = {0.00205864372747571571,
                                          0.00172418025417940784};
   ErrorMeasures ExpectedTrDel2Errors  = {0.00334357193650093847,
                                          0.00290978146207349032};
   ErrorMeasures ExpectedTrDel4Errors  = {0.00508833446725232875,
                                          0.00523080740758275625};

   KOKKOS_FUNCTION Real vectorX(Real X, Real Y) const {
      return std::sin(2 * Pi * X / Lx) * std::cos(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real vectorY(Real X, Real Y) const {
      return std::cos(2 * Pi * X / Lx) * std::sin(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real divergence(Real X, Real Y) const {
      return 2 * Pi * (1. / Lx + 1. / Ly) * std::cos(2 * Pi * X / Lx) *
             std::cos(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real scalar(Real X, Real Y) const {
      return std::sin(2 * Pi * X / Lx) * std::sin(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real gradX(Real X, Real Y) const {
      return 2 * Pi / Lx * std::cos(2 * Pi * X / Lx) *
             std::sin(2 * Pi * Y / Ly);
   }
   KOKKOS_FUNCTION Real gradY(Real X, Real Y) const {
      return 2 * Pi / Ly * std::sin(2 * Pi * X / Lx) *
             std::cos(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real curl(Real X, Real Y) const {
      return 2 * Pi * (-1. / Lx + 1. / Ly) * std::sin(2 * Pi * X / Lx) *
             std::sin(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real laplaceVecX(Real X, Real Y) const {
      return -4 * Pi * Pi * (1. / Lx / Lx + 1. / Ly / Ly) *
             std::sin(2 * Pi * X / Lx) * std::cos(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real laplaceVecY(Real X, Real Y) const {
      return -4 * Pi * Pi * (1. / Lx / Lx + 1. / Ly / Ly) *
             std::cos(2 * Pi * X / Lx) * std::sin(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real layerThick(Real X, Real Y) const {
      return 2. + std::sin(2 * Pi * X / Lx) * std::cos(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real planetaryVort(Real X, Real Y) const {
      return std::cos(2 * Pi * X / Lx) * std::cos(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real normRelVort(Real X, Real Y) const {
      return curl(X, Y) / layerThick(X, Y);
   }

   KOKKOS_FUNCTION Real normPlanetVort(Real X, Real Y) const {
      return planetaryVort(X, Y) / layerThick(X, Y);
   }

   KOKKOS_FUNCTION Real tracerFluxDiv(Real X, Real Y) const {
      return (2 * Pi / (Lx * Ly)) *
             (std::cos(2 * Pi * X / Lx) *
              (2 * (Lx + Ly) * std::cos(2 * Pi * Y / Ly) +
               (Lx + 2 * Ly) * std::sin(2 * Pi * X / Lx) *
                   std::pow(std::cos(2 * Pi * Y / Ly), 2) -
               Lx * std::sin(2 * Pi * X / Lx) *
                   std::pow(std::sin(2 * Pi * Y / Ly), 2)));
   }

   KOKKOS_FUNCTION Real scalarA(Real X, Real Y) const {
      return std::cos(2 * Pi * X / Lx) * std::sin(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real scalarB(Real X, Real Y) const {
      return 2. + std::cos(2 * Pi * X / Lx) * std::cos(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real tracerDiff(Real X, Real Y) const {
      return -4 * Pi * Pi * std::sin(2 * Pi * Y / Ly) *
             (2 * (1 / Lx / Lx + 1 / Ly / Ly) * std::cos(2 * Pi * X / Lx) +
              (1 / Ly / Ly +
               (1 / Lx / Lx + 1 / Ly / Ly) * std::cos(4 * Pi * X / Lx)) *
                  std::cos(2 * Pi * Y / Ly));
   }

   KOKKOS_FUNCTION Real scalarC(Real X, Real Y) const {
      return std::pow(std::cos(2 * Pi * X / Lx), 2) -
             std::pow(std::sin(2 * Pi * Y / Ly), 2);
   }

   KOKKOS_FUNCTION Real tracerHyperDiff(Real X, Real Y) const {
      return -8 * Pi * Pi *
             (std::cos(4 * Pi * X / Lx) / Lx / Lx +
              std::cos(4 * Pi * Y / Ly) / Ly / Ly);
   }

}; // end TestSetupPlane

struct TestSetupSphere {
   // radius of spherical mesh
   // TODO: get this from the mesh
   Real Radius = 6371220;

   Real Pi = M_PI;

   ErrorMeasures ExpectedDivErrors     = {0.0136595773989796766,
                                          0.00367052484586384131};
   ErrorMeasures ExpectedPVErrors      = {0.0219217796608757037,
                                          0.0122537418367830303};
   ErrorMeasures ExpectedGradErrors    = {0.00187912292540623471,
                                          0.00149841802817334935};
   ErrorMeasures ExpectedLaplaceErrors = {0.281930203304510130,
                                          0.270530313560271740};
   ErrorMeasures ExpectedTrHAdvErrors  = {0.0132310202299444034,
                                          0.0038523368564029538};
   ErrorMeasures ExpectedTrDel2Errors  = {0.0486107109846934185,
                                          0.00507514214194892694};
   ErrorMeasures ExpectedTrDel4Errors  = {0.000819552466009620408,
                                          0.00064700084412871962};

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

int testThickFluxDiv(int NVertLevels, Real RTol) {

   int Err = 0;
   TestSetup Setup;

   const auto Mesh = HorzMesh::getDefault();

   // Compute exact result
   Array2DReal ExactThickFluxDiv("ExactThickFluxDiv", Mesh->NCellsOwned,
                                 NVertLevels);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return -Setup.divergence(X, Y); },
       ExactThickFluxDiv, Geom, Mesh, OnCell, NVertLevels, ExchangeHalos::No);

   // Set input array
   Array2DReal ThickFluxEdge("ThickFluxEdge", Mesh->NEdgesSize, NVertLevels);

   // TODO(mwarusz) temporary fix for this test
   Array2DReal OnesEdge("OnesEdge", Mesh->NEdgesSize, NVertLevels);
   deepCopy(OnesEdge, 1);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.vectorX(X, Y);
          VecField[1] = Setup.vectorY(X, Y);
       },
       ThickFluxEdge, EdgeComponent::Normal, Geom, Mesh, NVertLevels);

   // Compute numerical result
   Array2DReal NumThickFluxDiv("NumThickFluxDiv", Mesh->NCellsOwned,
                               NVertLevels);
   ThicknessFluxDivOnCell ThickFluxDivOnC(Mesh);
   parallelFor(
       {Mesh->NCellsOwned, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
          ThickFluxDivOnC(NumThickFluxDiv, ICell, KLevel, OnesEdge,
                          ThickFluxEdge);
       });

   // Compute errors
   ErrorMeasures TFDivErrors;
   Err += computeErrors(TFDivErrors, NumThickFluxDiv, ExactThickFluxDiv, Mesh,
                        OnCell, NVertLevels);

   // Check error values
   Err += checkErrors("TendencyTermsTest", "ThickFluxDiv", TFDivErrors,
                      Setup.ExpectedDivErrors, RTol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: ThickFluxDiv PASS");
   }

   return Err;
} // end testThickFluxDiv

int testPotVortHAdv(int NVertLevels, Real RTol) {

   int Err = 0;
   TestSetup Setup;

   const auto Mesh = HorzMesh::getDefault();

   // Compute exact result
   Array2DReal ExactPotVortHAdv("ExactPotVortHAdv", Mesh->NEdgesOwned,
                                NVertLevels);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = (Setup.normRelVort(X, Y) + Setup.normPlanetVort(X, Y)) *
                        Setup.layerThick(X, Y) * Setup.vectorX(X, Y);
          VecField[1] = (Setup.normRelVort(X, Y) + Setup.normPlanetVort(X, Y)) *
                        Setup.layerThick(X, Y) * Setup.vectorY(X, Y);
       },
       ExactPotVortHAdv, EdgeComponent::Tangential, Geom, Mesh, NVertLevels,
       ExchangeHalos::No);

   // Set input arrays
   Array2DReal NormRelVortEdge("NormRelVortEdge", Mesh->NEdgesSize,
                               NVertLevels);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.normRelVort(X, Y); },
       NormRelVortEdge, Geom, Mesh, OnEdge, NVertLevels);

   Array2DReal NormPlanetVortEdge("NormPlanetVortEdge", Mesh->NEdgesSize,
                                  NVertLevels);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.normPlanetVort(X, Y); },
       NormPlanetVortEdge, Geom, Mesh, OnEdge, NVertLevels);

   Array2DReal LayerThickEdge("LayerThickEdge", Mesh->NEdgesSize, NVertLevels);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.layerThick(X, Y); },
       LayerThickEdge, Geom, Mesh, OnEdge, NVertLevels);

   Array2DReal NormVelEdge("NormVelEdge", Mesh->NEdgesSize, NVertLevels);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.vectorX(X, Y);
          VecField[1] = Setup.vectorY(X, Y);
       },
       NormVelEdge, EdgeComponent::Normal, Geom, Mesh, NVertLevels);

   // Compute numerical result
   Array2DReal NumPotVortHAdv("NumPotVortHAdv", Mesh->NEdgesOwned, NVertLevels);

   PotentialVortHAdvOnEdge PotVortHAdvOnE(Mesh);
   parallelFor(
       {Mesh->NEdgesOwned, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
          PotVortHAdvOnE(NumPotVortHAdv, IEdge, KLevel, NormRelVortEdge,
                         NormPlanetVortEdge, LayerThickEdge, NormVelEdge);
       });

   // Compute errors
   ErrorMeasures PotVortHAdvErrors;
   Err += computeErrors(PotVortHAdvErrors, NumPotVortHAdv, ExactPotVortHAdv,
                        Mesh, OnEdge, NVertLevels);

   // Check error values
   Err += checkErrors("TendencyTermsTest", "PotVortHAdv", PotVortHAdvErrors,
                      Setup.ExpectedPVErrors, RTol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: PotVortHAdv PASS");
   }

   return Err;
} // end testPotVortHAdv

int testKEGrad(int NVertLevels, Real RTol) {

   int Err = 0;
   TestSetup Setup;

   const auto Mesh = HorzMesh::getDefault();

   // Compute exact result
   Array2DReal ExactKEGrad("ExactKEGrad", Mesh->NEdgesOwned, NVertLevels);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = -Setup.gradX(X, Y);
          VecField[1] = -Setup.gradY(X, Y);
       },
       ExactKEGrad, EdgeComponent::Normal, Geom, Mesh, NVertLevels,
       ExchangeHalos::No);

   // Set input array
   Array2DReal KECell("KECell", Mesh->NCellsSize, NVertLevels);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.scalar(X, Y); }, KECell,
       Geom, Mesh, OnCell, NVertLevels);

   // Compute numerical result
   Array2DReal NumKEGrad("NumKEGrad", Mesh->NEdgesOwned, NVertLevels);

   KEGradOnEdge KEGradOnE(Mesh);
   parallelFor(
       {Mesh->NEdgesOwned, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
          KEGradOnE(NumKEGrad, IEdge, KLevel, KECell);
       });

   // Compute errors
   ErrorMeasures KEGradErrors;
   Err += computeErrors(KEGradErrors, NumKEGrad, ExactKEGrad, Mesh, OnEdge,
                        NVertLevels);

   // Check error values
   Err += checkErrors("TendencyTermsTest", "KEGrad", KEGradErrors,
                      Setup.ExpectedGradErrors, RTol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: KEGrad PASS");
   }

   return Err;
} // end testKEGrad

int testSSHGrad(int NVertLevels, Real RTol) {

   int Err = 0;
   TestSetup Setup;

   const auto Mesh = HorzMesh::getDefault();

   // Compute exact result
   Array2DReal ExactSSHGrad("ExactSSHGrad", Mesh->NEdgesOwned, NVertLevels);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = -9.80665_Real * Setup.gradX(X, Y);
          VecField[1] = -9.80665_Real * Setup.gradY(X, Y);
       },
       ExactSSHGrad, EdgeComponent::Normal, Geom, Mesh, NVertLevels,
       ExchangeHalos::No);

   // Set input array
   Array2DReal SSHCell("SSHCell", Mesh->NCellsSize, NVertLevels);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.scalar(X, Y); }, SSHCell,
       Geom, Mesh, OnCell, NVertLevels);

   // Compute numerical result
   Array2DReal NumSSHGrad("NumSSHGrad", Mesh->NEdgesOwned, NVertLevels);

   SSHGradOnEdge SSHGradOnE(Mesh);
   parallelFor(
       {Mesh->NEdgesOwned, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
          SSHGradOnE(NumSSHGrad, IEdge, KLevel, SSHCell);
       });

   // Compute errors
   ErrorMeasures SSHGradErrors;
   Err += computeErrors(SSHGradErrors, NumSSHGrad, ExactSSHGrad, Mesh, OnEdge,
                        NVertLevels);

   // Check error values
   Err += checkErrors("TendencyTermsTest", "SSHGrad", SSHGradErrors,
                      Setup.ExpectedGradErrors, RTol);

   return Err;
} // end testSSHGrad

int testVelDiff(int NVertLevels, Real RTol) {

   int Err = 0;
   TestSetup Setup;

   const auto Mesh = HorzMesh::getDefault();

   Config *OmegaConfig = Config::getOmegaConfig();
   Config TendConfig("Tendencies");
   Err = OmegaConfig->get(TendConfig);
   if (Err != 0) {
      LOG_CRITICAL("Tendencies: Tendencies group not found in Config");
   }

   VelocityDiffusionOnEdge VelDiffOnE(Mesh);
   Err = TendConfig.get("ViscDel2", VelDiffOnE.ViscDel2);
   if (Err != 0) {
      LOG_CRITICAL("Tendencies: ViscDel2 not found in TendConfig");
   }
   const Real ViscDel2 = VelDiffOnE.ViscDel2;

   // Compute exact result
   Array2DReal ExactVelDiff("ExactVelDiff", Mesh->NEdgesOwned, NVertLevels);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = ViscDel2 * Setup.laplaceVecX(X, Y);
          VecField[1] = ViscDel2 * Setup.laplaceVecY(X, Y);
       },
       ExactVelDiff, EdgeComponent::Normal, Geom, Mesh, NVertLevels,
       ExchangeHalos::No);

   // Set input arrays
   Array2DReal DivCell("DivCell", Mesh->NCellsSize, NVertLevels);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.divergence(X, Y); },
       DivCell, Geom, Mesh, OnCell, NVertLevels);

   Array2DReal RVortVertex("RVortVertex", Mesh->NVerticesSize, NVertLevels);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.curl(X, Y); }, RVortVertex,
       Geom, Mesh, OnVertex, NVertLevels);

   // Compute numerical result
   Array2DReal NumVelDiff("NumVelDiff", Mesh->NEdgesOwned, NVertLevels);

   parallelFor(
       {Mesh->NEdgesOwned, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
          VelDiffOnE(NumVelDiff, IEdge, KLevel, DivCell, RVortVertex);
       });

   // Compute errors
   ErrorMeasures VelDiffErrors;
   Err += computeErrors(VelDiffErrors, NumVelDiff, ExactVelDiff, Mesh, OnEdge,
                        NVertLevels);

   // Check error values
   Err += checkErrors("TendencyTermsTest", "VelDiff", VelDiffErrors,
                      Setup.ExpectedLaplaceErrors, RTol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: VelDiff PASS");
   }

   return Err;
} // end testVelDiff

int testVelHyperDiff(int NVertLevels, Real RTol) {

   int Err = 0;
   TestSetup Setup;

   const auto Mesh = HorzMesh::getDefault();

   Config *OmegaConfig = Config::getOmegaConfig();
   Config TendConfig("Tendencies");
   Err = OmegaConfig->get(TendConfig);
   if (Err != 0) {
      LOG_CRITICAL("Tendencies: Tendencies group not found in Config");
   }

   VelocityHyperDiffOnEdge VelHyperDiffOnE(Mesh);
   Err = TendConfig.get("ViscDel4", VelHyperDiffOnE.ViscDel4);
   if (Err != 0) {
      LOG_CRITICAL("Tendencies: ViscDel4 not found in TendConfig");
   }
   const Real ViscDel4 = VelHyperDiffOnE.ViscDel4;

   // Compute exact result
   Array2DReal ExactVelHyperDiff("ExactVelHyperDiff", Mesh->NEdgesOwned,
                                 NVertLevels);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = -ViscDel4 * Setup.laplaceVecX(X, Y);
          VecField[1] = -ViscDel4 * Setup.laplaceVecY(X, Y);
       },
       ExactVelHyperDiff, EdgeComponent::Normal, Geom, Mesh, NVertLevels,
       ExchangeHalos::No);

   // Set input arrays
   Array2DReal DivCell("DivCell", Mesh->NCellsSize, NVertLevels);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.divergence(X, Y); },
       DivCell, Geom, Mesh, OnCell, NVertLevels);

   Array2DReal RVortVertex("RVortVertex", Mesh->NVerticesSize, NVertLevels);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.curl(X, Y); }, RVortVertex,
       Geom, Mesh, OnVertex, NVertLevels);

   // Compute numerical result
   Array2DReal NumVelHyperDiff("NumVelHyperDiff", Mesh->NEdgesOwned,
                               NVertLevels);

   parallelFor(
       {Mesh->NEdgesOwned, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
          VelHyperDiffOnE(NumVelHyperDiff, IEdge, KLevel, DivCell, RVortVertex);
       });

   // Compute errors
   ErrorMeasures VelHyperDiffErrors;
   Err += computeErrors(VelHyperDiffErrors, NumVelHyperDiff, ExactVelHyperDiff,
                        Mesh, OnEdge, NVertLevels);

   // Check error values
   Err += checkErrors("TendencyTermsTest", "VelHyperDiff", VelHyperDiffErrors,
                      Setup.ExpectedLaplaceErrors, RTol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: VelHyperDiff PASS");
   }

   return Err;
} // end testVelHyperDiff

int testTracerHorzAdvOnCell(int NVertLevels, int NTracers, Real RTol) {

   I4 Err = 0;
   TestSetup Setup;

   const auto Mesh = HorzMesh::getDefault();

   // Compute exact result
   Array3DReal ExactTrFluxDiv("ExactTrFluxDiv", NTracers, Mesh->NCellsOwned,
                              NVertLevels);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.tracerFluxDiv(X, Y); },
       ExactTrFluxDiv, Geom, Mesh, OnCell, NVertLevels, NTracers,
       ExchangeHalos::No);

   // Set input arrays
   Array2DReal NormalVelocity("NormalVelocity", Mesh->NEdgesSize, NVertLevels);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.vectorX(X, Y);
          VecField[1] = Setup.vectorY(X, Y);
       },
       NormalVelocity, EdgeComponent::Normal, Geom, Mesh, NVertLevels);

   Array3DReal HTrOnEdge("HTrOnEdge", NTracers, Mesh->NEdgesSize, NVertLevels);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return -Setup.layerThick(X, Y); },
       HTrOnEdge, Geom, Mesh, OnEdge, NVertLevels, NTracers);

   // Compute numerical result
   Array3DReal NumTrFluxDiv("NumTrFluxDiv", NTracers, Mesh->NCellsOwned,
                            NVertLevels);
   TracerHorzAdvOnCell TrHorzAdvOnC(Mesh);
   parallelFor(
       {NTracers, Mesh->NCellsOwned, NVertLevels},
       KOKKOS_LAMBDA(int L, int ICell, int KLevel) {
          TrHorzAdvOnC(NumTrFluxDiv, L, ICell, KLevel, NormalVelocity,
                       HTrOnEdge);
       });

   ErrorMeasures TrHAdvErrors;
   Err += computeErrors(TrHAdvErrors, NumTrFluxDiv, ExactTrFluxDiv, Mesh,
                        OnCell, NVertLevels, NTracers);

   Err += checkErrors("TendencyTermsTest", "TracerHorzAdv", TrHAdvErrors,
                      Setup.ExpectedTrHAdvErrors, RTol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: TracerHorzAdv PASS");
   }

   return Err;
} // end testTracerHorzAdvOnCell

int testTracerDiffOnCell(int NVertLevels, int NTracers, Real RTol) {

   I4 Err = 0;
   TestSetup Setup;

   const auto Mesh = HorzMesh::getDefault();

   // Compute exact result
   Array3DReal ExactTracerDiff("ExactTracerDiff", NTracers, Mesh->NCellsOwned,
                               NVertLevels);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.tracerDiff(X, Y); },
       ExactTracerDiff, Geom, Mesh, OnCell, NVertLevels, NTracers,
       ExchangeHalos::No);

   // Set input arrays
   Array3DReal TracerCell("TracerCell", NTracers, Mesh->NCellsSize,
                          NVertLevels);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.scalarA(X, Y); },
       TracerCell, Geom, Mesh, OnCell, NVertLevels, NTracers);

   Array2DReal LayerThickEdge("LayerThickEdge", Mesh->NEdgesSize, NVertLevels);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.scalarB(X, Y); },
       LayerThickEdge, Geom, Mesh, OnEdge, NVertLevels);

   // Compute numerical result
   Array3DReal NumTracerDiff("NumTracerDiff", NTracers, Mesh->NCellsOwned,
                             NVertLevels);
   TracerDiffOnCell TrDiffOnC(Mesh);
   TrDiffOnC.EddyDiff2 = 1._Real;

   parallelFor(
       {NTracers, Mesh->NCellsOwned, NVertLevels},
       KOKKOS_LAMBDA(int L, int ICell, int KLevel) {
          TrDiffOnC(NumTracerDiff, L, ICell, KLevel, TracerCell,
                    LayerThickEdge);
       });

   ErrorMeasures TrDiffErrors;
   Err += computeErrors(TrDiffErrors, NumTracerDiff, ExactTracerDiff, Mesh,
                        OnCell, NVertLevels, NTracers);

   Err += checkErrors("TendencyTermsTest", "TracerDiff", TrDiffErrors,
                      Setup.ExpectedTrDel2Errors, RTol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: TracerDiff PASS");
   }

   return Err;
} // end testTracerDiffOnCell

int testTracerHyperDiffOnCell(int NVertLevels, int NTracers, Real RTol) {

   I4 Err = 0;
   TestSetup Setup;

   const auto Mesh = HorzMesh::getDefault();

   // Compute exact result
   Array3DReal ExactTracerHyperDiff("ExactTracerHyperDiff", NTracers,
                                    Mesh->NCellsOwned, NVertLevels);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return -Setup.tracerHyperDiff(X, Y); },
       ExactTracerHyperDiff, Geom, Mesh, OnCell, NVertLevels, NTracers,
       ExchangeHalos::No);

   // Set input arrays
   Array3DReal TrDel2Cell("TracerCell", NTracers, Mesh->NCellsSize,
                          NVertLevels);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.scalarC(X, Y); },
       TrDel2Cell, Geom, Mesh, OnCell, NVertLevels, NTracers);

   // Compute numerical result
   Array3DReal NumTracerHyperDiff("NumTracerHyperDiff", NTracers,
                                  Mesh->NCellsOwned, NVertLevels);
   TracerHyperDiffOnCell TrHypDiffOnC(Mesh);
   TrHypDiffOnC.EddyDiff4 = 1._Real;
   parallelFor(
       {NTracers, Mesh->NCellsOwned, NVertLevels},
       KOKKOS_LAMBDA(int L, int ICell, int KLevel) {
          TrHypDiffOnC(NumTracerHyperDiff, L, ICell, KLevel, TrDel2Cell);
       });

   ErrorMeasures TrHyperDiffErrors;
   Err +=
       computeErrors(TrHyperDiffErrors, NumTracerHyperDiff,
                     ExactTracerHyperDiff, Mesh, OnCell, NVertLevels, NTracers);

   Err += checkErrors("TendencyTermsTest", "TracerHyperDiff", TrHyperDiffErrors,
                      Setup.ExpectedTrDel4Errors, RTol);

   if (Err == 0) {
      LOG_INFO("TendencyTermsTest: TracerHyperDiff PASS");
   }

   return Err;
} // end testTracerHyperDiffOnCell

int initTendTest(const std::string &mesh) {

   I4 Err = 0;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // Initialize logging
   initLogging(DefEnv);

   // Open config file
   Config("Omega");
   Err = Config::readAll("omega.yml");
   if (Err != 0) {
      LOG_CRITICAL("TendencyTermsTest: Error reading config file");
      return Err;
   }

   I4 IOErr = IO::init(DefComm);
   if (IOErr != 0) {
      Err++;
      LOG_ERROR("TendencyTermsTest: error initializing parallel IO");
   }

   int DecompErr = Decomp::init(mesh);
   if (DecompErr != 0) {
      Err++;
      LOG_ERROR("TendencyTermsTest: error initializing default decomposition");
   }

   int HaloErr = Halo::init();
   if (HaloErr != 0) {
      Err++;
      LOG_ERROR("TendencyTermsTest: error initializing default halo");
   }

   int MeshErr = HorzMesh::init();
   if (MeshErr != 0) {
      Err++;
      LOG_ERROR("TendencyTermsTest: error initializing default mesh");
   }

   return Err;
} // end initTendTest

void finalizeTendTest() {

   HorzMesh::clear();
   Dimension::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();

} // end finalizeTendTest

int tendencyTermsTest(const std::string &mesh = DefaultMeshFile) {

   int Err = initTendTest(mesh);
   if (Err != 0) {
      LOG_CRITICAL("TendencyTermsTest: Error initializing");
   }

   const auto &Mesh = HorzMesh::getDefault();
   int NVertLevels  = 16;
   int NTracers     = 3;

   const Real RTol = sizeof(Real) == 4 ? 2e-2 : 1e-10;

   Err += testThickFluxDiv(NVertLevels, RTol);

   Err += testPotVortHAdv(NVertLevels, RTol);

   Err += testKEGrad(NVertLevels, RTol);

   Err += testSSHGrad(NVertLevels, RTol);

   Err += testVelDiff(NVertLevels, RTol);

   Err += testVelHyperDiff(NVertLevels, RTol);

   Err += testTracerHorzAdvOnCell(NVertLevels, NTracers, RTol);

   Err += testTracerDiffOnCell(NVertLevels, NTracers, RTol);

   Err += testTracerHyperDiffOnCell(NVertLevels, NTracers, RTol);

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

   RetErr = tendencyTermsTest();

   Kokkos::finalize();
   MPI_Finalize();

   return RetErr;

} // end of main
//===-----------------------------------------------------------------------===/
