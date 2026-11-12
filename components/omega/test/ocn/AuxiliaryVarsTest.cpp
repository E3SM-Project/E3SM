#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Field.h"
#include "GlobalConstants.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "VertCoord.h"
#include "auxiliaryVars/KineticAuxVars.h"
#include "auxiliaryVars/PseudoThicknessAuxVars.h"
#include "auxiliaryVars/TracerAuxVars.h"
#include "auxiliaryVars/VelocityDel2AuxVars.h"
#include "auxiliaryVars/VorticityAuxVars.h"
#include "mpi.h"

#include <cmath>
#include <iomanip>

using namespace OMEGA;

struct TestSetupPlane {

   Real Lx = 1;
   Real Ly = SqrtThree / 2;

   ErrorMeasures ExpectedKineticEnergyErrors = {0.00994439065100057897,
                                                0.00703403756741667954};
   ErrorMeasures ExpectedVelocityDivErrors   = {0.00124886886594453264,
                                                0.00124886886590973452};

   ErrorMeasures ExpectedFluxThickErrors = {0.0218166134247192549,
                                            0.0171404379252105554};
   ErrorMeasures ExpectedMeanThickErrors = {0.000890795148016506602,
                                            0.000741722075349612398};

   ErrorMeasures ExpectedRelVortVertexErrors        = {0.161365663569687623,
                                                       0.161348016897141511};
   ErrorMeasures ExpectedNormRelVortVertexErrors    = {0.185771689108325755,
                                                       0.170080698606596442};
   ErrorMeasures ExpectedNormPlanetVortVertexErrors = {0.000831626192159380336,
                                                       0.000562164971653627546};

   ErrorMeasures ExpectedNormRelVortEdgeErrors    = {0.0119295506805566498,
                                                     0.00779991259802507997};
   ErrorMeasures ExpectedNormPlanetVortEdgeErrors = {0.00223924332422219697,
                                                     0.0015382243254998785};

   ErrorMeasures ExpectedDel2Errors        = {0.00113090174765806731,
                                              0.00134324628763670241};
   ErrorMeasures ExpectedDel2DivErrors     = {0.002495925826729385,
                                              0.00249592582669975289};
   ErrorMeasures ExpectedDel2RelVortErrors = {0.0104455692965114266,
                                              0.0104135556263709097};

   ErrorMeasures ExpectedHTracerErrors    = {0.017402432114157595,
                                             0.00813360234680596434};
   ErrorMeasures ExpectedDel2TracerErrors = {0.0033346711042859123,
                                             0.0029202923731303323};

   KOKKOS_FUNCTION Real pseudoThickness(Real X, Real Y) const {
      return 2 + std::cos(TwoPi * X / Lx) * std::cos(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real velocityX(Real X, Real Y) const {
      return std::sin(TwoPi * X / Lx) * std::cos(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real velocityY(Real X, Real Y) const {
      return std::cos(TwoPi * X / Lx) * std::sin(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real divergence(Real X, Real Y) const {
      return TwoPi * (1. / Lx + 1. / Ly) * std::cos(TwoPi * X / Lx) *
             std::cos(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real relativeVorticity(Real X, Real Y) const {
      return TwoPi * (-1. / Lx + 1. / Ly) * std::sin(TwoPi * X / Lx) *
             std::sin(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real velocityDel2X(Real X, Real Y) const {
      return -TwoPi * TwoPi * (1 / (Lx * Lx) + 1 / (Ly * Ly)) * velocityX(X, Y);
   }

   KOKKOS_FUNCTION Real velocityDel2Y(Real X, Real Y) const {
      return -TwoPi * TwoPi * (1 / (Lx * Lx) + 1 / (Ly * Ly)) * velocityY(X, Y);
   }

   KOKKOS_FUNCTION Real velocityDel2Div(Real X, Real Y) const {
      return -TwoPi * TwoPi * (1 / (Lx * Lx) + 1 / (Ly * Ly)) *
             divergence(X, Y);
   }

   KOKKOS_FUNCTION Real velocityDel2Curl(Real X, Real Y) const {
      return -TwoPi * TwoPi * (1 / (Lx * Lx) + 1 / (Ly * Ly)) *
             relativeVorticity(X, Y);
   }

   KOKKOS_FUNCTION Real planetaryVorticity(Real X, Real Y) const {
      return std::sin(TwoPi * X / Lx) * std::sin(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real normalizedRelativeVorticity(Real X, Real Y) const {
      return relativeVorticity(X, Y) / pseudoThickness(X, Y);
   }

   KOKKOS_FUNCTION Real normalizedPlanetaryVorticity(Real X, Real Y) const {
      return planetaryVorticity(X, Y) / pseudoThickness(X, Y);
   }

   KOKKOS_FUNCTION Real kineticEnergy(Real X, Real Y) const {
      return (velocityX(X, Y) * velocityX(X, Y) +
              velocityY(X, Y) * velocityY(X, Y)) /
             2;
   }

   KOKKOS_FUNCTION Real tracer(Real X, Real Y) const {
      return 2 - std::cos(TwoPi * X / Lx) * std::cos(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real thickTracer(Real X, Real Y) const {
      return 4 - std::pow(std::cos(TwoPi * X / Lx), 2) *
                     std::pow(std::cos(TwoPi * Y / Ly), 2);
   }

   KOKKOS_FUNCTION Real del2Tracer(Real X, Real Y) const {
      return TwoPi * Pi *
             (4 * (1 / Lx / Lx + 1 / Ly / Ly) * std::cos(TwoPi * X / Lx) *
                  std::cos(TwoPi * Y / Ly) +
              std::pow(std::cos(TwoPi * X / Lx), 2) *
                  (1 / Lx / Lx +
                   (2 / Ly / Ly + 1 / Lx / Lx) * std::cos(2 * TwoPi * Y / Ly)) -
              (2 / Lx / Lx) * std::pow(std::sin(TwoPi * X / Lx), 2) *
                  std::pow(std::cos(TwoPi * Y / Ly), 2));
   }
};

struct TestSetupSphere {

   // radius of spherical mesh
   Real Radius = REarth;

   ErrorMeasures ExpectedKineticEnergyErrors = {0.014358161778181201,
                                                0.0068116910597281113};
   ErrorMeasures ExpectedVelocityDivErrors   = {0.013659556526123424,
                                                0.0036698023569595927};

   ErrorMeasures ExpectedFluxThickErrors = {0.015981055597025,
                                            0.010364584620791126};
   ErrorMeasures ExpectedMeanThickErrors = {0.00079985423942030413,
                                            0.00040654771018863028};

   ErrorMeasures ExpectedRelVortVertexErrors        = {0.027146182264485331,
                                                       0.025202560969432641};
   ErrorMeasures ExpectedNormRelVortVertexErrors    = {0.034882855311945148,
                                                       0.025950288113940474};
   ErrorMeasures ExpectedNormPlanetVortVertexErrors = {0.0045134270290698875,
                                                       0.0010179210697515548};

   ErrorMeasures ExpectedNormRelVortEdgeErrors    = {0.012528744550304859,
                                                     0.0030749833896684914};
   ErrorMeasures ExpectedNormPlanetVortEdgeErrors = {0.0049522329721762544,
                                                     0.0008554077431009539};

   ErrorMeasures ExpectedDel2Errors        = {0.0036017312925143057,
                                              0.003134053694265801};
   ErrorMeasures ExpectedDel2DivErrors     = {0.017767406015402661,
                                              0.007518774954098988};
   ErrorMeasures ExpectedDel2RelVortErrors = {0.091446507985509698,
                                              0.024667998533570906};

   ErrorMeasures ExpectedHTracerErrors    = {0.01603249913425972,
                                             0.00546762028673672059};
   ErrorMeasures ExpectedDel2TracerErrors = {0.0081209686316016313,
                                             0.0049181207204490136};

   KOKKOS_FUNCTION Real pseudoThickness(Real Lon, Real Lat) const {
      return (2 + std::cos(Lon) * std::pow(std::cos(Lat), 4));
   }

   KOKKOS_FUNCTION Real velocityX(Real Lon, Real Lat) const {
      return -std::pow(std::sin(Lon), 2) * std::pow(std::cos(Lat), 3);
   }

   KOKKOS_FUNCTION Real velocityY(Real Lon, Real Lat) const {
      return -4 * std::sin(Lon) * std::cos(Lon) * std::pow(std::cos(Lat), 3) *
             std::sin(Lat);
   }

   KOKKOS_FUNCTION Real relativeVorticity(Real Lon, Real Lat) const {
      return -4 * std::pow(std::cos(Lon), 2) * std::pow(std::cos(Lat), 2) *
             std::sin(Lat) / Radius;
   }

   KOKKOS_FUNCTION Real divergence(Real Lon, Real Lat) const {
      return std::sin(Lon) * std::cos(Lon) * std::pow(std::cos(Lat), 2) *
             (20 * std::pow(std::sin(Lat), 2) - 6) / Radius;
   }

   KOKKOS_FUNCTION Real velocityDel2X(Real Lon, Real Lat) const {
      return 1 / (Radius * Radius) *
                 (std::pow(std::cos(Lon), 2) - std::pow(std::sin(Lon), 2)) *
                 std::cos(Lat) * (20 * std::pow(std::sin(Lat), 2) - 6) +
             4 / (Radius * Radius) * std::pow(cos(Lon), 2) *
                 (std::pow(cos(Lat), 3) -
                  2 * std::cos(Lat) * std::pow(sin(Lat), 2));
   }

   KOKKOS_FUNCTION Real velocityDel2Y(Real Lon, Real Lat) const {
      return 1 / (Radius * Radius) * std::sin(Lon) * std::cos(Lon) *
                 std::sin(Lat) * std::cos(Lat) *
                 (80 * std::pow(std::cos(Lat), 2) - 28) +
             8 / (Radius * Radius) * std::sin(Lon) * std::cos(Lon) *
                 std::sin(Lat) * std::cos(Lat);
   }

   KOKKOS_FUNCTION Real velocityDel2Div(Real Lon, Real Lat) const {
      return 1 / (Radius * Radius * Radius) *
             (-2 * std::sin(Lon) * std::cos(Lon) *
                  (28 * std::pow(sin(Lat), 2) - 8) +
              std::sin(Lon) * std::cos(Lon) *
                  ((std::pow(cos(Lat), 2) - 2 * std::pow(sin(Lat), 2)) *
                       (80 * std::pow(cos(Lat), 2) - 20) -
                   160 * std::pow(sin(Lat) * cos(Lat), 2)));
   }

   KOKKOS_FUNCTION Real velocityDel2Curl(Real Lon, Real Lat) const {
      return 1 / (Radius * Radius * Radius) *
             (-std::sin(Lat) * (std::pow(std::cos(Lat), 2) *
                                    (56 * std::pow(cos(Lon), 2) - 40) -
                                2 * (std::pow(cos(Lon), 2) *
                                         (28 * std::pow(sin(Lat), 2) - 8) -
                                     20 * std::pow(sin(Lat), 2) + 6)) +
              std::sin(Lat) * (80 * std::pow(cos(Lat), 2) - 20) *
                  (std::pow(cos(Lon), 2) - std::pow(sin(Lon), 2)));
   }

   KOKKOS_FUNCTION Real planetaryVorticity(Real Lon, Real Lat) const {
      return std::sin(Lat);
   }

   KOKKOS_FUNCTION Real normalizedRelativeVorticity(Real Lon, Real Lat) const {
      return relativeVorticity(Lon, Lat) / pseudoThickness(Lon, Lat);
   }

   KOKKOS_FUNCTION Real normalizedPlanetaryVorticity(Real Lon, Real Lat) const {
      return planetaryVorticity(Lon, Lat) / pseudoThickness(Lon, Lat);
   }

   KOKKOS_FUNCTION Real kineticEnergy(Real Lon, Real Lat) const {
      return (velocityX(Lon, Lat) * velocityX(Lon, Lat) +
              velocityY(Lon, Lat) * velocityY(Lon, Lat)) /
             2;
   }

   KOKKOS_FUNCTION Real tracer(Real Lon, Real Lat) const {
      return (2 - std::cos(Lon) * std::pow(std::cos(Lat), 4));
   }

   KOKKOS_FUNCTION Real thickTracer(Real Lon, Real Lat) const {
      return (4 - std::pow(std::cos(Lon), 2) * std::pow(std::cos(Lat), 8));
   }

   KOKKOS_FUNCTION Real del2Tracer(Real Lon, Real Lat) const {
      return 1 / (Radius * Radius) *
             (10 * std::cos(Lon) * std::pow(std::cos(Lat), 2) *
                  (-1 + 2 * std::cos(2 * Lat)) +
              std::pow(std::cos(Lon), 2) * std::pow(std::cos(Lat), 6) *
                  (-13 + 18 * std::cos(2 * Lat)) -
              std::pow(std::cos(Lat), 6) * std::pow(std::sin(Lon), 2));
   }
};

#ifdef AUXVARS_TEST_PLANE
constexpr Geometry Geom          = Geometry::Planar;
constexpr char DefaultMeshFile[] = "OmegaPlanarMesh.nc";
using TestSetup                  = TestSetupPlane;
#else
constexpr Geometry Geom          = Geometry::Spherical;
constexpr char DefaultMeshFile[] = "OmegaSphereMesh.nc";
using TestSetup                  = TestSetupSphere;
#endif

constexpr int NVertLayers = 16;
constexpr int NTracers    = 3;

int initState(const Array2DReal &PseudoThickCell,
              const Array2DReal &NormalVelEdge, HorzMesh *Mesh) {
   int Err = 0;

   TestSetup Setup;

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.pseudoThickness(X, Y); },
       PseudoThickCell, Geom, Mesh, OnCell);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.velocityX(X, Y);
          VecField[1] = Setup.velocityY(X, Y);
       },
       NormalVelEdge, EdgeComponent::Normal, Geom, Mesh);

   // need to override FVertex with prescribed values
   const auto &FVertex = Mesh->FVertex;

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.planetaryVorticity(X, Y); },
       FVertex, Geom, Mesh, OnVertex);

   return Err;
}

int testKineticAuxVars(const Array2DReal &PseudoThicknessCell,
                       const Array2DReal &NormalVelocityEdge, Real RTol) {
   int Err = 0;
   TestSetup Setup;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   // Compute exact result

   Array2DReal ExactKineticEnergyCell("ExactKineticEnergyCell",
                                      Mesh->NCellsOwned, NVertLayers);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.kineticEnergy(X, Y); },
       ExactKineticEnergyCell, Geom, Mesh, OnCell, ExchangeHalos::No);

   Array2DReal ExactVelocityDivCell("ExactVelocityDivCell", Mesh->NCellsOwned,
                                    NVertLayers);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.divergence(X, Y); },
       ExactVelocityDivCell, Geom, Mesh, OnCell, ExchangeHalos::No);

   // Compute numerical result

   KineticAuxVars KineticAux("", Mesh, VCoord);

   parallelForOuter(
       LaunchConfig({Mesh->NCellsOwned}, TeamScratch<Real>(2 * NVertLayers)),
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          KineticAux.computeVarsOnCell(Team, ICell, NormalVelocityEdge);
       });
   const auto &NumKineticEnergyCell = KineticAux.KineticEnergyCell;
   const auto &NumVelocityDivCell   = KineticAux.VelocityDivCell;

   // Compute error measures and check error values

   ErrorMeasures KineticEnergyErrors;
   Err += computeErrors(KineticEnergyErrors, NumKineticEnergyCell,
                        ExactKineticEnergyCell, Mesh, OnCell);
   Err += checkErrors("AuxVarsTest", "KineticEnergy", KineticEnergyErrors,
                      Setup.ExpectedKineticEnergyErrors, RTol);

   ErrorMeasures VelocityDivErrors;
   Err += computeErrors(VelocityDivErrors, NumVelocityDivCell,
                        ExactVelocityDivCell, Mesh, OnCell);
   Err += checkErrors("AuxVarsTest", "VelocityDiv", VelocityDivErrors,
                      Setup.ExpectedVelocityDivErrors, RTol);

   if (Err == 0) {
      LOG_INFO("AuxVarsTest: KineticAuxVars PASS");
   }

   return Err;
}

int testPseudoThicknessAuxVars(const Array2DReal &PseudoThickCell,
                               const Array2DReal &NormalVelEdge, Real RTol) {
   int Err = 0;
   TestSetup Setup;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   // Compute exact result

   Array2DReal ExactThickEdge("ExactThickEdge", Mesh->NEdgesOwned, NVertLayers);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.pseudoThickness(X, Y); },
       ExactThickEdge, Geom, Mesh, OnEdge, ExchangeHalos::No);

   // Compute numerical result

   PseudoThicknessAuxVars PseudoThicknessAux("", Mesh, VCoord);
   PseudoThicknessAux.FluxThickEdgeChoice = FluxThickEdgeOption::Upwind;
   parallelForOuter(
       {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          PseudoThicknessAux.computeVarsOnEdge(Team, IEdge, PseudoThickCell,
                                               NormalVelEdge);
       });

   const auto &NumFluxPseudoThickEdge = PseudoThicknessAux.FluxPseudoThickEdge;
   const auto &NumMeanPseudoThickEdge = PseudoThicknessAux.MeanPseudoThickEdge;

   // Compute error measures and check error values

   ErrorMeasures FluxThickErrors;
   Err += computeErrors(FluxThickErrors, NumFluxPseudoThickEdge, ExactThickEdge,
                        Mesh, OnEdge);
   Err += checkErrors("AuxVarsTest", "FluxThick", FluxThickErrors,
                      Setup.ExpectedFluxThickErrors, RTol);

   ErrorMeasures MeanThickErrors;
   Err += computeErrors(MeanThickErrors, NumMeanPseudoThickEdge, ExactThickEdge,
                        Mesh, OnEdge);
   Err += checkErrors("AuxVarsTest", "MeanThick", MeanThickErrors,
                      Setup.ExpectedMeanThickErrors, RTol);

   if (Err == 0) {
      LOG_INFO("AuxVarsTest: PseudoThicknessAuxVars PASS");
   }

   return Err;
}

int testVorticityAuxVars(const Array2DReal &PseudoThickCell,
                         const Array2DReal &NormalVelEdge, Real RTol) {
   TestSetup Setup;
   int Err = 0;

   const auto Decomp = Decomp::getDefault();
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();
   VorticityAuxVars VorticityAux("", Mesh, VCoord);

   // Compute exact results for vertex variables

   Array2DReal ExactRelVortVertex("ExactRelVortVertex", Mesh->NVerticesOwned,
                                  NVertLayers);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.relativeVorticity(X, Y); },
       ExactRelVortVertex, Geom, Mesh, OnVertex, ExchangeHalos::No);

   Array2DReal ExactNormRelVortVertex("ExactNormRelVortVertex",
                                      Mesh->NVerticesOwned, NVertLayers);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) {
          return Setup.normalizedRelativeVorticity(X, Y);
       },
       ExactNormRelVortVertex, Geom, Mesh, OnVertex, ExchangeHalos::No);

   Array2DReal ExactNormPlanetVortVertex("ExactNormPlanetVortVertex",
                                         Mesh->NVerticesOwned, NVertLayers);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) {
          return Setup.normalizedPlanetaryVorticity(X, Y);
       },
       ExactNormPlanetVortVertex, Geom, Mesh, OnVertex, ExchangeHalos::No);

   // Compute numerical results for vertex variables

   parallelForOuter(
       LaunchConfig({Decomp->NVerticesHaloH(0)},
                    TeamScratch<Real>(2 * VCoord->NVertLayers)),
       KOKKOS_LAMBDA(int IVertex, const TeamMember &Team) {
          VorticityAux.computeVarsOnVertex(Team, IVertex, PseudoThickCell,
                                           NormalVelEdge);
       });

   const auto &NumRelVortVertex        = VorticityAux.RelVortVertex;
   const auto &NumNormRelVortVertex    = VorticityAux.NormRelVortVertex;
   const auto &NumNormPlanetVortVertex = VorticityAux.NormPlanetVortVertex;

   // Compute error measures and check errors for vertex variables

   ErrorMeasures RelVortVertexErrors;
   Err += computeErrors(RelVortVertexErrors, NumRelVortVertex,
                        ExactRelVortVertex, Mesh, OnVertex);
   Err += checkErrors("AuxVarsTest", "RelVortVertex", RelVortVertexErrors,
                      Setup.ExpectedRelVortVertexErrors, RTol);

   ErrorMeasures NormRelVortVertexErrors;
   Err += computeErrors(NormRelVortVertexErrors, NumNormRelVortVertex,
                        ExactNormRelVortVertex, Mesh, OnVertex);
   Err +=
       checkErrors("AuxVarsTest", "NormRelVortVertex", NormRelVortVertexErrors,
                   Setup.ExpectedNormRelVortVertexErrors, RTol);

   ErrorMeasures NormPlanetVortVertexErrors;
   Err += computeErrors(NormPlanetVortVertexErrors, NumNormPlanetVortVertex,
                        ExactNormPlanetVortVertex, Mesh, OnVertex);
   Err += checkErrors("AuxVarsTest", "NormPlanetVortVertex",
                      NormPlanetVortVertexErrors,
                      Setup.ExpectedNormPlanetVortVertexErrors, RTol);

   // Compute exact results for edge variables

   Array2DReal ExactNormRelVortEdge("ExactNormRelVortEdge", Mesh->NEdgesOwned,
                                    NVertLayers);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) {
          return Setup.normalizedRelativeVorticity(X, Y);
       },
       ExactNormRelVortEdge, Geom, Mesh, OnEdge, ExchangeHalos::No);

   Array2DReal ExactNormPlanetVortEdge("ExactNormPlanetVortEdge",
                                       Mesh->NEdgesOwned, NVertLayers);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) {
          return Setup.normalizedPlanetaryVorticity(X, Y);
       },
       ExactNormPlanetVortEdge, Geom, Mesh, OnEdge, ExchangeHalos::No);

   // Compute numerical results for vertex variables

   parallelForOuter(
       {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          VorticityAux.computeVarsOnEdge(Team, IEdge);
       });
   const auto &NumNormRelVortEdge    = VorticityAux.NormRelVortEdge;
   const auto &NumNormPlanetVortEdge = VorticityAux.NormPlanetVortEdge;

   // Compute error measures and check errors for edge variables

   ErrorMeasures NormRelVortEdgeErrors;
   Err += computeErrors(NormRelVortEdgeErrors, NumNormRelVortEdge,
                        ExactNormRelVortEdge, Mesh, OnEdge);
   Err += checkErrors("AuxVarsTest", "NormRelVortEdge", NormRelVortEdgeErrors,
                      Setup.ExpectedNormRelVortEdgeErrors, RTol);

   ErrorMeasures NormPlanetVortEdgeErrors;
   Err += computeErrors(NormPlanetVortEdgeErrors, NumNormPlanetVortEdge,
                        ExactNormPlanetVortEdge, Mesh, OnEdge);
   Err += checkErrors("AuxVarsTest", "NormPlanetVortEdge",
                      NormPlanetVortEdgeErrors,
                      Setup.ExpectedNormPlanetVortEdgeErrors, RTol);

   if (Err == 0) {
      LOG_INFO("AuxVarsTest: VorticityAuxVars PASS");
   }

   return Err;
}

int testVelocityDel2AuxVars(Real RTol) {
   TestSetup Setup;
   int Err = 0;

   const auto Decomp = Decomp::getDefault();
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();
   VelocityDel2AuxVars VelocityDel2Aux("", Mesh, VCoord);

   // Use analytical expressions to compute inputs

   Array2DReal ExactVelocityDivCell("ExactVelocityDivCell", Mesh->NCellsSize,
                                    NVertLayers);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.divergence(X, Y); },
       ExactVelocityDivCell, Geom, Mesh, OnCell);

   Array2DReal ExactRelVortVertex("ExactRelVortVertex", Mesh->NVerticesSize,
                                  NVertLayers);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.relativeVorticity(X, Y); },
       ExactRelVortVertex, Geom, Mesh, OnVertex);

   // Compute exact Del2

   Array2DReal ExactDel2Edge("ExactDel2Edge", Mesh->NEdgesOwned, NVertLayers);
   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.velocityDel2X(X, Y);
          VecField[1] = Setup.velocityDel2Y(X, Y);
       },
       ExactDel2Edge, EdgeComponent::Normal, Geom, Mesh, ExchangeHalos::No);

   // Compute numerical Del2

   parallelForOuter(
       {Decomp->NEdgesHaloH(1)},
       KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          VelocityDel2Aux.computeVarsOnEdge(Team, IEdge, ExactVelocityDivCell,
                                            ExactRelVortVertex);
       });
   const auto &NumDel2Edge = VelocityDel2Aux.Del2Edge;

   // Compute error measures and check errors for Del2

   ErrorMeasures Del2Errors;
   Err += computeErrors(Del2Errors, NumDel2Edge, ExactDel2Edge, Mesh, OnEdge);

   Err += checkErrors("AuxVarsTest", "Del2", Del2Errors,
                      Setup.ExpectedDel2Errors, RTol);

   // Compute exact Del2Div

   Array2DReal ExactDel2DivCell("ExactDel2DivCell", Mesh->NCellsOwned,
                                NVertLayers);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.velocityDel2Div(X, Y); },
       ExactDel2DivCell, Geom, Mesh, OnCell, ExchangeHalos::No);

   // Compute numerical Del2Div

   parallelForOuter(
       LaunchConfig({Mesh->NCellsOwned}, TeamScratch<Real>(NVertLayers)),
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          VelocityDel2Aux.computeVarsOnCell(Team, ICell);
       });
   const auto &NumDel2DivCell = VelocityDel2Aux.Del2DivCell;

   // Compute error measures and check errors for Del2Div

   ErrorMeasures Del2DivErrors;
   Err += computeErrors(Del2DivErrors, NumDel2DivCell, ExactDel2DivCell, Mesh,
                        OnCell);
   Err += checkErrors("AuxVarsTest", "Del2Div", Del2DivErrors,
                      Setup.ExpectedDel2DivErrors, RTol);

   // Compute exact Del2RelVort

   Array2DReal ExactDel2RelVortVertex("ExactDel2RelVortVertex",
                                      Mesh->NVerticesOwned, NVertLayers);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.velocityDel2Curl(X, Y); },
       ExactDel2RelVortVertex, Geom, Mesh, OnVertex, ExchangeHalos::No);

   // Compute numerical Del2RelVort

   parallelForOuter(
       LaunchConfig({Mesh->NVerticesOwned}, TeamScratch<Real>(NVertLayers)),
       KOKKOS_LAMBDA(int IVertex, const TeamMember &Team) {
          VelocityDel2Aux.computeVarsOnVertex(Team, IVertex);
       });
   const auto &NumDel2RelVortVertex = VelocityDel2Aux.Del2RelVortVertex;

   // Compute error measures and check errors for Del2RelVort

   ErrorMeasures Del2RelVortErrors;
   Err += computeErrors(Del2RelVortErrors, NumDel2RelVortVertex,
                        ExactDel2RelVortVertex, Mesh, OnVertex);
   Err += checkErrors("AuxVarsTest", "Del2RelVort", Del2RelVortErrors,
                      Setup.ExpectedDel2RelVortErrors, RTol);

   if (Err == 0) {
      LOG_INFO("AuxVarsTest: VelocityDel2AuxVars PASS");
   }

   return Err;
}

int testTracerAuxVars(const Array2DReal &PseudoThickCell,
                      const Array2DReal &NormalVelEdge, Real RTol) {

   TestSetup Setup;
   int Err = 0;

   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   TracerAuxVars TracerAux("", Mesh, VCoord, NTracers);

   // Set input arrays

   Array3DReal TracersOnCell("TracersOnCell", NTracers, Mesh->NCellsSize,
                             NVertLayers);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.tracer(X, Y); },
       TracersOnCell, Geom, Mesh, OnCell);

   Array2DReal PseudoThickEdge("PseudoThickEdge", Mesh->NEdgesSize,
                               NVertLayers);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.pseudoThickness(X, Y); },
       PseudoThickEdge, Geom, Mesh, OnEdge);

   // Compute exact Del2TracerCell

   Array3DReal ExactDel2TrCell("ExactDel2TrCell", NTracers, Mesh->NCellsOwned,
                               NVertLayers);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.del2Tracer(X, Y); },
       ExactDel2TrCell, Geom, Mesh, OnCell, ExchangeHalos::No);

   // Compute numerical Del2TracerCell

   parallelForOuter(
       LaunchConfig({NTracers, Mesh->NCellsOwned},
                    TeamScratch<Real>(NVertLayers)),
       KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
          TracerAux.computeVarsOnCells(Team, L, ICell, PseudoThickEdge,
                                       TracersOnCell);
       });

   // Compute error measures and check errors for Del2TracersCell

   const auto &NumDel2TrCell = TracerAux.Del2TracersCell;

   ErrorMeasures Del2TracerErrors;
   Err += computeErrors(Del2TracerErrors, NumDel2TrCell, ExactDel2TrCell, Mesh,
                        OnCell);
   Err += checkErrors("AuxVarsTest", "Del2Tracers", Del2TracerErrors,
                      Setup.ExpectedDel2TracerErrors, RTol);

   if (Err == 0) {
      LOG_INFO("AuxVarsTest: TracerAuxVars PASS");
   }

   return Err;
}
//------------------------------------------------------------------------------
// The initialization routine for aux vars testing
int initAuxVarsTest(const std::string &mesh) {
   int Err = 0;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // initialize logging
   initLogging(DefEnv);

   // Open config file
   Config("Omega");
   OMEGA::Config::readAll("omega.yml");

   IO::init(DefComm);
   Decomp::init(mesh);

   int HaloErr = Halo::init();
   if (HaloErr != 0) {
      Err++;
      LOG_ERROR("AuxVarsTest: error initializing default halo");
   }

   // Create dummy clock for IO
   Calendar::init("No Leap");
   TimeInstant StartTime(0, 1, 1, 0, 0, 0.0);
   TimeInterval TimeStep(1, TimeUnits::Hours);
   Clock ModelClockTmp(StartTime, TimeStep);
   Clock *ModelClock = &ModelClockTmp;

   // Read horizontal mesh
   Field::init(ModelClock);
   IOStream::init(ModelClock);
   HorzMesh::init(ModelClock);

   // initialize vertical coordinate, do not read stream and use local
   // NVertLayers value
   VertCoord::init(false, NVertLayers);

   return Err;
}

void finalizeAuxVarsTest() {
   VertCoord::clear();
   HorzMesh::clear();
   Field::clear();
   Dimension::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

int auxVarsTest(const std::string &mesh = DefaultMeshFile) {
   int Err = initAuxVarsTest(mesh);
   if (Err != 0) {
      LOG_CRITICAL("AuxVarsTest: Error initializing");
   }

   const auto &Mesh = HorzMesh::getDefault();

   Array2DReal PseudoThickCell("PseudoThickCell", Mesh->NCellsSize,
                               NVertLayers);
   Array2DReal NormalVelEdge("NormalVelEdge", Mesh->NEdgesSize, NVertLayers);
   Err += initState(PseudoThickCell, NormalVelEdge, Mesh);

   const Real RTol = sizeof(Real) == 4 ? 1e-2 : 2e-4;

   Err += testKineticAuxVars(PseudoThickCell, NormalVelEdge, RTol);

   Err += testPseudoThicknessAuxVars(PseudoThickCell, NormalVelEdge, RTol);

   Err += testVorticityAuxVars(PseudoThickCell, NormalVelEdge, RTol);

   Err += testVelocityDel2AuxVars(RTol);

   Err += testTracerAuxVars(PseudoThickCell, NormalVelEdge, RTol);

   if (Err == 0) {
      LOG_INFO("AuxVarsTest: Successful completion");
   }
   finalizeAuxVarsTest();

   return Err;
}

int main(int argc, char *argv[]) {

   int RetVal = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   RetVal += auxVarsTest();

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
