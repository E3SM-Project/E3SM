#include "DataTypes.h"
#include "Decomp.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOField.h"
#include "Logging.h"
#include "MachEnv.h"
#include "MetaData.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "auxiliaryVars/KineticAuxVars.h"
#include "auxiliaryVars/LayerThicknessAuxVars.h"
#include "auxiliaryVars/VorticityAuxVars.h"
#include "mpi.h"

#include <cmath>
#include <iomanip>

using namespace OMEGA;

struct TestSetupPlane {
   Real Pi = M_PI;

   Real Lx = 1;
   Real Ly = std::sqrt(3) / 2;

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

   KOKKOS_FUNCTION Real layerThickness(Real X, Real Y) const {
      return 2 + std::cos(2 * Pi * X / Lx) * std::cos(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real velocityX(Real X, Real Y) const {
      return std::sin(2 * Pi * X / Lx) * std::cos(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real velocityY(Real X, Real Y) const {
      return std::cos(2 * Pi * X / Lx) * std::sin(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real divergence(Real X, Real Y) const {
      return 2 * Pi * (1. / Lx + 1. / Ly) * std::cos(2 * Pi * X / Lx) *
             std::cos(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real relativeVorticity(Real X, Real Y) const {
      return 2 * Pi * (-1. / Lx + 1. / Ly) * std::sin(2 * Pi * X / Lx) *
             std::sin(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real planetaryVorticity(Real X, Real Y) const {
      return std::sin(2 * Pi * X / Lx) * std::sin(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real normalizedRelativeVorticity(Real X, Real Y) const {
      return relativeVorticity(X, Y) / layerThickness(X, Y);
   }

   KOKKOS_FUNCTION Real normalizedPlanetaryVorticity(Real X, Real Y) const {
      return planetaryVorticity(X, Y) / layerThickness(X, Y);
   }

   KOKKOS_FUNCTION Real kineticEnergy(Real X, Real Y) const {
      return (velocityX(X, Y) * velocityX(X, Y) +
              velocityY(X, Y) * velocityY(X, Y)) /
             2;
   }
};

struct TestSetupSphere {
   Real Radius = 6371220;

   ErrorMeasures ExpectedKineticEnergyErrors = {0.0143579382532765844,
                                                0.00681096618897046764};
   ErrorMeasures ExpectedVelocityDivErrors   = {0.0136595773989793799,
                                                0.00367052484586382699};

   ErrorMeasures ExpectedFluxThickErrors = {0.0159821090867812224,
                                            0.010364511516135164};
   ErrorMeasures ExpectedMeanThickErrors = {0.000800109287518277435,
                                            0.000406527457820634436};

   ErrorMeasures ExpectedRelVortVertexErrors        = {0.0271404735181343393,
                                                       0.0252023166109219786};
   ErrorMeasures ExpectedNormRelVortVertexErrors    = {0.0348741350737879693,
                                                       0.0259506101504540822};
   ErrorMeasures ExpectedNormPlanetVortVertexErrors = {0.00451268952953497778,
                                                       0.00101771171197261793};

   ErrorMeasures ExpectedNormRelVortEdgeErrors    = {0.0125376497261775952,
                                                     0.00307521304930552519};
   ErrorMeasures ExpectedNormPlanetVortEdgeErrors = {0.00495174534686814403,
                                                     0.000855432390947949515};

   KOKKOS_FUNCTION Real layerThickness(Real Lon, Real Lat) const {
      return (2 + std::cos(Lon) * std::pow(std::cos(Lat), 4));
   }

   KOKKOS_FUNCTION Real velocityX(Real Lon, Real Lat) const {
      return -Radius * std::pow(std::sin(Lon), 2) * std::pow(std::cos(Lat), 3);
   }

   KOKKOS_FUNCTION Real velocityY(Real Lon, Real Lat) const {
      return -4 * Radius * std::sin(Lon) * std::cos(Lon) *
             std::pow(std::cos(Lat), 3) * std::sin(Lat);
   }

   KOKKOS_FUNCTION Real relativeVorticity(Real Lon, Real Lat) const {
      return -4 * std::pow(std::cos(Lon), 2) * std::pow(std::cos(Lat), 2) *
             std::sin(Lat);
   }

   KOKKOS_FUNCTION Real divergence(Real Lon, Real Lat) const {
      return std::sin(Lon) * std::cos(Lon) * std::pow(std::cos(Lat), 2) *
             (20 * std::pow(std::sin(Lat), 2) - 6);
   }

   KOKKOS_FUNCTION Real planetaryVorticity(Real Lon, Real Lat) const {
      return std::sin(Lat);
   }

   KOKKOS_FUNCTION Real normalizedRelativeVorticity(Real Lon, Real Lat) const {
      return relativeVorticity(Lon, Lat) / layerThickness(Lon, Lat);
   }

   KOKKOS_FUNCTION Real normalizedPlanetaryVorticity(Real Lon, Real Lat) const {
      return planetaryVorticity(Lon, Lat) / layerThickness(Lon, Lat);
   }

   KOKKOS_FUNCTION Real kineticEnergy(Real Lon, Real Lat) const {
      return (velocityX(Lon, Lat) * velocityX(Lon, Lat) +
              velocityY(Lon, Lat) * velocityY(Lon, Lat)) /
             2;
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

constexpr int NVertLevels = 16;

int initState(const Array2DReal &LayerThickCell,
              const Array2DReal &NormalVelEdge, HorzMesh *Mesh) {
   int Err = 0;

   TestSetup Setup;

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.layerThickness(X, Y); },
       LayerThickCell, Geom, Mesh, OnCell, NVertLevels);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.velocityX(X, Y);
          VecField[1] = Setup.velocityY(X, Y);
       },
       NormalVelEdge, EdgeComponent::Normal, Geom, Mesh, NVertLevels);

   // need to override FVertex with prescribed values
   // cannot use setScalar because it doesn't support setting 1D arrays
   const auto &FVertex = Mesh->FVertex;

   auto XVertex = createDeviceMirrorCopy(Mesh->XVertexH);
   auto YVertex = createDeviceMirrorCopy(Mesh->YVertexH);

   auto LonVertex = createDeviceMirrorCopy(Mesh->LonVertexH);
   auto LatVertex = createDeviceMirrorCopy(Mesh->LatVertexH);

   parallelFor(
       {Mesh->NVerticesOwned}, KOKKOS_LAMBDA(int IVertex) {
          if (Geom == Geometry::Planar) {
             const Real XV    = XVertex(IVertex);
             const Real YV    = YVertex(IVertex);
             FVertex(IVertex) = Setup.planetaryVorticity(XV, YV);
          } else {
             const Real XV    = LonVertex(IVertex);
             const Real YV    = LatVertex(IVertex);
             FVertex(IVertex) = Setup.planetaryVorticity(XV, YV);
          }
       });

   auto MyHalo    = Halo::getDefault();
   auto &FVertexH = Mesh->FVertexH;
   deepCopy(FVertexH, FVertex);
   Err += MyHalo->exchangeFullArrayHalo(FVertexH, OnVertex);
   deepCopy(FVertex, FVertexH);

   return Err;
}

int testKineticAuxVars(const Array2DReal &LayerThicknessCell,
                       const Array2DReal &NormalVelocityEdge, Real RTol) {
   int Err = 0;
   TestSetup Setup;

   const auto Mesh = HorzMesh::getDefault();

   // Compute exact result

   Array2DReal ExactKineticEnergyCell("ExactKineticEnergyCell",
                                      Mesh->NCellsOwned, NVertLevels);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.kineticEnergy(X, Y); },
       ExactKineticEnergyCell, Geom, Mesh, OnCell, NVertLevels, false);

   Array2DReal ExactVelocityDivCell("ExactVelocityDivCell", Mesh->NCellsOwned,
                                    NVertLevels);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.divergence(X, Y); },
       ExactVelocityDivCell, Geom, Mesh, OnCell, NVertLevels, false);

   // Compute numerical result

   KineticAuxVars KineticAux(Mesh, NVertLevels);
   parallelFor(
       {Mesh->NCellsOwned, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
          KineticAux.computeVarsOnCell(ICell, KLevel, NormalVelocityEdge);
       });
   const auto &NumKineticEnergyCell = KineticAux.KineticEnergyCell;
   const auto &NumVelocityDivCell   = KineticAux.VelocityDivCell;

   // Compute error measures and check error values

   ErrorMeasures KineticEnergyErrors;
   Err += computeErrors(KineticEnergyErrors, NumKineticEnergyCell,
                        ExactKineticEnergyCell, Mesh, OnCell, NVertLevels);
   Err += checkErrors("AuxVarsTest", "KineticEnergy", KineticEnergyErrors,
                      Setup.ExpectedKineticEnergyErrors, RTol);

   ErrorMeasures VelocityDivErrors;
   Err += computeErrors(VelocityDivErrors, NumVelocityDivCell,
                        ExactVelocityDivCell, Mesh, OnCell, NVertLevels);
   Err += checkErrors("AuxVarsTest", "VelocityDiv", VelocityDivErrors,
                      Setup.ExpectedVelocityDivErrors, RTol);

   if (Err == 0) {
      LOG_INFO("AuxVarsTest: KineticAuxVars PASS");
   }

   return Err;
}

int testLayerThicknessAuxVars(const Array2DReal &LayerThickCell,
                              const Array2DReal &NormalVelEdge, Real RTol) {
   int Err = 0;
   TestSetup Setup;

   const auto Mesh = HorzMesh::getDefault();

   // Compute exact result

   Array2DReal ExactThickEdge("ExactThickEdge", Mesh->NEdgesOwned, NVertLevels);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.layerThickness(X, Y); },
       ExactThickEdge, Geom, Mesh, OnEdge, NVertLevels, false);

   // Compute numerical result

   LayerThicknessAuxVars LayerThicknessAux(Mesh, NVertLevels);
   parallelFor(
       {Mesh->NEdgesOwned, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
          LayerThicknessAux.computeVarsOnEdge(IEdge, KLevel, LayerThickCell,
                                              NormalVelEdge);
       });

   const auto &NumFluxLayerThickEdge = LayerThicknessAux.FluxLayerThickEdge;
   const auto &NumMeanLayerThickEdge = LayerThicknessAux.MeanLayerThickEdge;

   // Compute error measures and check error values

   ErrorMeasures FluxThickErrors;
   Err += computeErrors(FluxThickErrors, NumFluxLayerThickEdge, ExactThickEdge,
                        Mesh, OnEdge, NVertLevels);
   Err += checkErrors("AuxVarsTest", "FluxThick", FluxThickErrors,
                      Setup.ExpectedFluxThickErrors, RTol);

   ErrorMeasures MeanThickErrors;
   Err += computeErrors(MeanThickErrors, NumMeanLayerThickEdge, ExactThickEdge,
                        Mesh, OnEdge, NVertLevels);
   Err += checkErrors("AuxVarsTest", "MeanThick", MeanThickErrors,
                      Setup.ExpectedMeanThickErrors, RTol);

   if (Err == 0) {
      LOG_INFO("AuxVarsTest: ThickAuxVars PASS");
   }

   return Err;
}

int testVorticityAuxVars(const Array2DReal &LayerThickCell,
                         const Array2DReal &NormalVelEdge, Real RTol) {
   TestSetup Setup;
   int Err = 0;

   const auto Decomp = Decomp::getDefault();
   const auto Mesh   = HorzMesh::getDefault();

   VorticityAuxVars VorticityAux(Mesh, NVertLevels);

   // Compute exact results for vertex variables

   Array2DReal ExactRelVortVertex("ExactRelVortVertex", Mesh->NVerticesOwned,
                                  NVertLevels);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.relativeVorticity(X, Y); },
       ExactRelVortVertex, Geom, Mesh, OnVertex, NVertLevels, false);

   Array2DReal ExactNormRelVortVertex("ExactNormRelVortVertex",
                                      Mesh->NVerticesOwned, NVertLevels);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) {
          return Setup.normalizedRelativeVorticity(X, Y);
       },
       ExactNormRelVortVertex, Geom, Mesh, OnVertex, NVertLevels, false);

   Array2DReal ExactNormPlanetVortVertex("ExactNormPlanetVortVertex",
                                         Mesh->NVerticesOwned, NVertLevels);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) {
          return Setup.normalizedPlanetaryVorticity(X, Y);
       },
       ExactNormPlanetVortVertex, Geom, Mesh, OnVertex, NVertLevels, false);

   // Compute numerical results for vertex variables

   parallelFor(
       {Decomp->NVerticesHaloH(0), NVertLevels},
       KOKKOS_LAMBDA(int IVertex, int KLevel) {
          VorticityAux.computeVarsOnVertex(IVertex, KLevel, LayerThickCell,
                                           NormalVelEdge);
       });

   const auto &NumRelVortVertex        = VorticityAux.RelVortVertex;
   const auto &NumNormRelVortVertex    = VorticityAux.NormRelVortVertex;
   const auto &NumNormPlanetVortVertex = VorticityAux.NormPlanetVortVertex;

   // Compute error measures and check errors for vertex variables

   ErrorMeasures RelVortVertexErrors;
   Err += computeErrors(RelVortVertexErrors, NumRelVortVertex,
                        ExactRelVortVertex, Mesh, OnVertex, NVertLevels);
   Err += checkErrors("AuxVarsTest", "RelVortVertex", RelVortVertexErrors,
                      Setup.ExpectedRelVortVertexErrors, RTol);

   ErrorMeasures NormRelVortVertexErrors;
   Err += computeErrors(NormRelVortVertexErrors, NumNormRelVortVertex,
                        ExactNormRelVortVertex, Mesh, OnVertex, NVertLevels);
   Err +=
       checkErrors("AuxVarsTest", "NormRelVortVertex", NormRelVortVertexErrors,
                   Setup.ExpectedNormRelVortVertexErrors, RTol);

   ErrorMeasures NormPlanetVortVertexErrors;
   Err += computeErrors(NormPlanetVortVertexErrors, NumNormPlanetVortVertex,
                        ExactNormPlanetVortVertex, Mesh, OnVertex, NVertLevels);
   Err += checkErrors("AuxVarsTest", "NormPlanetVortVertex",
                      NormPlanetVortVertexErrors,
                      Setup.ExpectedNormPlanetVortVertexErrors, RTol);

   // Compute exact results for edge variables

   Array2DReal ExactNormRelVortEdge("ExactNormRelVortEdge", Mesh->NEdgesOwned,
                                    NVertLevels);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) {
          return Setup.normalizedRelativeVorticity(X, Y);
       },
       ExactNormRelVortEdge, Geom, Mesh, OnEdge, NVertLevels, false);

   Array2DReal ExactNormPlanetVortEdge("ExactNormPlanetVortEdge",
                                       Mesh->NEdgesOwned, NVertLevels);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) {
          return Setup.normalizedPlanetaryVorticity(X, Y);
       },
       ExactNormPlanetVortEdge, Geom, Mesh, OnEdge, NVertLevels, false);

   // Compute numerical results for vertex variables

   parallelFor(
       {Mesh->NEdgesOwned, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
          VorticityAux.computeVarsOnEdge(IEdge, KLevel);
       });
   const auto &NumNormRelVortEdge    = VorticityAux.NormRelVortEdge;
   const auto &NumNormPlanetVortEdge = VorticityAux.NormPlanetVortEdge;

   // Compute error measures and check errors for edge variables

   ErrorMeasures NormRelVortEdgeErrors;
   Err += computeErrors(NormRelVortEdgeErrors, NumNormRelVortEdge,
                        ExactNormRelVortEdge, Mesh, OnEdge, NVertLevels);
   Err += checkErrors("AuxVarsTest", "NormRelVortEdge", NormRelVortEdgeErrors,
                      Setup.ExpectedNormRelVortEdgeErrors, RTol);

   ErrorMeasures NormPlanetVortEdgeErrors;
   Err += computeErrors(NormPlanetVortEdgeErrors, NumNormPlanetVortEdge,
                        ExactNormPlanetVortEdge, Mesh, OnEdge, NVertLevels);
   Err += checkErrors("AuxVarsTest", "NormPlanetVortEdge",
                      NormPlanetVortEdgeErrors,
                      Setup.ExpectedNormPlanetVortEdgeErrors, RTol);

   if (Err == 0) {
      LOG_INFO("AuxVarsTest: VorticityAuxVars PASS");
   }

   return Err;
}

//------------------------------------------------------------------------------
// The initialization routine for aux vars testing
int initAuxVarsTest(const std::string &mesh) {
   int Err = 0;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefaultEnv();
   MPI_Comm DefComm = DefEnv->getComm();

   int IOErr = IO::init(DefComm);
   if (IOErr != 0) {
      Err++;
      LOG_ERROR("AuxVarsTest: error initializing parallel IO");
   }

   int DecompErr = Decomp::init(mesh);
   if (DecompErr != 0) {
      Err++;
      LOG_ERROR("AuxVarsTest: error initializing default decomposition");
   }

   int HaloErr = Halo::init();
   if (HaloErr != 0) {
      Err++;
      LOG_ERROR("AuxVarsTest: error initializing default halo");
   }

   int MeshErr = HorzMesh::init();
   if (MeshErr != 0) {
      Err++;
      LOG_ERROR("AuxVarsTest: error initializing default mesh");
   }

   const auto &Mesh = HorzMesh::getDefault();
   MetaDim::create("NCells", Mesh->NCellsSize);
   MetaDim::create("NVertices", Mesh->NVerticesSize);
   MetaDim::create("NEdges", Mesh->NEdgesSize);
   MetaDim::create("NVertLevels", NVertLevels);
   MetaGroup::create("auxiliaryVars");

   return Err;
}

void finalizeAuxVarsTest() {
   IOField::clear();
   HorzMesh::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

void auxVarsTest(const std::string &mesh = DefaultMeshFile) {
   int Err = initAuxVarsTest(mesh);
   if (Err != 0) {
      LOG_CRITICAL("AuxVarsTest: Error initializing");
   }

   const auto &Mesh = HorzMesh::getDefault();

   Array2DReal LayerThickCell("LayerThickCell", Mesh->NCellsSize, NVertLevels);
   Array2DReal NormalVelEdge("NormalVelEdge", Mesh->NEdgesSize, NVertLevels);
   Err += initState(LayerThickCell, NormalVelEdge, Mesh);

   const Real RTol = sizeof(Real) == 4 ? 1e-2 : 1e-10;

   Err += testKineticAuxVars(LayerThickCell, NormalVelEdge, RTol);
   Err += testLayerThicknessAuxVars(LayerThickCell, NormalVelEdge, RTol);
   Err += testVorticityAuxVars(LayerThickCell, NormalVelEdge, RTol);

   if (Err == 0) {
      LOG_INFO("AuxVarsTest: Successful completion");
   }
   finalizeAuxVarsTest();
}

int main(int argc, char *argv[]) {
   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);

   auxVarsTest();

   Kokkos::finalize();
   MPI_Finalize();
} // end of main
//===-----------------------------------------------------------------------===/
