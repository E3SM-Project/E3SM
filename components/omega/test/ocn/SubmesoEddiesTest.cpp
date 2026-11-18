#include "SubmesoEddies.h"
#include "FillValues.h"
#include "GlobalConstants.h"
#include "IO.h"
#include "IOStream.h"
#include "OceanDriver.h"
#include "OceanTestCommon.h"
#include "Pacer.h"
#include "TimeStepper.h"
#include "auxiliaryVars/MixedLayerAuxVars.h"

#include "mpi.h"
#include <limits>
#include <utility>

using namespace OMEGA;

constexpr Geometry Geom = Geometry::Spherical;

KOKKOS_FUNCTION Real verticalSpacing(Real Lon, Real Lat, Real Scale) {
   return Scale * 10 *
          (0.5 + 0.3 * Kokkos::sin(Lon) * Kokkos::sin(Lon) * Kokkos::cos(Lat));
}

KOKKOS_FUNCTION Real surfaceHeight(Real Lon, Real Lat) {
   return -20 * Kokkos::sin(Lon) * Kokkos::cos(Lat);
}

KOKKOS_FUNCTION Real geomHeightMid(int K, Real Lon, Real Lat, Real Scale) {
   const Real Z0 = surfaceHeight(Lon, Lat);
   const Real DZ = verticalSpacing(Lon, Lat, Scale);
   return Z0 - (K + 1) * DZ + DZ / 2;
}

KOKKOS_FUNCTION Real geomHeightInterface(int K, Real Lon, Real Lat,
                                         Real Scale) {
   const Real Z0 = surfaceHeight(Lon, Lat);
   const Real DZ = verticalSpacing(Lon, Lat, Scale);
   return Z0 - K * DZ;
}

constexpr Real H0    = 2e4;
constexpr Real Alpha = 0.4;

KOKKOS_FUNCTION Real specVol(Real Lon, Real Lat, Real Z) {
   return 1 / RhoSw *
          (1 + Alpha * Kokkos::cos(Lon) * Kokkos::pow(Kokkos::cos(Lat), 4)) *
          Kokkos::exp(-(Z * Z) / (H0 * H0));
}

KOKKOS_FUNCTION Real gradBuoyancyX(Real Lon, Real Lat, Real Z) {
   return Alpha * Gravity / REarth * Kokkos::sin(Lon) *
          Kokkos::pow(Kokkos::cos(Lat), 3) * Kokkos::exp(-(Z * Z) / (H0 * H0));
}

KOKKOS_FUNCTION Real gradBuoyancyY(Real Lon, Real Lat, Real Z) {
   return Alpha * Gravity / REarth * 4 * Kokkos::cos(Lon) *
          Kokkos::pow(Kokkos::cos(Lat), 3) * Kokkos::sin(Lat) *
          Kokkos::exp(-(Z * Z) / (H0 * H0));
}

// Buoyancy gradient averaged between Z1 and Z2
KOKKOS_FUNCTION Real meanGradBuoyancyX(Real Lon, Real Lat, Real Z1, Real Z2) {
   const Real Tmp = Alpha * Gravity / REarth * Kokkos::sin(Lon) *
                    Kokkos::pow(Kokkos::cos(Lat), 3);
   return Tmp * H0 * Kokkos::sqrt(Pi) / 2 *
          (Kokkos::erf(Z2 / H0) - Kokkos::erf(Z1 / H0)) / (Z2 - Z1);
}

KOKKOS_FUNCTION Real meanGradBuoyancyY(Real Lon, Real Lat, Real Z1, Real Z2) {
   const Real Tmp = Alpha * Gravity / REarth * 4 * Kokkos::cos(Lon) *
                    Kokkos::pow(Kokkos::cos(Lat), 3) * Kokkos::sin(Lat);
   return Tmp * H0 * Kokkos::sqrt(Pi) / 2 *
          (Kokkos::erf(Z2 / H0) - Kokkos::erf(Z1 / H0)) / (Z2 - Z1);
}

KOKKOS_FUNCTION Real bruntVaisalaFreqSq(Real Lon, Real Lat, Real Z) {
   return Gravity * (-2 * Z / (H0 * H0));
}

// Brunt-Vaisala frequency averaged between Z1 and Z2
KOKKOS_FUNCTION Real meanBVML(Real Lon, Real Lat, Real Z1, Real Z2) {
   const Real Tmp = Gravity * (2 / (H0 * H0));
   return 2._Real / 3 * Kokkos::sqrt(Tmp) *
          (Z2 * Kokkos::sqrt(-Z2) - Z1 * Kokkos::sqrt(-Z1)) / (Z2 - Z1);
}

KOKKOS_FUNCTION Real shapeFunctionDeriv(Real Z, Real H) {
   const Real Tmp  = (2 * Z) / H + 1;
   const Real Fac1 = (1 - Tmp * Tmp);
   const Real Fac2 = (1 + 5 * Tmp * Tmp / 21);

   if (Fac1 * Fac2 >= 0) {
      return -4 / H * Tmp * Fac2 + 20 / (21 * H) * Tmp * Fac1;
   } else {
      return 0;
   }
}

Error initSubmesoEddiesTest(MPI_Comm Comm, const std::string MeshFile,
                            int NVertLayers) {
   Error Err;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // Initialize the Logging system
   initLogging(DefEnv);

   LOG_INFO("------ SubmesoEddies Unit Tests ------");

   // Open config file
   Config("Omega");
   Config::readAll("omega.yml");

   // First step of time stepper initialization needed for IOstream
   TimeStepper::init1();

   // Get model clock
   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   // Initialize the IO system
   IO::init(DefComm);

   // Create the default decomposition (initializes the decomposition)
   Decomp::init(MeshFile);

   // Initialize streams
   IOStream::init();

   // Initialize the default halo
   Halo::init();

   // Initialize the default mesh
   HorzMesh::init(ModelClock);

   // Initialize the default vertical coordinate
   if (NVertLayers == 0) {
      VertCoord::init();
   } else {
      VertCoord::init(false, NVertLayers);
   }

   return Err;
}

void finalizeSubmesoEddiesTest() {

   IOStream::finalize();
   Tracers::clear();
   VertAdv::clear();
   VertCoord::clear();
   TimeStepper::clear();
   HorzMesh::clear();
   Field::clear();
   Dimension::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

std::pair<Array2DReal, Array2DReal> setupVerticalCoord() {
   auto *Mesh   = HorzMesh::getDefault();
   auto *VCoord = VertCoord::getDefault();

   const auto &LonCellH = Mesh->LonCellH;
   const auto &LatCellH = Mesh->LatCellH;
   auto LonCell         = createDeviceMirrorCopy(LonCellH);
   auto LatCell         = createDeviceMirrorCopy(LatCellH);

   const int NVertLayers    = VCoord->NVertLayers;
   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;
   auto &GeomZMid           = VCoord->GeomZMid;
   auto &GeomZInterface     = VCoord->GeomZInterface;

   const Real Scale = 64.0 / NVertLayers;

   // Set vertical coordinates
   parallelForOuter(
       {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          const Real Lon = LonCell(ICell);
          const Real Lat = LatCell(ICell);

          const int KMin = MinLayerCell(ICell);
          const int KMax = MaxLayerCell(ICell);

          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                 GeomZMid(ICell, K) = geomHeightMid(K, Lon, Lat, Scale);
              });

          parallelForInner(
              Team, Range{KMin, KMax + 1}, INNER_LAMBDA(int K) {
                 GeomZInterface(ICell, K) =
                     geomHeightInterface(K, Lon, Lat, Scale);
              });
       });

   // for this test we also need vertical coordinates on edges
   Array2DReal GeomZMidEdge("GeomZMidEdge", Mesh->NEdgesSize,
                            VCoord->NVertLayers);
   Array2DReal GeomZInterfaceEdge("GeomZInterfaceEdge", Mesh->NEdgesSize,
                                  VCoord->NVertLayersP1);

   const auto &LonEdgeH = Mesh->LonEdgeH;
   const auto &LatEdgeH = Mesh->LatEdgeH;
   auto LonEdge         = createDeviceMirrorCopy(LonEdgeH);
   auto LatEdge         = createDeviceMirrorCopy(LatEdgeH);

   const auto &MinLayerEdgeTop = VCoord->MinLayerEdgeTop;
   const auto &MaxLayerEdgeBot = VCoord->MaxLayerEdgeBot;

   parallelForOuter(
       {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          const Real Lon = LonEdge(IEdge);
          const Real Lat = LatEdge(IEdge);

          const int KMin = MinLayerEdgeTop(IEdge);
          const int KMax = MaxLayerEdgeBot(IEdge);

          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                 GeomZMidEdge(IEdge, K) = geomHeightMid(K, Lon, Lat, Scale);
              });

          parallelForInner(
              Team, Range{KMin, KMax + 1}, INNER_LAMBDA(int K) {
                 GeomZInterfaceEdge(IEdge, K) =
                     geomHeightInterface(K, Lon, Lat, Scale);
              });
       });

   return {GeomZMidEdge, GeomZInterfaceEdge};
}

Array2DReal computePseudoThickOnEdges() {
   auto *Mesh              = HorzMesh::getDefault();
   const auto &CellsOnEdge = Mesh->CellsOnEdge;

   auto *VCoord                = VertCoord::getDefault();
   const auto &MinLayerCell    = VCoord->MinLayerCell;
   const auto &MaxLayerCell    = VCoord->MaxLayerCell;
   const auto &MinLayerEdgeBot = VCoord->MinLayerEdgeBot;
   const auto &MaxLayerEdgeTop = VCoord->MaxLayerEdgeTop;
   const auto &GeomZInterface  = VCoord->GeomZInterface;

   Array2DReal PseudoThick("PseudoThick", Mesh->NCellsSize,
                           VCoord->NVertLayers);

   parallelForOuter(
       {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          const int KMin = MinLayerCell(ICell);
          const int KMax = MaxLayerCell(ICell);
          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                 // The unit tests assume that pseudo thickness is equal to
                 // geometric thickness
                 PseudoThick(ICell, K) =
                     GeomZInterface(ICell, K) - GeomZInterface(ICell, K + 1);
              });
       });

   Array2DReal MeanPseudoThickEdge("MeanPseudoThickEdge", Mesh->NEdgesSize,
                                   VCoord->NVertLayers);

   parallelForOuter(
       {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          const int KMin = MinLayerEdgeBot(IEdge);
          const int KMax = MaxLayerEdgeTop(IEdge);

          const int JCell0 = CellsOnEdge(IEdge, 0);
          const int JCell1 = CellsOnEdge(IEdge, 1);

          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                 MeanPseudoThickEdge(IEdge, K) =
                     (PseudoThick(JCell0, K) + PseudoThick(JCell1, K)) / 2;
              });
       });

   return MeanPseudoThickEdge;
}

Array2DReal computeExactGradInterface(const Array2DReal &GeomZInterfaceEdge) {
   auto *Mesh                  = HorzMesh::getDefault();
   auto *VCoord                = VertCoord::getDefault();
   const auto &MinLayerEdgeBot = VCoord->MinLayerEdgeBot;
   const auto &MaxLayerEdgeTop = VCoord->MaxLayerEdgeTop;

   Array1DI4 MaxLayerEdgeTopP1("MaxLayerEdgeTopP1", Mesh->NEdgesSize);

   parallelFor(
       {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge) {
          MaxLayerEdgeTopP1(IEdge) = MaxLayerEdgeTop(IEdge) + 1;
       });

   Array2DReal ExactGradInterface("ExactGradInterface", Mesh->NEdgesSize,
                                  VCoord->NVertLayersP1);
   deepCopy(ExactGradInterface, FillValueReal);

   setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], int IEdge, int K, Real Lon, Real Lat,
                     Real Z) {
          VecField[0] = gradBuoyancyX(Lon, Lat, Z);
          VecField[1] = gradBuoyancyY(Lon, Lat, Z);
       },
       ExactGradInterface, EdgeComponent::Normal, Geom, Mesh, MinLayerEdgeBot,
       MaxLayerEdgeTopP1, GeomZInterfaceEdge, ExchangeHalos::No,
       CartProjection::No);

   // Boundary edges should have fill value
   parallelFor(
       {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge) {
          int MinLyrEdgeBot = MinLayerEdgeBot(IEdge);
          int MaxLyrEdgeTop = MaxLayerEdgeTop(IEdge);
          if (MaxLyrEdgeTop < MinLyrEdgeBot) {
             ExactGradInterface(IEdge, MinLyrEdgeBot) = FillValueReal;
          }
       });

   return ExactGradInterface;
}

Array1DReal computeExactGradML(const Array2DReal &GeomZInterfaceEdge,
                               const Array1DReal &DenMixLayerDepth) {
   auto *Mesh              = HorzMesh::getDefault();
   const auto &CellsOnEdge = Mesh->CellsOnEdge;

   auto *VCoord                = VertCoord::getDefault();
   const auto &MinLayerEdgeBot = VCoord->MinLayerEdgeBot;

   Array1DReal ExactGradML("ExactGradML", Mesh->NEdgesSize);

   setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], int IEdge, Real Lon, Real Lat) {
          const int JCell0 = CellsOnEdge(IEdge, 0);
          const int JCell1 = CellsOnEdge(IEdge, 1);
          const Real MLDepthEdge =
              Kokkos::min(DenMixLayerDepth(JCell0), DenMixLayerDepth(JCell1));

          const Real Z0 = GeomZInterfaceEdge(IEdge, MinLayerEdgeBot(IEdge));

          VecField[0] = meanGradBuoyancyX(Lon, Lat, Z0, Z0 - MLDepthEdge);
          VecField[1] = meanGradBuoyancyY(Lon, Lat, Z0, Z0 - MLDepthEdge);
       },
       ExactGradML, EdgeComponent::Normal, Geom, Mesh, ExchangeHalos::No,
       CartProjection::No);

   return ExactGradML;
}

Array1DReal computeExactBVML(const Array2DReal &GeomZInterfaceEdge,
                             const Array1DReal &DenMixLayerDepth) {
   auto *Mesh              = HorzMesh::getDefault();
   const auto &CellsOnEdge = Mesh->CellsOnEdge;

   auto *VCoord                = VertCoord::getDefault();
   const auto &MinLayerEdgeBot = VCoord->MinLayerEdgeBot;

   Array1DReal ExactBVML("ExactBVML", Mesh->NEdgesSize);

   setScalar(
       KOKKOS_LAMBDA(int IEdge, Real Lon, Real Lat) {
          const int JCell0 = CellsOnEdge(IEdge, 0);
          const int JCell1 = CellsOnEdge(IEdge, 1);
          const Real MLDepthEdge =
              Kokkos::min(DenMixLayerDepth(JCell0), DenMixLayerDepth(JCell1));

          const Real Z0 = GeomZInterfaceEdge(IEdge, MinLayerEdgeBot(IEdge));

          return meanBVML(Lon, Lat, Z0, Z0 - MLDepthEdge);
       },
       ExactBVML, Geom, Mesh, OnEdge, ExchangeHalos::No);

   return ExactBVML;
}

// If SetupExact is true then this test is set up such that ML depth is located
// exactly on layer midpoints. Since no interpolation is necessary,
// very low error norms are expected.
Error testDenMixedLayerDepth(bool SetupExact) {
   Error Err;

   auto *Mesh = HorzMesh::getDefault();

   auto *VCoord             = VertCoord::getDefault();
   const int NVertLayers    = VCoord->NVertLayers;
   const auto &MinLayerCell = VCoord->MinLayerCell;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;
   auto &GeomZMid           = VCoord->GeomZMid;
   auto &GeomZInterface     = VCoord->GeomZInterface;

   MixedLayerAuxVars MixedLayerAux("", Mesh, VCoord);

   Array2DReal ReferenceSpecVol("ReferenceSpecVol", Mesh->NCellsAll,
                                VCoord->NVertLayers);
   Array1DReal ExactMixLayerDepth("ExactMixLayerDepth", Mesh->NCellsAll);

   const Real DenThreshold = MixedLayerAux.DenThreshold;

   parallelForOuter(
       {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          const int KMin = MinLayerCell(ICell);
          const int KMax = MaxLayerCell(ICell);

          const Real MaxDepth =
              GeomZInterface(ICell, KMin) - GeomZMid(ICell, KMax);

          const int KDen = Kokkos::max(10, ICell % NVertLayers);

          const bool KDenValid = KDen <= KMax;

          const Real DepthKDen =
              KDenValid ? GeomZInterface(ICell, KMin) - GeomZMid(ICell, KDen)
                        : MaxDepth;

          const Real ExactDepth = SetupExact ? DepthKDen : 0.95 * DepthKDen;

          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                 const Real Depth =
                     GeomZInterface(ICell, KMin) - GeomZMid(ICell, K);

                 Real Dens;
                 if (SetupExact) {
                    Dens = K < KDen ? RhoSw : RhoSw + (1 + 1e-6) * DenThreshold;
                 } else {
                    Dens = RhoSw + DenThreshold * (Depth / ExactDepth) *
                                       (Depth / ExactDepth);
                 }

                 ReferenceSpecVol(ICell, K) = 1._Real / Dens;
              });

          ExactMixLayerDepth(ICell) = KDenValid ? ExactDepth : MaxDepth;
       });

   parallelForOuter(
       {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          MixedLayerAux.computeVarsOnCell(Team, ICell, ReferenceSpecVol);
       });
   const auto &NumMixLayerDepth = MixedLayerAux.DenMixLayerDepth;

   ErrorMeasures MixLayerDepthErrors;
   computeErrors(MixLayerDepthErrors, NumMixLayerDepth, ExactMixLayerDepth,
                 Mesh, OnCell);

   if (SetupExact &&
       (MixLayerDepthErrors.L2 > 1e-7 || MixLayerDepthErrors.LInf > 1e-7)) {
      Err += Error(ErrorCode::Fail, "denMixedLayerDepth Exact FAIL");
   }

   if (!SetupExact &&
       (MixLayerDepthErrors.L2 > 5e-2 || MixLayerDepthErrors.LInf > 5e-2)) {
      Err += Error(ErrorCode::Fail, "denMixedLayerDepth NonExact FAIL");
   }

   return Err;
}

Error testBuoyancyGrad(const Array2DReal &MeanPseudoThickEdge,
                       const Array2DReal &GeomZInterfaceEdge) {
   Error Err;

   auto *Mesh = HorzMesh::getDefault();

   auto *VCoord             = VertCoord::getDefault();
   auto &GeomZMid           = VCoord->GeomZMid;
   auto &GeomZInterface     = VCoord->GeomZInterface;
   const auto &MaxLayerCell = VCoord->MaxLayerCell;

   // Compute specific volume
   Array2DReal SpecVol("SpecVol", Mesh->NCellsSize, VCoord->NVertLayers);
   setScalar(
       KOKKOS_LAMBDA(int ICell, int K, Real Lon, Real Lat, Real Z) {
          return specVol(Lon, Lat, Z);
       },
       SpecVol, Geom, Mesh, OnCell, VCoord->MinLayerCell, VCoord->MaxLayerCell,
       GeomZMid);

   // Compute squared Brunt-Vaisala frequency
   Array1DI4 MaxLayerCellP1("MaxLayerCellP1", Mesh->NCellsSize);
   parallelFor(
       {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell) {
          MaxLayerCellP1(ICell) = MaxLayerCell(ICell) + 1;
       });

   Array2DReal BVFreqSq("BVFreqSq", Mesh->NCellsSize, VCoord->NVertLayersP1);
   setScalar(
       KOKKOS_LAMBDA(int ICell, int K, Real Lon, Real Lat, Real Z) {
          return bruntVaisalaFreqSq(Lon, Lat, Z);
       },
       BVFreqSq, Geom, Mesh, OnCell, VCoord->MinLayerCell, MaxLayerCellP1,
       GeomZInterface);

   // Compute numerical buoyancy gradient
   auto *SubEddies = SubmesoEddies::getInstance();

   SubEddies->computeBuoyGrad(SpecVol, MeanPseudoThickEdge, GeomZMid, BVFreqSq);

   // Compute exact buoyancy gradient
   Array2DReal ExactGradInterface =
       computeExactGradInterface(GeomZInterfaceEdge);

   // Compute errors and check they are reasonable
   ErrorMeasures InterfaceGradErrors;
   computeErrors(InterfaceGradErrors, SubEddies->GradBuoyEdgeInterface,
                 ExactGradInterface, Mesh, OnEdge);

   const Real MaxL2Error = 4.1e-4;
   if (InterfaceGradErrors.L2 > MaxL2Error) {
      Err += Error(ErrorCode::Fail, "buoyancyGrad L2 FAIL, {:e} > {:e}",
                   InterfaceGradErrors.L2, MaxL2Error);
   }

   const Real MaxLInfError = 5.5e-4;
   if (InterfaceGradErrors.LInf > MaxLInfError) {
      Err += Error(ErrorCode::Fail, "buoyancyGrad LInf FAIL, {:e} > {:e}",
                   InterfaceGradErrors.LInf, MaxLInfError);
   }

   return Err;
}

Error testEddyVelocity(const Array2DReal &GeomZInterfaceEdge,
                       const Array2DReal &MeanPseudoThickEdge) {
   Error Err;

   auto *Mesh               = HorzMesh::getDefault();
   const auto &CellsOnEdge  = Mesh->CellsOnEdge;
   const auto &CellsOnCell  = Mesh->CellsOnCell;
   const auto &NEdgesOnCell = Mesh->NEdgesOnCell;
   const auto &DcEdge       = Mesh->DcEdge;
   const int NCellsAll      = Mesh->NCellsAll;

   auto *VCoord                = VertCoord::getDefault();
   const auto &MinLayerCell    = VCoord->MinLayerCell;
   const auto &MaxLayerCell    = VCoord->MaxLayerCell;
   const auto &MinLayerEdgeBot = VCoord->MinLayerEdgeBot;
   const auto &MaxLayerEdgeTop = VCoord->MaxLayerEdgeTop;
   auto &GeomZMid              = VCoord->GeomZMid;
   auto &GeomZInterface        = VCoord->GeomZInterface;

   auto *SubEddies = SubmesoEddies::getInstance();

   const Real LfMin = 3500;
   SubEddies->LfMin = LfMin;
   const Real DsMax = SubEddies->DsMax;

   Array1DReal DenMixLayerDepth("DenMixLayerDepth", Mesh->NCellsSize);
   Array1DI4 DenMixLayerIndex("DenMixLayerIndex", Mesh->NCellsSize);
   const auto &TimeScale = SubEddies->TimeScale;
   const Real Ce         = SubEddies->Ce;

   // Pick initial mixed layer indices arbitrarily
   Array1DI4 DenMixLayerIndexTmp("DenMixLayerIndexTmp", Mesh->NCellsSize);
   parallelFor(
       {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell) {
          const int KMax = MaxLayerCell(ICell);

          int MLIndex                = Kokkos::min(KMax, 15 + ICell % 10);
          DenMixLayerIndexTmp(ICell) = MLIndex;
       });

   // Make sure there are no large differences in mixed layer depth between
   // neighbouring cells
   parallelFor(
       {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell) {
          const int KMin = MinLayerCell(ICell);

          int MLIndex = DenMixLayerIndexTmp(ICell);
          for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
             const int JCell = CellsOnCell(ICell, J);
             if (JCell != NCellsAll) {
                MLIndex = Kokkos::min(MLIndex, DenMixLayerIndexTmp(JCell));
             }
          }
          DenMixLayerIndex(ICell) = MLIndex;
          DenMixLayerDepth(ICell) =
              GeomZInterface(ICell, KMin) - GeomZMid(ICell, MLIndex);
       });

   // Compute exact buoyancy gradient
   SubEddies->GradBuoyEdgeInterface =
       computeExactGradInterface(GeomZInterfaceEdge);

   Array2DReal BVFreqSq("BVFreqSq", Mesh->NCellsSize, VCoord->NVertLayersP1);

   // Compute exact Brunt-Vaisala frequency
   Array1DI4 MaxLayerCellP1("MaxLayerCellP1", Mesh->NCellsSize);
   parallelFor(
       {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell) {
          MaxLayerCellP1(ICell) = MaxLayerCell(ICell) + 1;
       });

   setScalar(
       KOKKOS_LAMBDA(int ICell, int K, Real Lon, Real Lat, Real Z) {
          return bruntVaisalaFreqSq(Lon, Lat, Z);
       },
       BVFreqSq, Geom, Mesh, OnCell, VCoord->MinLayerCell, MaxLayerCellP1,
       GeomZInterface);

   // Compute numerical eddy velocity
   SubEddies->computeEddyVelocity(DenMixLayerDepth, DenMixLayerIndex, BVFreqSq,
                                  MeanPseudoThickEdge);

   const auto &EddyVelocity = SubEddies->EddyVelocity;

   // Compute exact mixed layer average of buoyancy gradient
   Array1DReal MeanBuoyGrad =
       computeExactGradML(GeomZInterfaceEdge, DenMixLayerDepth);

   // Compute exact mixed layer average of Brunt-Vaisala frequency
   Array1DReal MeanBV = computeExactBVML(GeomZInterfaceEdge, DenMixLayerDepth);

   // Compute exact eddy velocity
   Array2DReal ExactEddyVelocity("ExactEddyVelocity", Mesh->NEdgesSize,
                                 VCoord->NVertLayers);
   deepCopy(ExactEddyVelocity, FillValueReal);

   parallelForOuter(
       {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          const int MinLyrEdgeBot = MinLayerEdgeBot(IEdge);
          const int MaxLyrEdgeTop = MaxLayerEdgeTop(IEdge);

          const int JCell0 = CellsOnEdge(IEdge, 0);
          const int JCell1 = CellsOnEdge(IEdge, 1);

          const Real TScale = TimeScale(IEdge);

          const Real MLDepthEdge =
              Kokkos::min(DenMixLayerDepth(JCell0), DenMixLayerDepth(JCell1));

          const Real Ds = Kokkos::min(DcEdge(IEdge), DsMax);

          const Real GradBuoyML = MeanBuoyGrad(IEdge);
          const Real BVFreqML   = MeanBV(IEdge);

          const Real Lf1 =
              Kokkos::abs(GradBuoyML) * MLDepthEdge / (TScale * TScale);
          const Real Lf2 = BVFreqML * MLDepthEdge / TScale;

          const Real Lf = Kokkos::max(LfMin, Kokkos::max(Lf1, Lf2));

          const Real Factor =
              Ce * Ds / Lf * MLDepthEdge * MLDepthEdge * GradBuoyML / TScale;

          if (MaxLyrEdgeTop >= MinLyrEdgeBot) {

             parallelForInner(
                 Team, Range{MinLyrEdgeBot, MaxLyrEdgeTop},
                 INNER_LAMBDA(int K) {
                    const Real ZEdge =
                        0.5_Real *
                        (GeomZMid(JCell0, K) -
                         GeomZInterface(JCell0, MinLayerCell(JCell0)) +
                         GeomZMid(JCell1, K) -
                         GeomZInterface(JCell1, MinLayerCell(JCell1)));

                    ExactEddyVelocity(IEdge, K) =
                        -Factor * shapeFunctionDeriv(ZEdge, MLDepthEdge);
                 });
          }
       });

   // Compute errors and check that they are reasonable
   ErrorMeasures EddyVelocityErrors;
   computeErrors(EddyVelocityErrors, EddyVelocity, ExactEddyVelocity, Mesh,
                 OnEdge);

   const Real MaxL2Error = 2.4e-1;
   if (EddyVelocityErrors.L2 > MaxL2Error) {
      Err += Error(ErrorCode::Fail, "eddyVelocity L2 FAIL, {:e} > {:e}",
                   EddyVelocityErrors.L2, MaxL2Error);
   }

   const Real MaxLInfError = 5.3e-1;
   if (EddyVelocityErrors.LInf > MaxLInfError) {
      Err += Error(ErrorCode::Fail, "eddyVelocity LInf FAIL, {:e} > {:e}",
                   EddyVelocityErrors.LInf, MaxLInfError);
   }

   return Err;
}

Error testSubmesoEddies(MPI_Comm Comm, std::string MeshFile, int NVertLayers) {
   Error Err;

   // Initialize Omega modules needed for this test
   Err += initSubmesoEddiesTest(Comm, MeshFile, NVertLayers);

   // Setup vertical coordinates
   auto [GeomZMidEdge, GeomZInterfaceEdge] = setupVerticalCoord();

   // Initialize submesoscale eddy parametrization
   SubmesoEddies::init();

   // Test retrieval
   if (!SubmesoEddies::getInstance()) {
      ABORT_ERROR("SubmesoEddiesTest: SubmesoEddies retrieval FAIL");
   }

   // Test mixed layer depth computation
   for (bool SetupExact : {true, false}) {
      Err += testDenMixedLayerDepth(SetupExact);
   }

   // Compute pseudo-thickness on edges for the subsequent tests
   auto MeanPseudoThickEdge = computePseudoThickOnEdges();

   // Test buoyancy gradient computation
   Err += testBuoyancyGrad(MeanPseudoThickEdge, GeomZInterfaceEdge);

   // Test eddy velocity computation
   Err += testEddyVelocity(GeomZInterfaceEdge, MeanPseudoThickEdge);

   // Destroy submesoscale eddy parametrization
   SubmesoEddies::destroyInstance();

   finalizeSubmesoEddiesTest();

   return Err;
}

int main(int argc, char *argv[]) {
   Error Err;

   const MPI_Comm Comm = MPI_COMM_WORLD;

   MPI_Init(&argc, &argv);
   Pacer::initialize(Comm);
   Pacer::setPrefix("Omega:");

   try {
      Kokkos::initialize(argc, argv);
      { Err += testSubmesoEddies(Comm, "OmegaMesh.nc", 0); }
      Kokkos::finalize();
   } catch (const std::exception &Ex) {
      Err += Error(ErrorCode::Fail, Ex.what() + std::string(": FAIL"));
   } catch (...) {
      Err += Error(ErrorCode::Fail, "Unknown: FAIL");
   }

   CHECK_ERROR_ABORT(Err, "Submeso Eddies Unit Tests FAIL");

   Pacer::finalize();
   MPI_Finalize();

   return 0;

} // end of main
//===-----------------------------------------------------------------------===/
