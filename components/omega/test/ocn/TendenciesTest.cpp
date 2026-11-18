#include "Tendencies.h"
#include "AuxiliaryState.h"
#include "Config.h"
#include "CustomTendencyTerms.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Eos.h"
#include "Error.h"
#include "Field.h"
#include "Forcing.h"
#include "GlobalConstants.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "KPPMix.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "PGrad.h"
#include "Pacer.h"
#include "TimeStepper.h"
#include "VertCoord.h"
#include "VertMix.h"
#include "mpi.h"

#include <cmath>
#include <iomanip>

using namespace OMEGA;

struct TestSetup {
   // radius of spherical mesh
   Real Radius = REarth;

   KOKKOS_FUNCTION Real pseudoThickness(Real Lon, Real Lat) const {
      return (2 + std::cos(Lon) * std::pow(std::cos(Lat), 4));
   }

   KOKKOS_FUNCTION Real velocityX(Real Lon, Real Lat) const {
      return -Radius * std::pow(std::sin(Lon), 2) * std::pow(std::cos(Lat), 3);
   }

   KOKKOS_FUNCTION Real velocityY(Real Lon, Real Lat) const {
      return -4 * Radius * std::sin(Lon) * std::cos(Lon) *
             std::pow(std::cos(Lat), 3) * std::sin(Lat);
   }

   KOKKOS_FUNCTION Real tracer(Real Lon, Real Lat) const {
      return (2 - std::cos(Lon) * std::pow(std::cos(Lat), 4));
   }
};

constexpr Geometry Geom   = Geometry::Spherical;
constexpr int NVertLayers = 60;

int testSfcTracerForcing();
int testSfcThicknessForcing();
int testKPPNonLocalVerticalMixGate();

int initState() {
   int Err = 0;

   TestSetup Setup;
   auto *Mesh   = HorzMesh::getDefault();
   auto *VCoord = VertCoord::getDefault();
   auto *State  = OceanState::getDefault();

   // Define tendency fields
   int NDims = 2;
   std::vector<std::string> DimNamesThickness(NDims);
   DimNamesThickness[0] = "NCells";

   Array2DReal PseudoThickCell = State->getPseudoThickness(0);
   Array2DReal NormalVelEdge   = State->getNormalVelocity(0);

   Array3DReal TracersArray = Tracers::getAll(0);
   const auto &TracersCell  = TracersArray;

   int NTracers = Tracers::getNumTracers();

   deepCopy(PseudoThickCell, NAN);
   deepCopy(NormalVelEdge, NAN);
   deepCopy(TracersCell, NAN);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.pseudoThickness(X, Y); },
       PseudoThickCell, Geom, Mesh, OnCell, VCoord->MinLayerCell,
       VCoord->MaxLayerCell, ExchangeHalos::Yes, SetBoundary::Yes);

   deepCopy(VCoord->SurfacePressure, 0.0_Real);
   VCoord->computePressure(PseudoThickCell, VCoord->SurfacePressure);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.tracer(X, Y); },
       TracersCell, Geom, Mesh, OnCell, VCoord->MinLayerCell,
       VCoord->MaxLayerCell, ExchangeHalos::Yes, SetBoundary::Yes);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real Lon, Real Lat) {
          VecField[0] = Setup.velocityX(Lon, Lat);
          VecField[1] = Setup.velocityY(Lon, Lat);
       },
       NormalVelEdge, EdgeComponent::Normal, Geom, Mesh,
       VCoord->MinLayerEdgeTop, VCoord->MaxLayerEdgeBot, ExchangeHalos::Yes,
       CartProjection::No, SetBoundary::Yes);

   return Err;
}

//------------------------------------------------------------------------------
// The initialization routine for tendencies testing
int initTendenciesTest(const std::string &mesh) {
   int Err = 0;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   initLogging(DefEnv);

   // Open config file
   Config("Omega");
   Config::readAll("omega.yml");

   Config *OmegaConfig = Config::getOmegaConfig();
   Config TendConfig("Tendencies");
   OmegaConfig->get(TendConfig);
   TendConfig.set("FrazilTendencyEnable", true);

   Config TimeIntConfig("TimeIntegration");
   OmegaConfig->get(TimeIntConfig);
   TimeIntConfig.set("TimeStepper", std::string("RungeKutta4"));

   // Initialize time stepping and model clock
   TimeStepper::init1();
   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   IO::init(DefComm);
   Decomp::init(mesh);

   // Initialize streams
   Field::init(ModelClock);
   IOStream::init(ModelClock);

   int HaloErr = Halo::init();
   if (HaloErr != 0) {
      Err++;
      LOG_ERROR("TendenciesTest: error initializing default halo");
   }

   // Read mesh
   HorzMesh::init(ModelClock);
   VertCoord::init();

   // SurfacePressure holds the attachData fill value until this defaults it
   // to zero (no initial-state read in this test) and exchanges halos.
   VertCoord::getDefault()->initSurfacePressure(Halo::getDefault());

   Tracers::init();
   Frazil::init();
   VertAdv::init();
   PressureGrad::init();
   Eos::init();
   Forcing::init();
   VertMix::init();

   int StateErr = OceanState::init();
   if (StateErr != 0) {
      Err++;
      LOG_ERROR("TendenciesTest: error initializing default state");
   }

   AuxiliaryState::init();

   return Err;
}

int testTendencies() {
   int Err = 0;
   Error Err1;

   // test initialization
   Tendencies::init();

   // test retrievel of default
   Tendencies *DefTendencies = Tendencies::getDefault();

   if (DefTendencies) {
      LOG_INFO("TendenciesTest: Default tendencies retrieval PASS");
   } else {
      Err++;
      LOG_INFO("TendenciesTest: Default tendencies retrieval FAIL");
      return -1;
   }

   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   const auto VAdv     = VertAdv::getDefault();
   const auto PGrad    = PressureGrad::getDefault();
   const auto EqState  = Eos::getInstance();
   const auto VMix     = VertMix::getInstance();
   VCoord->NVertLayers = 12;

   // test creation of another tendencies

   TimeInterval ZeroTimeStep; // Zero-length time step placeholder
   Config *Options = Config::getOmegaConfig();
   Config TendConfig("Tendencies");
   Err1             = Options->get(TendConfig);
   int NTracersTest = 3;

   Tendencies::create("TestTendencies", Mesh, VCoord, VAdv, PGrad, EqState,
                      VMix, NTracersTest, ZeroTimeStep, &TendConfig);

   // test retrievel of another tendencies
   if (Tendencies::get("TestTendencies")) {
      LOG_INFO("TendenciesTest: Non-default tendencies retrieval PASS");
   } else {
      Err++;
      LOG_INFO("TendenciesTest: Non-default tendencies retrieval FAIL");
   }

   // test erase
   Tendencies::erase("TestTendencies");

   if (Tendencies::get("TestTendencies")) {
      Err++;
      LOG_INFO("TendenciesTest: Non-default tendencies erase FAIL");
   } else {
      LOG_INFO("TendenciesTest: Non-default tendencies erase PASS");
   }

   VCoord->NVertLayers = NVertLayers;

   // put NANs in every tendency variables
   deepCopy(DefTendencies->PseudoThicknessTend, NAN);
   deepCopy(DefTendencies->NormalVelocityTend, NAN);
   deepCopy(DefTendencies->TracerTend, NAN);

   // compute tendencies
   const auto *State       = OceanState::getDefault();
   const auto *AuxState    = AuxiliaryState::getDefault();
   auto *DefForcing        = Forcing::getDefault();
   Array3DReal TracerArray = Tracers::getAll(0);
   int ThickTimeLevel      = 0;
   int VelTimeLevel        = 0;
   int TracerTimeLevel     = 0;
   TimeInstant Time;
   TimeInterval Interval(1., TimeUnits::Seconds);

   if (DefForcing == nullptr) {
      LOG_ERROR("TendenciesTest: Default forcing retrieval FAIL");
      Tendencies::clear();
      return -1;
   }

   auto &ZonalStressCell  = DefForcing->SfcStressForcing.ZonalStressCell;
   auto &MeridStressCell  = DefForcing->SfcStressForcing.MeridStressCell;
   auto &NormalStressEdge = DefForcing->SfcStressForcing.NormalStressEdge;

   const bool OrigSfcStressEnabled = DefTendencies->SfcStressForcing.Enabled;
   Array2DReal BaselineNormalVelocityTend(
       "BaselineNormalVelocityTend", Mesh->NEdgesSize, VCoord->NVertLayers);

   DefTendencies->SfcStressForcing.Enabled = false;
   DefTendencies->FrazilTerm.Enabled = false; // frazil needs to be off for now
   DefTendencies->computeAllTendencies(State, AuxState, TracerArray,
                                       ThickTimeLevel, VelTimeLevel,
                                       TracerTimeLevel, Time, Interval);
   deepCopy(BaselineNormalVelocityTend, DefTendencies->NormalVelocityTend);

   deepCopy(DefTendencies->PseudoThicknessTend, NAN);
   deepCopy(DefTendencies->NormalVelocityTend, NAN);
   deepCopy(DefTendencies->TracerTend, NAN);

   DefTendencies->SfcStressForcing.Enabled = true;
   deepCopy(ZonalStressCell, 1._Real);
   deepCopy(MeridStressCell, 0.5_Real);
   deepCopy(NormalStressEdge, NAN);
   DefForcing->computeAll();

   DefTendencies->computeAllTendencies(State, AuxState, TracerArray,
                                       ThickTimeLevel, VelTimeLevel,
                                       TracerTimeLevel, Time, Interval);

   Array2DReal NormalVelocityTendDiff("NormalVelocityTendDiff",
                                      Mesh->NEdgesSize, VCoord->NVertLayers);
   deepCopy(NormalVelocityTendDiff, 0._Real);

   OMEGA_SCOPE(LocNormalVelocityTendDiff, NormalVelocityTendDiff);
   OMEGA_SCOPE(LocBaselineNormalVelocityTend, BaselineNormalVelocityTend);
   OMEGA_SCOPE(LocNormalVelocityTend, DefTendencies->NormalVelocityTend);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   parallelForOuter(
       "TendenciesTest:NormalVelocityTendDiff", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          const int KMin   = MinLayerEdgeBot(IEdge);
          const int KMax   = MaxLayerEdgeTop(IEdge);
          const int KRange = vertRangeChunked(KMin, KMax);

          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 for (int K = KChunk; K <= KMax; K += VecLength) {
                    if (K >= KMin) {
                       LocNormalVelocityTendDiff(IEdge, K) =
                           Kokkos::abs(LocNormalVelocityTend(IEdge, K) -
                                       LocBaselineNormalVelocityTend(IEdge, K));
                    }
                 }
              });
       });

   const Real NormVelTendDelta =
       sum(NormalVelocityTendDiff, Mesh->NEdgesOwned, VCoord->MinLayerEdgeBot,
           VCoord->MaxLayerEdgeTop);
   constexpr Real DeltaATol = 1e-12_Real;
   if (!Kokkos::isfinite(NormVelTendDelta) ||
       isApprox(NormVelTendDelta, 0._Real, 0._Real, DeltaATol)) {
      Err++;
      LOG_ERROR("TendenciesTest: SfcStress forcing did not change "
                "NormalVelocityTend");
   }

   const Real NormVelTendSum =
       sum(DefTendencies->NormalVelocityTend, Mesh->NEdgesOwned,
           VCoord->MinLayerEdgeBot, VCoord->MaxLayerEdgeTop);
   if (!Kokkos::isfinite(NormVelTendSum) || NormVelTendSum == 0) {
      Err++;
      LOG_ERROR("TendenciesTest: NormVelTendSum FAIL");
   }

   DefTendencies->SfcStressForcing.Enabled = OrigSfcStressEnabled;

   // Test surface tracer forcing with enthalpy terms (TEOS-10 CtFrz path)
   const int TracerForcingErr = testSfcTracerForcing();
   Err += TracerForcingErr;

   // Test surface thickness forcing with freshwater terms
   const int ThicknessForcingErr = testSfcThicknessForcing();
   Err += ThicknessForcingErr;

   // Verify frazil tendencies are plumbed into pseudo-thickness and tracer
   // tendencies by comparing runs with FrazilTerm disabled/enabled.
   I4 TempTracerIndex = -1;
   I4 SaltTracerIndex = -1;
   Tracers::getIndex(TempTracerIndex, "Temperature");
   Tracers::getIndex(SaltTracerIndex, "Salinity");

   if (TempTracerIndex < 0 || SaltTracerIndex < 0) {
      Err++;
      LOG_ERROR("TendenciesTest: missing Temperature/Salinity tracer indices");
   } else {
      auto MinLayerCellH = createHostMirrorCopy(VCoord->MinLayerCell);
      auto MaxLayerCellH = createHostMirrorCopy(VCoord->MaxLayerCell);
      auto TracerArrayH  = createHostMirrorCopy(TracerArray);

      // Set Salinity in active layers to [31, 33], top-to-bottom ramp.
      for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
         const I4 KMin = MinLayerCellH(ICell);
         const I4 KMax = MaxLayerCellH(ICell);

         for (I4 K = KMin; K <= KMax; ++K) {
            Real frac = 0.0_Real;
            if (KMax > KMin) {
               frac =
                   static_cast<Real>(K - KMin) / static_cast<Real>(KMax - KMin);
            }
            TracerArrayH(SaltTracerIndex, ICell, K) =
                Real(31.0) + Real(2.0) * frac;
         }
      }

      // Optional: keep Temperature warm everywhere first, then add one cold
      // point.
      for (I4 ICell = 0; ICell < Mesh->NCellsAll; ++ICell) {
         const I4 KMin = MinLayerCellH(ICell);
         const I4 KMax = MaxLayerCellH(ICell);
         for (I4 K = KMin; K <= KMax; ++K) {
            TracerArrayH(TempTracerIndex, ICell, K) = 8.0_Real;
         }
      }
      const I4 FrazilTestCell  = 0;
      const I4 FrazilTestLayer = MinLayerCellH(FrazilTestCell);
      TracerArrayH(TempTracerIndex, FrazilTestCell, FrazilTestLayer) =
          -2.0_Real;
      deepCopy(TracerArray, TracerArrayH);

      LOG_INFO("TendenciesTest: minLayerCell = {}", FrazilTestLayer);

      Array2DReal BaselinePseudoThicknessTend(
          "BaselinePseudoThicknessTend", Mesh->NCellsSize, VCoord->NVertLayers);
      Array3DReal BaselineTracerTend("BaselineTracerTend",
                                     Tracers::getNumTracers(), Mesh->NCellsSize,
                                     VCoord->NVertLayers);

      const bool OrigFrazilEnabled = DefTendencies->FrazilTerm.Enabled;

      DefTendencies->FrazilTerm.Enabled = false;
      DefTendencies->computeAllTendencies(State, AuxState, TracerArray,
                                          ThickTimeLevel, VelTimeLevel,
                                          TracerTimeLevel, Time, Interval);
      deepCopy(BaselinePseudoThicknessTend, DefTendencies->PseudoThicknessTend);
      deepCopy(BaselineTracerTend, DefTendencies->TracerTend);

      deepCopy(DefTendencies->PseudoThicknessTend, 0.0_Real);
      deepCopy(DefTendencies->NormalVelocityTend, 0.0_Real);
      deepCopy(DefTendencies->TracerTend, 0.0_Real);

      DefTendencies->FrazilTerm.Enabled = true;
      DefTendencies->computeAllTendencies(State, AuxState, TracerArray,
                                          ThickTimeLevel, VelTimeLevel,
                                          TracerTimeLevel, Time, Interval);

      Array2DReal PseudoThicknessTendDiff(
          "PseudoThicknessTendDiff", Mesh->NCellsSize, VCoord->NVertLayers);
      Array2DReal TempTracerTendDiff("TempTracerTendDiff", Mesh->NCellsSize,
                                     VCoord->NVertLayers);
      Array2DReal SaltTracerTendDiff("SaltTracerTendDiff", Mesh->NCellsSize,
                                     VCoord->NVertLayers);

      OMEGA_SCOPE(LocPseudoThicknessTendDiff, PseudoThicknessTendDiff);
      OMEGA_SCOPE(LocTempTracerTendDiff, TempTracerTendDiff);
      OMEGA_SCOPE(LocSaltTracerTendDiff, SaltTracerTendDiff);
      OMEGA_SCOPE(LocPseudoThicknessTend, DefTendencies->PseudoThicknessTend);
      OMEGA_SCOPE(LocBaselinePseudoThicknessTend, BaselinePseudoThicknessTend);
      OMEGA_SCOPE(LocTracerTend, DefTendencies->TracerTend);
      OMEGA_SCOPE(LocBaselineTracerTend, BaselineTracerTend);
      OMEGA_SCOPE(LocMinLayerCell, VCoord->MinLayerCell);
      OMEGA_SCOPE(LocMaxLayerCell, VCoord->MaxLayerCell);

      parallelForOuter(
          "TendenciesTest:FrazilTendDiff", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
             const int KMin = LocMinLayerCell(ICell);
             const int KMax = LocMaxLayerCell(ICell);

             parallelForInner(
                 Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                    LocPseudoThicknessTendDiff(ICell, K) =
                        Kokkos::abs(LocPseudoThicknessTend(ICell, K) -
                                    LocBaselinePseudoThicknessTend(ICell, K));
                    LocTempTracerTendDiff(ICell, K) = Kokkos::abs(
                        LocTracerTend(TempTracerIndex, ICell, K) -
                        LocBaselineTracerTend(TempTracerIndex, ICell, K));
                    LocSaltTracerTendDiff(ICell, K) = Kokkos::abs(
                        LocTracerTend(SaltTracerIndex, ICell, K) -
                        LocBaselineTracerTend(SaltTracerIndex, ICell, K));
                 });
          });

      const Real PseudoDelta = sum(PseudoThicknessTendDiff, Mesh->NCellsOwned,
                                   VCoord->MinLayerCell, VCoord->MaxLayerCell);
      const Real TempDelta   = sum(TempTracerTendDiff, Mesh->NCellsOwned,
                                   VCoord->MinLayerCell, VCoord->MaxLayerCell);
      const Real SaltDelta   = sum(SaltTracerTendDiff, Mesh->NCellsOwned,
                                   VCoord->MinLayerCell, VCoord->MaxLayerCell);

      constexpr Real FrazilDeltaATol = 1e-12_Real;
      if (!Kokkos::isfinite(PseudoDelta) ||
          isApprox(PseudoDelta, 0._Real, 0._Real, FrazilDeltaATol)) {
         Err++;
         LOG_ERROR("TendenciesTest: Frazil did not change "
                   "PseudoThicknessTend");
      }
      if (!Kokkos::isfinite(TempDelta) ||
          isApprox(TempDelta, 0._Real, 0._Real, FrazilDeltaATol)) {
         Err++;
         LOG_ERROR("TendenciesTest: Frazil did not change Temperature tracer "
                   "tendency");
      }
      if (!Kokkos::isfinite(SaltDelta) ||
          isApprox(SaltDelta, 0._Real, 0._Real, FrazilDeltaATol)) {
         Err++;
         LOG_ERROR("TendenciesTest: Frazil did not change Salinity tracer "
                   "tendency");
      }
      // add a log info if no errors?
      DefTendencies->FrazilTerm.Enabled = OrigFrazilEnabled;
   }

   // Verify that the KPP non-local tracer flux is gated by vertical mixing.
   Err += testKPPNonLocalVerticalMixGate();

   // check that everything got computed correctly
   int NCellsOwned = Mesh->NCellsOwned;
   int NTracers    = Tracers::getNumTracers();

   const Real PseudoThickTendSum =
       sum(DefTendencies->PseudoThicknessTend, NCellsOwned,
           VCoord->MinLayerCell, VCoord->MaxLayerCell);
   if (!Kokkos::isfinite(PseudoThickTendSum) || PseudoThickTendSum == 0) {
      Err++;
      LOG_ERROR("TendenciesTest: PseudoThickTend FAIL");
   }

   const Real TracerTendSum =
       sum(DefTendencies->TracerTend, NTracers, NCellsOwned,
           VCoord->MinLayerCell, VCoord->MaxLayerCell);
   if (!Kokkos::isfinite(TracerTendSum) || TracerTendSum == 0) {
      Err++;
      LOG_ERROR("TendenciesTest: TracerTendSum FAIL");
   }

   Tendencies::clear();
   return Err;
}

int testSfcTracerForcing() {
   int Err = 0;

   auto *VCoord        = VertCoord::getDefault();
   auto *DefTendencies = Tendencies::getDefault();
   auto *State         = OceanState::getDefault();
   auto *AuxState      = AuxiliaryState::getDefault();
   auto *DefForcing    = Forcing::getDefault();
   auto *EosInst       = Eos::getInstance();

   Array3DReal TracerArray = Tracers::getAll(0);

   const I4 TempIndex = Tracers::IndxTemp;
   const I4 SaltIndex = Tracers::IndxSalt;

   if (TempIndex < 0 || SaltIndex < 0) {
      LOG_ERROR("TendenciesTest: Invalid tracer indices for SfcTracerForcing");
      return -1;
   }

   deepCopy(DefTendencies->TracerTend, 0._Real);

   // Set up single test cell at top layer
   const I4 ICellTest = 0;
   const I4 KTop      = VCoord->MinLayerCellH(ICellTest);

   if (KTop > VCoord->MaxLayerCellH(ICellTest)) {
      LOG_ERROR("TendenciesTest: Test cell has no layers");
      return -1;
   }

   // Known tracer values for testing
   const Real CtTopValue = 15.0_Real; // °C (conservative temperature)
   const Real SaTopValue = 35.0_Real; // g/kg (salinity)

   // Set tracer values at test cell
   OMEGA_SCOPE(LocTracerArray, TracerArray);
   Kokkos::parallel_for(
       "SetTestTracersForcing", 1, KOKKOS_LAMBDA(int i) {
          LocTracerArray(TempIndex, ICellTest, KTop) = CtTopValue;
          LocTracerArray(SaltIndex, ICellTest, KTop) = SaTopValue;
       });

   // Retrieve forcing field views
   auto &SensibleHeatFlux   = DefForcing->TracerForcing.SensibleHeatFluxCell;
   auto &LatentHeatFluxEvap = DefForcing->TracerForcing.LatentHeatFluxEvapCell;
   auto &LongWaveHeatFluxUp = DefForcing->TracerForcing.LongWaveHeatFluxUpCell;
   auto &LongWaveHeatFluxDown =
       DefForcing->TracerForcing.LongWaveHeatFluxDownCell;
   auto &SeaIceHeatFlux    = DefForcing->TracerForcing.SeaIceHeatFluxCell;
   auto &ShortWaveHeatFlux = DefForcing->TracerForcing.ShortWaveHeatFluxCell;
   auto &RainFlux          = DefForcing->TracerForcing.RainFluxCell;
   auto &RiverRunoffFlux   = DefForcing->TracerForcing.RiverRunoffFluxCell;
   auto &SnowFlux          = DefForcing->TracerForcing.SnowFluxCell;
   auto &IceRunoffFlux     = DefForcing->TracerForcing.IceRunoffFluxCell;
   auto &SeaIceSaltFlux    = DefForcing->TracerForcing.SeaIceSaltFluxCell;

   // Initialize all fluxes to zero
   deepCopy(SensibleHeatFlux, 0._Real);
   deepCopy(LatentHeatFluxEvap, 0._Real);
   deepCopy(LongWaveHeatFluxUp, 0._Real);
   deepCopy(LongWaveHeatFluxDown, 0._Real);
   deepCopy(SeaIceHeatFlux, 0._Real);
   deepCopy(ShortWaveHeatFlux, 0._Real);
   deepCopy(RainFlux, 0._Real);
   deepCopy(RiverRunoffFlux, 0._Real);
   deepCopy(SnowFlux, 0._Real);
   deepCopy(IceRunoffFlux, 0._Real);
   deepCopy(SeaIceSaltFlux, 0._Real);

   // Set test forcing values
   // Non-zero sensible heat: 100 W/m²
   const Real TestSensibleHeat = 100.0_Real;
   // Non-zero rain: 1e-8 kg/m²/s
   const Real TestRain = 1.0e-8_Real;
   // Non-zero snow: 5e-9 kg/m²/s
   const Real TestSnow = 5.0e-9_Real;
   // Sea ice salt flux: 1e-4 kg/m²/s
   const Real TestSeaIceSaltFlux = 1.0e-4_Real;

   OMEGA_SCOPE(LocSensibleHeatFlux, SensibleHeatFlux);
   OMEGA_SCOPE(LocRainFlux, RainFlux);
   OMEGA_SCOPE(LocSnowFlux, SnowFlux);
   OMEGA_SCOPE(LocSeaIceSaltFlux, SeaIceSaltFlux);
   Kokkos::parallel_for(
       "SetTestForcingTracer", 1, KOKKOS_LAMBDA(int i) {
          LocSensibleHeatFlux(ICellTest) = TestSensibleHeat;
          LocRainFlux(ICellTest)         = TestRain;
          LocSnowFlux(ICellTest)         = TestSnow;
          LocSeaIceSaltFlux(ICellTest)   = TestSeaIceSaltFlux;
       });

   DefForcing->computeAll();

   // Disable all tendencies except SfcTracerForcing
   const bool OrigSfcStressEnabled = DefTendencies->SfcStressForcing.Enabled;
   const bool OrigSfcThicknessEnabled =
       DefTendencies->SfcThicknessForcing.Enabled;
   const bool OrigSfcTracerEnabled = DefTendencies->SfcTracerForcing.Enabled;
   const bool OrigPseudoThicknessDiv =
       DefTendencies->PseudoThicknessFluxDiv.Enabled;
   const bool OrigPotentialVortHAdv = DefTendencies->PotentialVortHAdv.Enabled;
   const bool OrigKEGrad            = DefTendencies->KEGrad.Enabled;
   const bool OrigVelocityDiffusion = DefTendencies->VelocityDiffusion.Enabled;
   const bool OrigVelocityHyperDiff = DefTendencies->VelocityHyperDiff.Enabled;
   const bool OrigTracerHorzAdv     = DefTendencies->TracerHorzAdv.Enabled;
   const bool OrigTracerDiffusion   = DefTendencies->TracerDiffusion.Enabled;
   const bool OrigTracerHyperDiff   = DefTendencies->TracerHyperDiff.Enabled;
   const bool OrigSurfaceTracerRestoring =
       DefTendencies->SurfaceTracerRestoring.Enabled;

   DefTendencies->SfcStressForcing.Enabled       = false;
   DefTendencies->SfcThicknessForcing.Enabled    = false;
   DefTendencies->SfcTracerForcing.Enabled       = false;
   DefTendencies->PseudoThicknessFluxDiv.Enabled = false;
   DefTendencies->PotentialVortHAdv.Enabled      = false;
   DefTendencies->KEGrad.Enabled                 = false;
   DefTendencies->VelocityDiffusion.Enabled      = false;
   DefTendencies->VelocityHyperDiff.Enabled      = false;
   DefTendencies->TracerHorzAdv.Enabled          = false;
   DefTendencies->TracerDiffusion.Enabled        = false;
   DefTendencies->TracerHyperDiff.Enabled        = false;
   DefTendencies->SurfaceTracerRestoring.Enabled = false;

   // Compute tendencies
   int ThickTimeLevel  = 0;
   int VelTimeLevel    = 0;
   int TracerTimeLevel = 0;
   TimeInstant Time;
   TimeInterval Interval(1., TimeUnits::Seconds);

   // because vertical advection tendencies are always on, we need to compute a
   // baseline first. the actual test is whether the total tendencies change
   // with the flag toggling.
   DefTendencies->computeAllTendencies(State, AuxState, TracerArray,
                                       ThickTimeLevel, VelTimeLevel,
                                       TracerTimeLevel, Time, Interval);

   HostArray3DReal TracerTendBaseH =
       createHostMirrorCopy(DefTendencies->TracerTend);
   deepCopy(TracerTendBaseH, DefTendencies->TracerTend);
   const Real BaselineTempTend = TracerTendBaseH(TempIndex, ICellTest, KTop);
   const Real BaselineSaltTend = TracerTendBaseH(SaltIndex, ICellTest, KTop);

   // Now enable SfcTracerForcing and compute again
   DefTendencies->SfcTracerForcing.Enabled = true;

   DefTendencies->computeAllTendencies(State, AuxState, TracerArray,
                                       ThickTimeLevel, VelTimeLevel,
                                       TracerTimeLevel, Time, Interval);

   // Build a reference expectations for temperature tendency. Production
   // code approximates the snow/ice enthalpy as the constant -LatIce (see
   // SfcTracerForcingOnCell), so no CtFrz term is expected here.
   const Real ExpectedTempTend =
       (TestSensibleHeat + TestRain * Cp0Sw * CtTopValue - TestSnow * LatIce) *
       HFluxFac;

   // SaltTend = SeaIceSaltFlux * SFluxFac
   const Real ExpectedSaltTend = TestSeaIceSaltFlux * SFluxFac;

   HostArray3DReal TracerTendH =
       createHostMirrorCopy(DefTendencies->TracerTend);
   deepCopy(TracerTendH, DefTendencies->TracerTend);
   const Real ComputedTempTend =
       TracerTendH(TempIndex, ICellTest, KTop) - BaselineTempTend;
   const Real ComputedSaltTend =
       TracerTendH(SaltIndex, ICellTest, KTop) - BaselineSaltTend;

   constexpr Real RelTol = 1.0e-10_Real;
   constexpr Real AbsTol = 1.0e-12_Real; // flux precision is ~e-15

   if (!isApprox(ComputedTempTend, ExpectedTempTend, RelTol, AbsTol)) {
      Err++;
      LOG_ERROR("TendenciesTest: SfcTracerForcing temp tendency FAIL");
      LOG_ERROR("  Expected: {},  Computed: {}, Diff: {}", ExpectedTempTend,
                ComputedTempTend,
                Kokkos::abs(ComputedTempTend - ExpectedTempTend));
   } else {
      LOG_INFO("TendenciesTest: SfcTracerForcing temp tendency PASS");
   }

   // Check salinity tendency
   if (!isApprox(ComputedSaltTend, ExpectedSaltTend, RelTol, AbsTol)) {
      Err++;
      LOG_ERROR("TendenciesTest: SfcTracerForcing salt tendency FAIL");
      LOG_INFO("  Expected: {},  Computed: {}, Diff: {}", ExpectedSaltTend,
               ComputedSaltTend,
               Kokkos::abs(ComputedSaltTend - ExpectedSaltTend));
   } else {
      LOG_INFO("TendenciesTest: SfcTracerForcing salt tendency PASS");
   }

   DefTendencies->SfcStressForcing.Enabled       = OrigSfcStressEnabled;
   DefTendencies->SfcThicknessForcing.Enabled    = OrigSfcThicknessEnabled;
   DefTendencies->SfcTracerForcing.Enabled       = OrigSfcTracerEnabled;
   DefTendencies->PseudoThicknessFluxDiv.Enabled = OrigPseudoThicknessDiv;
   DefTendencies->PotentialVortHAdv.Enabled      = OrigPotentialVortHAdv;
   DefTendencies->KEGrad.Enabled                 = OrigKEGrad;
   DefTendencies->VelocityDiffusion.Enabled      = OrigVelocityDiffusion;
   DefTendencies->VelocityHyperDiff.Enabled      = OrigVelocityHyperDiff;
   DefTendencies->TracerHorzAdv.Enabled          = OrigTracerHorzAdv;
   DefTendencies->TracerDiffusion.Enabled        = OrigTracerDiffusion;
   DefTendencies->TracerHyperDiff.Enabled        = OrigTracerHyperDiff;
   DefTendencies->SurfaceTracerRestoring.Enabled = OrigSurfaceTracerRestoring;

   return Err;
}

int testKPPNonLocalVerticalMixGate() {
   int Err = 0;

   auto *DefTendencies = Tendencies::getDefault();
   auto *State         = OceanState::getDefault();
   auto *AuxState      = AuxiliaryState::getDefault();
   auto *DefForcing    = Forcing::getDefault();
   auto *VMix          = VertMix::getInstance();
   auto *KPPInstance   = KPPMix::getInstance();
   if (!KPPInstance) {
      KPPMix::init();
      KPPInstance = KPPMix::getInstance();
   }
   if (!KPPInstance || !DefTendencies || !State || !AuxState || !DefForcing ||
       !VMix) {
      LOG_ERROR("TendenciesTest: KPP non-local gate setup failed");
      return 1;
   }

   const auto *VCoord      = VertCoord::getDefault();
   const auto *Mesh        = HorzMesh::getDefault();
   Array3DReal TracerArray = Tracers::getAll(0);
   const I4 TempIndex      = Tracers::IndxTemp;
   const I4 ICellTest      = 0;
   const I4 KTest          = VCoord->MinLayerCellH(ICellTest);
   if (TempIndex < 0 || KTest > VCoord->MaxLayerCellH(ICellTest)) {
      LOG_ERROR("TendenciesTest: invalid cell or temperature tracer for KPP "
                "gate test");
      return 1;
   }

   const bool OriginalVelMixEnabled = DefTendencies->VelVertMixTendencyEnable;
   const bool OriginalTracerMixEnabled =
       DefTendencies->TracerVertMixTendencyEnable;
   const bool OriginalVelSetupEnabled    = VMix->VelVertMixSetup.Enabled;
   const bool OriginalTracerSetupEnabled = VMix->TracerVertMixSetup.Enabled;
   const bool OriginalKPPEnabled         = KPPInstance->Enabled;
   const bool OriginalNonLocalEnabled =
       DefTendencies->KPPNonLocalTracerFlux.Enabled;
   const bool OriginalSfcTracerEnabled =
       DefTendencies->SfcTracerForcing.Enabled;

   DefTendencies->VelVertMixTendencyEnable      = false;
   DefTendencies->TracerVertMixTendencyEnable   = false;
   VMix->VelVertMixSetup.Enabled                = false;
   VMix->TracerVertMixSetup.Enabled             = false;
   DefTendencies->KPPNonLocalTracerFlux.Enabled = true;
   DefTendencies->SfcTracerForcing.Enabled      = true;
   KPPInstance->Enabled                         = true;

   DefForcing->resetArrays();
   deepCopy(DefForcing->TracerForcing.SnowFluxCell, 0.0_Real);
   deepCopy(DefForcing->TracerForcing.RainFluxCell, 0.0_Real);
   deepCopy(DefForcing->TracerForcing.EvaporationFluxCell, 0.0_Real);
   deepCopy(DefForcing->TracerForcing.SeaIceFreshWaterFluxCell, 0.0_Real);
   deepCopy(DefForcing->TracerForcing.IceRunoffFluxCell, 0.0_Real);
   deepCopy(DefForcing->TracerForcing.RiverRunoffFluxCell, 0.0_Real);
   deepCopy(DefForcing->TracerForcing.LatentHeatFluxEvapCell, 0.0_Real);
   deepCopy(DefForcing->TracerForcing.SensibleHeatFluxCell, 0.0_Real);
   deepCopy(DefForcing->TracerForcing.LongWaveHeatFluxUpCell, 0.0_Real);
   deepCopy(DefForcing->TracerForcing.LongWaveHeatFluxDownCell, 0.0_Real);
   deepCopy(DefForcing->TracerForcing.SeaIceHeatFluxCell, 0.0_Real);
   deepCopy(DefForcing->TracerForcing.ShortWaveHeatFluxCell, 0.0_Real);
   deepCopy(DefForcing->TracerForcing.SeaIceSaltFluxCell, 0.0_Real);
   deepCopy(DefForcing->TracerForcing.SurfaceTracerFluxCell, 0.0_Real);
   deepCopy(DefForcing->TracerForcing.SensibleHeatFluxCell, 100.0_Real);
   DefForcing->computeAll();

   Array2DReal VertNonLocalFlux = KPPInstance->VertNonLocalFlux;
   parallelFor(
       {Mesh->NCellsSize, VCoord->NVertLayersP1},
       KOKKOS_LAMBDA(I4 ICell, I4 K) { VertNonLocalFlux(ICell, K) = Real(K); });

   TimeInstant Time;
   TimeInterval Interval(1.0, TimeUnits::Seconds);
   AuxState->computeAll(State, TracerArray, 0, 0, Interval);

   const auto ComputeTempTendency = [&]() {
      deepCopy(DefTendencies->TracerTend, 0.0_Real);
      DefTendencies->computeTracerTendenciesOnly(State, AuxState, TracerArray,
                                                 0, 0, Time);
      auto TendencyHost = createHostMirrorCopy(DefTendencies->TracerTend);
      deepCopy(TendencyHost, DefTendencies->TracerTend);
      return TendencyHost(TempIndex, ICellTest, KTest);
   };

   DefTendencies->KPPNonLocalTracerFlux.Enabled = false;
   const Real BaselineTendency                  = ComputeTempTendency();
   DefTendencies->KPPNonLocalTracerFlux.Enabled = true;
   const Real BothDisabledTendency              = ComputeTempTendency();

   if (!isApprox(BothDisabledTendency, BaselineTendency, 1.0e-10_Real,
                 1.0e-12_Real)) {
      Err++;
      LOG_ERROR("TendenciesTest: KPP non-local flux was applied with both "
                "vertical mixing flags disabled");
   }

   DefTendencies->VelVertMixTendencyEnable = true;
   const Real VelocityOnlyTendency         = ComputeTempTendency();
   if (!isApprox(VelocityOnlyTendency, BaselineTendency, 1.0e-10_Real,
                 1.0e-12_Real)) {
      Err++;
      LOG_ERROR("TendenciesTest: KPP non-local flux was applied with only "
                "velocity vertical mixing enabled");
   }

   DefTendencies->VelVertMixTendencyEnable    = false;
   DefTendencies->TracerVertMixTendencyEnable = true;
   const Real EnabledTendency                 = ComputeTempTendency();
   const auto SurfaceFluxHost =
       createHostMirrorCopy(DefForcing->TracerForcing.SurfaceTracerFluxCell);
   deepCopy(SurfaceFluxHost, DefForcing->TracerForcing.SurfaceTracerFluxCell);
   const Real ExpectedTendencyDelta =
       SurfaceFluxHost(TempIndex, ICellTest) * (Real(KTest) - Real(KTest + 1));
   const Real ComputedTendencyDelta = EnabledTendency - BaselineTendency;
   if (isApprox(ExpectedTendencyDelta, 0.0_Real, 0.0_Real, 1.0e-12_Real) ||
       !isApprox(ComputedTendencyDelta, ExpectedTendencyDelta, 1.0e-10_Real,
                 1.0e-12_Real)) {
      Err++;
      LOG_ERROR("TendenciesTest: KPP non-local flux did not run with tracer "
                "vertical mixing enabled (expected {}, got {})",
                ExpectedTendencyDelta, ComputedTendencyDelta);
   }

   DefTendencies->VelVertMixTendencyEnable      = OriginalVelMixEnabled;
   DefTendencies->TracerVertMixTendencyEnable   = OriginalTracerMixEnabled;
   VMix->VelVertMixSetup.Enabled                = OriginalVelSetupEnabled;
   VMix->TracerVertMixSetup.Enabled             = OriginalTracerSetupEnabled;
   KPPInstance->Enabled                         = OriginalKPPEnabled;
   DefTendencies->KPPNonLocalTracerFlux.Enabled = OriginalNonLocalEnabled;
   DefTendencies->SfcTracerForcing.Enabled      = OriginalSfcTracerEnabled;

   if (Err == 0) {
      LOG_INFO("TendenciesTest: KPP non-local vertical-mixing gate PASS");
   }
   return Err;
}

int testSfcThicknessForcing() {
   int Err = 0;

   auto *VCoord        = VertCoord::getDefault();
   auto *DefTendencies = Tendencies::getDefault();
   auto *State         = OceanState::getDefault();
   auto *AuxState      = AuxiliaryState::getDefault();
   auto *DefForcing    = Forcing::getDefault();

   Array3DReal TracerArray = Tracers::getAll(0);

   deepCopy(DefTendencies->PseudoThicknessTend, 0._Real);

   // Set up single test cell at top layer
   const I4 ICellTest = 0;
   const I4 KTop      = VCoord->MinLayerCellH(ICellTest);

   if (KTop > VCoord->MaxLayerCellH(ICellTest)) {
      LOG_ERROR("TendenciesTest: Test cell has no layers for thickness test");
      return -1;
   }

   // Retrieve forcing field views for thickness
   auto &SnowFlux         = DefForcing->TracerForcing.SnowFluxCell;
   auto &RainFlux         = DefForcing->TracerForcing.RainFluxCell;
   auto &EvaporationFlux  = DefForcing->TracerForcing.EvaporationFluxCell;
   auto &SeaIceFreshWater = DefForcing->TracerForcing.SeaIceFreshWaterFluxCell;
   auto &IceRunoffFlux    = DefForcing->TracerForcing.IceRunoffFluxCell;
   auto &RiverRunoffFlux  = DefForcing->TracerForcing.RiverRunoffFluxCell;
   auto &SeaIceSaltFlux   = DefForcing->TracerForcing.SeaIceSaltFluxCell;

   // Initialize all fluxes to zero
   deepCopy(SnowFlux, 0._Real);
   deepCopy(RainFlux, 0._Real);
   deepCopy(EvaporationFlux, 0._Real);
   deepCopy(SeaIceFreshWater, 0._Real);
   deepCopy(IceRunoffFlux, 0._Real);
   deepCopy(RiverRunoffFlux, 0._Real);
   deepCopy(SeaIceSaltFlux, 0._Real);

   // Set test freshwater flux values
   // Rain: 1e-8 kg/m²/s
   const Real TestRain = 1.0e-8_Real;
   // Snow: 5e-9 kg/m²/s
   const Real TestSnow = 5.0e-9_Real;
   // Ice runoff: 2e-9 kg/m²/s
   const Real TestIceRunoff = 2.0e-9_Real;
   // River runoff: 3e-9 kg/m²/s
   const Real TestRiverRunoff = 3.0e-9_Real;
   // Sea ice freshwater: 1e-9 kg/m²/s
   const Real TestSeaIceFreshWater = 1.0e-9_Real;
   // Sea ice salt flux: 1e-4 kg/m²/s (affects thickness via salt)
   const Real TestSeaIceSaltFlux = 1.0e-4_Real;

   OMEGA_SCOPE(LocSnowFlux, SnowFlux);
   OMEGA_SCOPE(LocRainFlux, RainFlux);
   OMEGA_SCOPE(LocIceRunoffFlux, IceRunoffFlux);
   OMEGA_SCOPE(LocRiverRunoffFlux, RiverRunoffFlux);
   OMEGA_SCOPE(LocSeaIceFreshWater, SeaIceFreshWater);
   OMEGA_SCOPE(LocSeaIceSaltFlux, SeaIceSaltFlux);
   Kokkos::parallel_for(
       "SetTestForcingThickness", 1, KOKKOS_LAMBDA(int i) {
          LocRainFlux(ICellTest)         = TestRain;
          LocSnowFlux(ICellTest)         = TestSnow;
          LocIceRunoffFlux(ICellTest)    = TestIceRunoff;
          LocRiverRunoffFlux(ICellTest)  = TestRiverRunoff;
          LocSeaIceFreshWater(ICellTest) = TestSeaIceFreshWater;
          LocSeaIceSaltFlux(ICellTest)   = TestSeaIceSaltFlux;
       });
   DefForcing->computeAll();

   const bool OrigSfcStressEnabled = DefTendencies->SfcStressForcing.Enabled;
   const bool OrigSfcThicknessEnabled =
       DefTendencies->SfcThicknessForcing.Enabled;
   const bool OrigSfcTracerEnabled = DefTendencies->SfcTracerForcing.Enabled;
   const bool OrigPseudoThicknessDiv =
       DefTendencies->PseudoThicknessFluxDiv.Enabled;
   const bool OrigPotentialVortHAdv = DefTendencies->PotentialVortHAdv.Enabled;
   const bool OrigKEGrad            = DefTendencies->KEGrad.Enabled;
   const bool OrigVelocityDiffusion = DefTendencies->VelocityDiffusion.Enabled;
   const bool OrigVelocityHyperDiff = DefTendencies->VelocityHyperDiff.Enabled;
   const bool OrigTracerHorzAdv     = DefTendencies->TracerHorzAdv.Enabled;
   const bool OrigTracerDiffusion   = DefTendencies->TracerDiffusion.Enabled;
   const bool OrigTracerHyperDiff   = DefTendencies->TracerHyperDiff.Enabled;
   const bool OrigSurfaceTracerRestoring =
       DefTendencies->SurfaceTracerRestoring.Enabled;

   DefTendencies->SfcStressForcing.Enabled       = false;
   DefTendencies->SfcThicknessForcing.Enabled    = false;
   DefTendencies->SfcTracerForcing.Enabled       = false;
   DefTendencies->PseudoThicknessFluxDiv.Enabled = false;
   DefTendencies->PotentialVortHAdv.Enabled      = false;
   DefTendencies->KEGrad.Enabled                 = false;
   DefTendencies->VelocityDiffusion.Enabled      = false;
   DefTendencies->VelocityHyperDiff.Enabled      = false;
   DefTendencies->TracerHorzAdv.Enabled          = false;
   DefTendencies->TracerDiffusion.Enabled        = false;
   DefTendencies->TracerHyperDiff.Enabled        = false;
   DefTendencies->SurfaceTracerRestoring.Enabled = false;

   // Compute baseline tendencies (vertical advection is always on)
   int ThickTimeLevel = 0;
   int VelTimeLevel   = 0;
   TimeInstant Time;
   DefTendencies->computePseudoThicknessTendenciesOnly(
       State, AuxState, ThickTimeLevel, VelTimeLevel, Time);

   HostArray2DReal PseudoThicknessTendBaseH =
       createHostMirrorCopy(DefTendencies->PseudoThicknessTend);
   deepCopy(PseudoThicknessTendBaseH, DefTendencies->PseudoThicknessTend);
   const Real BaselineThickTend = PseudoThicknessTendBaseH(ICellTest, KTop);

   // Now enable SfcThicknessForcing and compute again
   DefTendencies->SfcThicknessForcing.Enabled = true;
   DefTendencies->computePseudoThicknessTendenciesOnly(
       State, AuxState, ThickTimeLevel, VelTimeLevel, Time);

   // Calculate expected thickness tendency
   // ThickTend = (Rain + Snow + IceRunoff + RiverRunoff + SeaIceFreshWater +
   // SeaIceSaltFlux) / RhoSw
   const Real ExpectedThickTend =
       (TestRain + TestSnow + TestIceRunoff + TestRiverRunoff +
        TestSeaIceFreshWater + TestSeaIceSaltFlux) /
       RhoSw;

   HostArray2DReal PseudoThicknessTendH =
       createHostMirrorCopy(DefTendencies->PseudoThicknessTend);
   deepCopy(PseudoThicknessTendH, DefTendencies->PseudoThicknessTend);
   const Real ComputedThickTend =
       PseudoThicknessTendH(ICellTest, KTop) - BaselineThickTend;

   constexpr Real RelTol = 1.0e-10_Real;
   constexpr Real AbsTol = 1.0e-12_Real;

   // Check thickness tendency
   if (!isApprox(ComputedThickTend, ExpectedThickTend, RelTol, AbsTol)) {
      Err++;
      LOG_ERROR("TendenciesTest: SfcThicknessForcing thickness tendency FAIL");
      LOG_INFO("  Expected: {},  Computed: {}, Diff: {}", ExpectedThickTend,
               ComputedThickTend,
               Kokkos::abs(ComputedThickTend - ExpectedThickTend));
   } else {
      LOG_INFO("TendenciesTest: SfcThicknessForcing thickness tendency PASS");
   }

   DefTendencies->SfcStressForcing.Enabled       = OrigSfcStressEnabled;
   DefTendencies->SfcThicknessForcing.Enabled    = OrigSfcThicknessEnabled;
   DefTendencies->SfcTracerForcing.Enabled       = OrigSfcTracerEnabled;
   DefTendencies->PseudoThicknessFluxDiv.Enabled = OrigPseudoThicknessDiv;
   DefTendencies->PotentialVortHAdv.Enabled      = OrigPotentialVortHAdv;
   DefTendencies->KEGrad.Enabled                 = OrigKEGrad;
   DefTendencies->VelocityDiffusion.Enabled      = OrigVelocityDiffusion;
   DefTendencies->VelocityHyperDiff.Enabled      = OrigVelocityHyperDiff;
   DefTendencies->TracerHorzAdv.Enabled          = OrigTracerHorzAdv;
   DefTendencies->TracerDiffusion.Enabled        = OrigTracerDiffusion;
   DefTendencies->TracerHyperDiff.Enabled        = OrigTracerHyperDiff;
   DefTendencies->SurfaceTracerRestoring.Enabled = OrigSurfaceTracerRestoring;

   return Err;
}

void finalizeTendenciesTest() {
   Forcing::clear();
   Tracers::clear();
   Frazil::clear();
   PressureGrad::clear();
   KPPMix::destroyInstance();
   VertMix::destroyInstance();
   Eos::destroyInstance();
   AuxiliaryState::clear();
   OceanState::clear();
   VertAdv::clear();
   VertCoord::clear();
   HorzMesh::clear();
   Field::clear();
   Dimension::clear();
   TimeStepper::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

int tendenciesTest(const std::string &MeshFile = "OmegaMesh.nc") {
   int Err = initTendenciesTest(MeshFile);
   if (Err != 0) {
      LOG_CRITICAL("TendenciesTest: Error initializing");
   }

   const auto &Mesh = HorzMesh::getDefault();

   Err += initState();

   Err += testTendencies();

   if (Err == 0) {
      LOG_INFO("TendenciesTest: Successful completion");
   }
   finalizeTendenciesTest();

   return Err;
}

int main(int argc, char *argv[]) {

   int RetVal = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   RetVal += tendenciesTest();

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
