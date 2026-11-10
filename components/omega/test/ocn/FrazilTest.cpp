//===-- Test driver for OMEGA FrazilFormation ---------------------------*- C++
//-*-===//
//
/// \file
/// \brief Minimal test driver for OMEGA frazil formation functor
//
//===-----------------------------------------------------------------------===/

#include "Frazil.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "TimeMgr.h"
#include "VertCoord.h"
#include "mpi.h"

using namespace OMEGA;

constexpr int NVertLayers = 60;

void initFrazilTest(const std::string &mesh) {
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   initLogging(DefEnv);
   LOG_INFO("------ Frazil Unit Tests ------");

   Config("Omega");
   Config::readAll("omega.yml");

   Config *OmegaConfig = Config::getOmegaConfig();
   Config TendConfig("Tendencies");
   OmegaConfig->get(TendConfig);
   TendConfig.set("FrazilTendencyEnable", true);

   Calendar::init("No Leap");
   TimeInstant StartTime(0, 1, 1, 0, 0, 0.0);
   TimeInterval TimeStep(1, TimeUnits::Hours);
   Clock ModelClockTmp(StartTime, TimeStep);
   Clock *ModelClock = &ModelClockTmp;

   IO::init(DefComm);
   Decomp::init(mesh);
   Field::init(ModelClock);
   IOStream::init(ModelClock);
   Halo::init();
   HorzMesh::init(ModelClock);
   VertCoord::init(false);
   Frazil::init();
}

void finalizeFrazilTest() {
   Frazil::clear();
   VertCoord::clear();
   HorzMesh::clear();
   Halo::clear();
   Decomp::clear();
   Field::clear();
   Dimension::clear();
   IOStream::finalize();
   MachEnv::removeAll();
}

// this test only excercises the frazil formation functor (no melt)
// in a cold case: the frazil terms should be
// - strictly positive for ice, liquid, and salt mass
//- strictly negative for ice and liquid energy
// - positive for T tendency and negative for S, H tendencies
void testFrazilFormationCold() {
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   VCoord->NVertLayers = NVertLayers;

   const Real AbsSalinityIn     = 35.0_Real;
   const Real ConservTempIn     = -2.0_Real;
   const Real PressureDbIn      = 100.0_Real;
   const Real PseudoThicknessIn = 10.0_Real;
   const Real RTol              = 1e-10_Real;

   (void)Mesh;

   FrazilFormation ComputeFrazilFormation;
   ComputeFrazilFormation.Phi              = 0.75_Real;
   ComputeFrazilFormation.LayerMassFracMax = 0.1_Real;

   Real AccMIce  = 0.0_Real;
   Real AccMLiq  = 0.0_Real;
   Real AccMSalt = 0.0_Real;
   Real AccELiq  = 0.0_Real;
   Real AccEIce  = 0.0_Real;

   Real HTend = 0.0_Real;
   Real TTend = 0.0_Real;
   Real STend = 0.0_Real;

   ComputeFrazilFormation(AbsSalinityIn, ConservTempIn, PressureDbIn,
                          PseudoThicknessIn, AccMIce, AccMLiq, AccMSalt,
                          AccELiq, AccEIce, HTend, TTend, STend);

   if (AccMIce <= 0.0_Real) {
      ABORT_ERROR(
          "FrazilFormationTestCold: accumulated ice mass is non-positive: {}",
          AccMIce);
   }
   if (AccMLiq <= 0.0_Real) {
      ABORT_ERROR("FrazilFormationTestCold: accumulated liquid mass is "
                  "non-positive: {}",
                  AccMLiq);
   }
   if (AccMSalt <= 0.0_Real) {
      ABORT_ERROR(
          "FrazilFormationTestCold: accumulated salt mass is non-positive: {}",
          AccMSalt);
   }

   if (AccELiq >= 0.0_Real) {
      ABORT_ERROR("FrazilFormationTestCold: accumulated liquid energy is "
                  "positive (exp. negative): {}",
                  AccELiq);
   }
   if (AccEIce >= 0.0_Real) {
      ABORT_ERROR("FrazilFormationTestCold: accumulated ice energy is positive "
                  "(exp. negative): {}",
                  AccEIce);
   }
   if (HTend >= 0.0_Real) {
      ABORT_ERROR(
          "FrazilFormationTestCold: HTend is positive (exp. negative): {}",
          HTend);
   }
   if (TTend <= 0.0_Real) {
      ABORT_ERROR(
          "FrazilFormationTestCold: TTend is negative (exp. positive): {}",
          TTend);
   }
   if (STend >= 0.0_Real) {
      ABORT_ERROR(
          "FrazilFormationTestCold: STend is positive (exp. negative): {}",
          STend);
   }
   LOG_INFO(
       "FrazilFormationTestCold: AccMIce = {}, AccMLiq = {}, AccMSalt = {}, "
       "AccELiq = {}, AccEIce = {}, HTend = {}, TTend = {}, STend = {}",
       AccMIce, AccMLiq, AccMSalt, AccELiq, AccEIce, HTend, TTend, STend);
}

// this test only excercises the frazil formation functor (no melt)
// in a warm case: the frazil FORMATION terms should all be zero
void testFrazilFormationWarm() {
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   VCoord->NVertLayers = NVertLayers;

   const Real AbsSalinityIn     = 35.0_Real;
   const Real ConservTempIn     = 10.0_Real;
   const Real PressureDbIn      = 100.0_Real;
   const Real PseudoThicknessIn = 10.0_Real;
   const Real RTol              = 1e-10_Real;

   (void)Mesh;

   FrazilFormation ComputeFrazilFormation;
   ComputeFrazilFormation.Phi              = 0.75_Real;
   ComputeFrazilFormation.LayerMassFracMax = 0.1_Real;

   Real AccMIce  = 0.0_Real;
   Real AccMLiq  = 0.0_Real;
   Real AccMSalt = 0.0_Real;
   Real AccELiq  = 0.0_Real;
   Real AccEIce  = 0.0_Real;

   Real HTend = 0.0_Real;
   Real TTend = 0.0_Real;
   Real STend = 0.0_Real;

   ComputeFrazilFormation(AbsSalinityIn, ConservTempIn, PressureDbIn,
                          PseudoThicknessIn, AccMIce, AccMLiq, AccMSalt,
                          AccELiq, AccEIce, HTend, TTend, STend);

   if (!isApprox(AccMIce, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFormationTest warm: expected zero AccMIce, got {}",
                  AccMIce);
   }

   if (!isApprox(AccMSalt, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFormationTest warm: expected zero AccMSalt, got {}",
                  AccMSalt);
   }
   if (!isApprox(AccMLiq, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFormationTest warm: expected zero AccMLiq, got {}",
                  AccMLiq);
   }

   if (!isApprox(AccELiq, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFormationTest warm: expected zero AccELiq, got {}",
                  AccELiq);
   }
   if (!isApprox(AccEIce, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFormationTest warm: expected zero AccEIce, got {}",
                  AccEIce);
   }
   if (!isApprox(HTend, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFormationTest warm: expected zero HTend, got {}",
                  HTend);
   }

   if (!isApprox(TTend, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFormationTest warm: expected zero TTend, got {}",
                  TTend);
   }

   if (!isApprox(STend, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFormationTest warm: expected zero STend, got {}",
                  STend);
   }
   LOG_INFO(
       "FrazilFormationTestWarm: AccMIce = {}, AccMLiq = {}, AccMSalt = {}, "
       "AccELiq = {}, AccEIce = {}, HTend = {}, TTend = {}, STend = {}",
       AccMIce, AccMLiq, AccMSalt, AccELiq, AccEIce, HTend, TTend, STend);
}

// this test exercises the frazil formation mass limiter in the TEOS path
void testFrazilFormationMassLimit() {

   const Real AbsSalinityIn     = 20.0_Real;
   const Real ConservTempIn     = -5.0_Real;
   const Real PressureDbIn      = 100.0_Real;
   const Real PseudoThicknessIn = 1.0_Real;
   const Real Phi               = 0.75_Real;
   const Real LayerMassFracMax  = 0.10_Real;
   const Real RTol              = 1e-10_Real;

   FrazilFormation ComputeFrazilFormation;
   ComputeFrazilFormation.Phi              = Phi;
   ComputeFrazilFormation.LayerMassFracMax = LayerMassFracMax;

   Real AccMIce  = 0.0_Real;
   Real AccMLiq  = 0.0_Real;
   Real AccMSalt = 0.0_Real;
   Real AccELiq  = 0.0_Real;
   Real AccEIce  = 0.0_Real;

   Real HTend = 0.0_Real;
   Real TTend = 0.0_Real;
   Real STend = 0.0_Real;

   ComputeFrazilFormation(AbsSalinityIn, ConservTempIn, PressureDbIn,
                          PseudoThicknessIn, AccMIce, AccMLiq, AccMSalt,
                          AccELiq, AccEIce, HTend, TTend, STend);

   const Real ExpectedIceMass =
       PseudoThicknessIn * (1.0_Real - Phi) * LayerMassFracMax;
   const Real ExpectedLiquidMass = Phi / (1.0_Real - Phi) * ExpectedIceMass;
   const Real ExpectedTotalMass  = ExpectedIceMass + ExpectedLiquidMass;

   if (!isApprox(AccMIce, ExpectedIceMass, RTol)) {
      ABORT_ERROR("FrazilFormationMassLimit: expected AccMIce={}, got {}",
                  ExpectedIceMass, AccMIce);
   }
   if (!isApprox(AccMLiq, ExpectedLiquidMass, RTol)) {
      ABORT_ERROR("FrazilFormationMassLimit: expected AccMLiq={}, got {}",
                  ExpectedLiquidMass, AccMLiq);
   }
   if (!isApprox(HTend, -ExpectedTotalMass, RTol)) {
      ABORT_ERROR("FrazilFormationMassLimit: expected HTend={}, got {}",
                  -ExpectedTotalMass, HTend);
   }
   LOG_INFO("FrazilFormationMassLimit: AccMIce = {}, AccMLiq = {}, "
            "HTend = {}",
            AccMIce, AccMLiq, HTend);
}

// this test exercises the Phi parameter in the TEOS frazil formation path
void testFrazilFormationPhi() {
   const Real AbsSalinityIn     = 35.0_Real;
   const Real ConservTempIn     = -2.0_Real;
   const Real PressureDbIn      = 100.0_Real;
   const Real PseudoThicknessIn = 1.0_Real;
   const Real LayerMassFracMax  = 0.10_Real;
   const Real Phi0              = 0.75_Real;
   const Real Phi1              = 0.85_Real;
   const Real RTol              = 1e-10_Real;

   FrazilFormation ComputeFrazilFormation;
   ComputeFrazilFormation.LayerMassFracMax = LayerMassFracMax;

   Real AccMIce0  = 0.0_Real;
   Real AccMLiq0  = 0.0_Real;
   Real AccMSalt0 = 0.0_Real;
   Real AccELiq0  = 0.0_Real;
   Real AccEIce0  = 0.0_Real;
   Real HTend0    = 0.0_Real;
   Real TTend0    = 0.0_Real;
   Real STend0    = 0.0_Real;

   ComputeFrazilFormation.Phi = Phi0;
   ComputeFrazilFormation(AbsSalinityIn, ConservTempIn, PressureDbIn,
                          PseudoThicknessIn, AccMIce0, AccMLiq0, AccMSalt0,
                          AccELiq0, AccEIce0, HTend0, TTend0, STend0);

   Real AccMIce1  = 0.0_Real;
   Real AccMLiq1  = 0.0_Real;
   Real AccMSalt1 = 0.0_Real;
   Real AccELiq1  = 0.0_Real;
   Real AccEIce1  = 0.0_Real;
   Real HTend1    = 0.0_Real;
   Real TTend1    = 0.0_Real;
   Real STend1    = 0.0_Real;

   ComputeFrazilFormation.Phi = Phi1;
   ComputeFrazilFormation(AbsSalinityIn, ConservTempIn, PressureDbIn,
                          PseudoThicknessIn, AccMIce1, AccMLiq1, AccMSalt1,
                          AccELiq1, AccEIce1, HTend1, TTend1, STend1);

   if (!isApprox(AccMIce0, AccMIce1, RTol)) {
      ABORT_ERROR("FrazilFormationPhi: AccMIce changed from {} to {}; "
                  "chosen state may be mass-limit capped",
                  AccMIce0, AccMIce1);
   }

   if ((AccMIce1 + AccMLiq1) <= (AccMIce0 + AccMLiq0)) {
      ABORT_ERROR("FrazilFormationPhi: expected total frazil mass to increase "
                  "from {}, got {}",
                  AccMIce0 + AccMLiq0, AccMIce1 + AccMLiq1);
   }
   if (AccMSalt1 <= AccMSalt0) {
      ABORT_ERROR("FrazilFormationPhi: expected AccMSalt to increase from {}, "
                  "got {}",
                  AccMSalt0, AccMSalt1);
   }
   if ((AccELiq1 + AccEIce1) >= (AccELiq0 + AccEIce0)) {
      ABORT_ERROR("FrazilFormationPhi: expected total frazil energy to become "
                  "more negative than {}, got {}",
                  AccELiq0 + AccEIce0, AccELiq1 + AccEIce1);
   }

   LOG_INFO("FrazilFormationPhi: Phi0 = {}, Phi1 = {}, TotalMass0 = {}, "
            "TotalMass1 = {}, TotalEnergy0 = {}, TotalEnergy1 = {}",
            Phi0, Phi1, AccMIce0 + AccMLiq0, AccMIce1 + AccMLiq1,
            AccELiq0 + AccEIce0, AccELiq1 + AccEIce1);
}

// this test only excercises the fixed-property frazil formation functor (no
// melt) in a warm case: the frazil FORMATION terms should all be zero
void testFixedPropertyFrazilFormationWarm() {

   const Real AbsSalinityIn     = 35.0_Real;
   const Real ConservTempIn     = 10.0_Real;
   const Real PressureDbIn      = 100.0_Real;
   const Real PseudoThicknessIn = 10.0_Real;
   const Real RTol              = 1e-10_Real;
   const Real LayerMassFracMax  = 0.10_Real;

   FixedPropertyFrazilFormation ComputeFrazilFormation;
   ComputeFrazilFormation.LayerMassFracMax = LayerMassFracMax;

   Real AccMIce  = 0.0_Real;
   Real AccMLiq  = 0.0_Real;
   Real AccMSalt = 0.0_Real;
   Real AccELiq  = 0.0_Real;
   Real AccEIce  = 0.0_Real;

   Real HTend = 0.0_Real;
   Real TTend = 0.0_Real;
   Real STend = 0.0_Real;

   Real CtFreezing =
       gsw_ct_freezing_poly(AbsSalinityIn, PressureDbIn, 0.0_Real);

   ComputeFrazilFormation(AbsSalinityIn, ConservTempIn, PressureDbIn,
                          PseudoThicknessIn, AccMIce, AccMSalt, AccEIce, HTend,
                          TTend, STend, CtFreezing);

   if (!isApprox(AccMIce, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFixedPropertyFormationTest warm: expected zero "
                  "AccMIce, got {}",
                  AccMIce);
   }

   if (!isApprox(AccMSalt, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFixedPropertyFormationTest warm: expected zero "
                  "AccMSalt, got {}",
                  AccMSalt);
   }

   if (!isApprox(AccEIce, 0.0_Real, RTol)) {
      ABORT_ERROR("FrazilFixedPropertyFormationTest warm: expected zero "
                  "AccEIce, got {}",
                  AccEIce);
   }
   if (!isApprox(HTend, 0.0_Real, RTol)) {
      ABORT_ERROR(
          "FrazilFixedPropertyFormationTest warm: expected zero HTend, got {}",
          HTend);
   }

   if (!isApprox(TTend, 0.0_Real, RTol)) {
      ABORT_ERROR(
          "FrazilFixedPropertyFormationTest warm: expected zero TTend, got {}",
          TTend);
   }

   if (!isApprox(STend, 0.0_Real, RTol)) {
      ABORT_ERROR(
          "FrazilFixedPropertyFormationTest warm: expected zero STend, got {}",
          STend);
   }
   LOG_INFO("FrazilFixedPropertyFormationTestWarm: AccMIce = {}, AccMLiq = {}, "
            "AccMSalt = {}, "
            "AccELiq = {}, AccEIce = {}, HTend = {}, TTend = {}, STend = {}",
            AccMIce, AccMLiq, AccMSalt, AccELiq, AccEIce, HTend, TTend, STend);
}

// this test only excercises the frazil formation functor (no melt)
// in a cold case: the frazil terms should be
// - strictly positive for ice, liquid, and salt mass
//- strictly negative for ice and liquid energy
// - positive for T tendency and negative for S, H tendencies
void testFixedPropertyFrazilFormationCold() {
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();

   VCoord->NVertLayers = NVertLayers;

   const Real AbsSalinityIn     = 35.0_Real;
   const Real ConservTempIn     = -2.0_Real;
   const Real PressureDbIn      = 100.0_Real;
   const Real PseudoThicknessIn = 10.0_Real;
   const Real RTol              = 1e-10_Real;
   const Real LayerMassFracMax  = 0.10_Real;
   (void)Mesh;

   FixedPropertyFrazilFormation ComputeFrazilFormation;
   ComputeFrazilFormation.LayerMassFracMax = LayerMassFracMax;

   Real AccMIce  = 0.0_Real;
   Real AccMLiq  = 0.0_Real;
   Real AccMSalt = 0.0_Real;
   Real AccELiq  = 0.0_Real;
   Real AccEIce  = 0.0_Real;

   Real HTend = 0.0_Real;
   Real TTend = 0.0_Real;
   Real STend = 0.0_Real;

   Real CtFreezing =
       gsw_ct_freezing_poly(AbsSalinityIn, PressureDbIn, 0.0_Real);

   ComputeFrazilFormation(AbsSalinityIn, ConservTempIn, PressureDbIn,
                          PseudoThicknessIn, AccMIce, AccMSalt, AccEIce, HTend,
                          TTend, STend, CtFreezing);

   if (AccMIce <= 0.0_Real) {
      ABORT_ERROR(
          "FrazilFixedPropertyFormationTestCold: accumulated ice mass is "
          "non-positive: {}",
          AccMIce);
   }

   if (AccMSalt <= 0.0_Real) {
      ABORT_ERROR(
          "FrazilFixedPropertyFormationTestCold: accumulated salt mass is "
          "non-positive: {}",
          AccMSalt);
   }

   if (AccEIce >= 0.0_Real) {
      ABORT_ERROR("FrazilFixedPropertyFormationTestCold: accumulated ice "
                  "energy is positive "
                  "(exp. negative): {}",
                  AccEIce);
   }
   if (HTend >= 0.0_Real) {
      ABORT_ERROR("FrazilFixedPropertyFormationTestCold: HTend is positive "
                  "(exp. negative): {}",
                  HTend);
   }
   if (TTend <= 0.0_Real) {
      ABORT_ERROR("FrazilFixedPropertyFormationTestCold: TTend is negative "
                  "(exp. positive): {}",
                  TTend);
   }
   if (STend >= 0.0_Real) {
      ABORT_ERROR("FrazilFixedPropertyFormationTestCold: STend is positive "
                  "(exp. negative): {}",
                  STend);
   }
   LOG_INFO("FrazilFixedPropertyFormationTestCold: AccMIce = {}, AccMLiq = {}, "
            "AccMSalt = {}, "
            "AccELiq = {}, AccEIce = {}, HTend = {}, TTend = {}, STend = {}",
            AccMIce, AccMLiq, AccMSalt, AccELiq, AccEIce, HTend, TTend, STend);
}

// this test exercises the frazil formation mass limiter in the fixed-property
// path
void testFixedPropertyFrazilFormationMassLimit() {
   const Real AbsSalinityIn     = 20.0_Real;
   const Real ConservTempIn     = -12.0_Real;
   const Real PressureDbIn      = 100.0_Real;
   const Real PseudoThicknessIn = 1.0_Real;
   const Real LayerMassFracMax  = 0.10_Real;
   const Real RTol              = 1e-10_Real;

   FixedPropertyFrazilFormation ComputeFrazilFormation;
   ComputeFrazilFormation.LayerMassFracMax = LayerMassFracMax;

   Real AccMIce  = 0.0_Real;
   Real AccMSalt = 0.0_Real;
   Real AccEIce  = 0.0_Real;

   Real HTend = 0.0_Real;
   Real TTend = 0.0_Real;
   Real STend = 0.0_Real;

   Real CtFreezing =
       gsw_ct_freezing_poly(AbsSalinityIn, PressureDbIn, 0.0_Real);

   ComputeFrazilFormation(AbsSalinityIn, ConservTempIn, PressureDbIn,
                          PseudoThicknessIn, AccMIce, AccMSalt, AccEIce, HTend,
                          TTend, STend, CtFreezing);

   const Real ExpectedIceThickness = PseudoThicknessIn * LayerMassFracMax;

   if (!isApprox(AccMIce, ExpectedIceThickness, RTol)) {
      ABORT_ERROR("FrazilFixedPropertyFormationMassLimit: expected AccMIce={}, "
                  "got {}",
                  ExpectedIceThickness, AccMIce);
   }
   if (!isApprox(HTend, -ExpectedIceThickness, RTol)) {
      ABORT_ERROR("FrazilFixedPropertyFormationMassLimit: expected HTend={}, "
                  "got {}",
                  -ExpectedIceThickness, HTend);
   }
   LOG_INFO("FrazilFixedPropertyFormationMassLimit: AccMIce = {}, HTend = {}",
            AccMIce, HTend);
}

// this test exercises the frazil melt mass limiter in the TEOS path
void testFrazilMeltMassLimit() {
   const Real AbsSalinityIn     = 32.0_Real;
   const Real ConservTempIn     = 35.0_Real;
   const Real PressureDbIn      = 100.0_Real;
   const Real PseudoThicknessIn = 1.0_Real;
   const Real LayerMassFracMax  = 0.10_Real;
   const Real RTol              = 1e-10_Real;

   FrazilMelt ComputeFrazilMelt;
   ComputeFrazilMelt.LayerMassFracMax = LayerMassFracMax;

   const Real AccMIce0  = 0.25_Real;
   const Real AccMLiq0  = 0.75_Real;
   const Real AccMSalt0 = 3.0_Real;
   const Real AccELiq0  = -10000.0_Real;
   const Real AccEIce0  = -100000.0_Real;

   Real AccMIce  = AccMIce0;
   Real AccMLiq  = AccMLiq0;
   Real AccMSalt = AccMSalt0;
   Real AccELiq  = AccELiq0;
   Real AccEIce  = AccEIce0;
   Real HTend    = 0.0_Real;
   Real TTend    = 0.0_Real;
   Real STend    = 0.0_Real;
   Real CtFreezing =
       gsw_ct_freezing_poly(AbsSalinityIn, PressureDbIn, 0.0_Real);

   ComputeFrazilMelt(AbsSalinityIn, ConservTempIn, PressureDbIn,
                     PseudoThicknessIn, AccMIce, AccMLiq, AccMSalt, AccELiq,
                     AccEIce, HTend, TTend, STend, CtFreezing);

   const Real ExpectedFraction =
       PseudoThicknessIn * LayerMassFracMax / (AccMIce0 + AccMLiq0);
   if (ExpectedFraction < 0.0_Real || ExpectedFraction > 1.0_Real) {
      ABORT_ERROR("FrazilMeltMassLimit: expected fraction {} is outside "
                  "[0, 1]",
                  ExpectedFraction);
   }
   if (!isApprox(HTend, ExpectedFraction * (AccMIce0 + AccMLiq0), RTol)) {
      ABORT_ERROR("FrazilMeltMassLimit: expected HTend={}, got {}",
                  ExpectedFraction * (AccMIce0 + AccMLiq0), HTend);
   }
   if (!isApprox(STend, ExpectedFraction * AccMSalt0, RTol)) {
      ABORT_ERROR("FrazilMeltMassLimit: expected STend={}, got {}",
                  ExpectedFraction * AccMSalt0, STend);
   }
   if (!isApprox(TTend, ExpectedFraction * (AccELiq0 + AccEIce0) / Cp0Sw,
                 RTol)) {
      ABORT_ERROR("FrazilMeltMassLimit: expected TTend={}, got {}",
                  ExpectedFraction * (AccELiq0 + AccEIce0) / Cp0Sw, TTend);
   }
   if (!isApprox(AccMIce, (1.0_Real - ExpectedFraction) * AccMIce0, RTol) ||
       !isApprox(AccMLiq, (1.0_Real - ExpectedFraction) * AccMLiq0, RTol) ||
       !isApprox(AccMSalt, (1.0_Real - ExpectedFraction) * AccMSalt0, RTol) ||
       !isApprox(AccELiq, (1.0_Real - ExpectedFraction) * AccELiq0, RTol) ||
       !isApprox(AccEIce, (1.0_Real - ExpectedFraction) * AccEIce0, RTol)) {
      ABORT_ERROR("FrazilMeltMassLimit: remaining reservoirs areinconsistent "
                  "with expected fraction {}",
                  ExpectedFraction);
   }
   LOG_INFO("FrazilMeltMassLimit: fraction = {}, HTend = {}, TTend = {}, "
            "STend = {}",
            ExpectedFraction, HTend, TTend, STend);
}

// this test exercises the frazil melt mass limiter in the fixed-property path
void testFixedPropertyFrazilMeltMassLimit() {
   const Real AbsSalinityIn     = 32.0_Real;
   const Real ConservTempIn     = 35.0_Real;
   const Real PressureDbIn      = 100.0_Real;
   const Real PseudoThicknessIn = 1.0_Real;
   const Real LayerMassFracMax  = 0.10_Real;
   const Real RTol              = 1e-10_Real;

   FixedPropertyFrazilMelt ComputeFrazilMelt;
   ComputeFrazilMelt.LayerMassFracMax = LayerMassFracMax;

   const Real SumIce0    = 0.5_Real;
   const Real SumSalt0   = 2.0_Real;
   const Real SumEnergy0 = -100000.0_Real;

   Real SumIce    = SumIce0;
   Real SumSalt   = SumSalt0;
   Real SumEnergy = SumEnergy0;
   Real HTend     = 0.0_Real;
   Real TTend     = 0.0_Real;
   Real STend     = 0.0_Real;
   Real CtFreezing =
       gsw_ct_freezing_poly(AbsSalinityIn, PressureDbIn, 0.0_Real);

   ComputeFrazilMelt(AbsSalinityIn, ConservTempIn, PressureDbIn,
                     PseudoThicknessIn, SumIce, SumSalt, SumEnergy, HTend,
                     TTend, STend, CtFreezing);

   const Real ExpectedFraction = PseudoThicknessIn * LayerMassFracMax / SumIce0;
   if (ExpectedFraction < 0.0_Real || ExpectedFraction > 1.0_Real) {
      ABORT_ERROR("FrazilFixedPropertyMeltMassLimit: expected fraction {} is "
                  "outside [0, 1]",
                  ExpectedFraction);
   }
   if (!isApprox(HTend, ExpectedFraction * SumIce0, RTol)) {
      ABORT_ERROR("FrazilFixedPropertyMeltMassLimit: expected HTend={}, got {}",
                  ExpectedFraction * SumIce0, HTend);
   }
   if (!isApprox(STend, ExpectedFraction * SumSalt0, RTol)) {
      ABORT_ERROR("FrazilFixedPropertyMeltMassLimit: expected STend={}, got {}",
                  ExpectedFraction * SumSalt0, STend);
   }
   if (!isApprox(TTend, ExpectedFraction * SumEnergy0 / Cp0Sw, RTol)) {
      ABORT_ERROR("FrazilFixedPropertyMeltMassLimit: expected TTend={}, got {}",
                  ExpectedFraction * SumEnergy0 / Cp0Sw, TTend);
   }
   if (!isApprox(SumIce, (1.0_Real - ExpectedFraction) * SumIce0, RTol) ||
       !isApprox(SumSalt, (1.0_Real - ExpectedFraction) * SumSalt0, RTol) ||
       !isApprox(SumEnergy, (1.0_Real - ExpectedFraction) * SumEnergy0, RTol)) {
      ABORT_ERROR("FrazilFixedPropertyMeltMassLimit: remaining reservoir is "
                  "inconsistent with expected fraction {}",
                  ExpectedFraction);
   }
   LOG_INFO("FrazilFixedPropertyMeltMassLimit: fraction = {}, HTend = {}, "
            "TTend = {}, STend = {}",
            ExpectedFraction, HTend, TTend, STend);
}

// this test exercises the frazil formation and melt functors
// in a column of water with both cold and warm layers.
// It turns to frazil column conservation check.
// In dev, there is extra verbose logging in the frazil code (TBRemoved).
void testComputeFrazilColumn() {
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();
   auto *TestFrazil  = Frazil::getDefault();

   if (!TestFrazil) {
      ABORT_ERROR("FrazilTestColumn: default frazil object is null");
   }

   const Real RTol               = 1e-10_Real;
   const Real AbsSalinityCold    = 35.0_Real;
   const Real PressureRef        = 100000.0_Real; // computeFrazil() expects Pa
   const Real PseudoThicknessRef = 10.0_Real;
   const Real ConservTempCold    = -2.0_Real;
   const Real ConservTempWarm    = 0.0_Real;
   const Real ConservTempWarm2   = -1.9_Real;

   Array2DReal AbsSalinity("AbsSalinity", Mesh->NCellsSize, NVertLayers);
   Array2DReal ConservTemp("ConservTemp", Mesh->NCellsSize, NVertLayers);
   Array2DReal Pressure("Pressure", Mesh->NCellsSize, NVertLayers);
   Array2DReal PseudoThickness("PseudoThickness", Mesh->NCellsSize,
                               NVertLayers);

   deepCopy(AbsSalinity, AbsSalinityCold);
   deepCopy(ConservTemp, ConservTempWarm);
   deepCopy(Pressure, PressureRef);
   deepCopy(PseudoThickness, PseudoThicknessRef);

   deepCopy(TestFrazil->AccMIce, 0.0_Real);
   deepCopy(TestFrazil->AccMLiq, 0.0_Real);
   deepCopy(TestFrazil->AccMSalt, 0.0_Real);
   deepCopy(TestFrazil->AccELiq, 0.0_Real);
   deepCopy(TestFrazil->AccEIce, 0.0_Real);
   deepCopy(TestFrazil->FrazilHTend, 0.0_Real);
   deepCopy(TestFrazil->FrazilTTend, 0.0_Real);
   deepCopy(TestFrazil->FrazilSTend, 0.0_Real);

   auto MinLayerCellH = createHostMirrorCopy(VCoord->MinLayerCell);
   auto MaxLayerCellH = createHostMirrorCopy(VCoord->MaxLayerCell);

   const I4 ICell = 0;
   const I4 KMin  = MinLayerCellH(ICell);
   const I4 KMax  = MaxLayerCellH(ICell);
   if ((KMax - KMin + 1) < 4) {
      ABORT_ERROR("FrazilTestColumn: cell {} has fewer than 4 active layers",
                  ICell);
   }

   const I4 KBottom0  = KMax;
   const I4 KBottom1  = KMax - 1;
   const I4 KWarm     = KMax - 2;
   const I4 KTopCold  = KMax - 3;
   const I4 KCold2    = KMin + 3;
   const I4 KCold3    = KMin + 2;
   const I4 KWarm2    = KMin + 1;
   const I4 KTopCold2 = KMin;

   auto ConservTempH               = createHostMirrorCopy(ConservTemp);
   ConservTempH(ICell, KBottom0)   = ConservTempCold;
   ConservTempH(ICell, KBottom1)   = ConservTempCold;
   ConservTempH(ICell, KWarm)      = ConservTempWarm;
   ConservTempH(ICell, KTopCold)   = ConservTempCold;
   ConservTempH(ICell, KCold2 + 2) = ConservTempCold;
   ConservTempH(ICell, KCold2 + 1) = ConservTempCold - .5_Real;
   ConservTempH(ICell, KCold2)     = ConservTempCold;
   ConservTempH(ICell, KCold3)     = ConservTempCold;
   ConservTempH(ICell, KWarm2)     = ConservTempWarm2;
   ConservTempH(ICell, KTopCold2)  = ConservTempCold;
   deepCopy(ConservTemp, ConservTempH);

   const bool SavedConservationCheck = TestFrazil->ConservationCheck;
   TestFrazil->ConservationCheck     = true;

   TestFrazil->computeFrazil(ConservTemp, AbsSalinity, Pressure,
                             PseudoThickness);
   TestFrazil->ConservationCheck = SavedConservationCheck;

   auto HTendH = createHostMirrorCopy(TestFrazil->FrazilHTend);
   auto TTendH = createHostMirrorCopy(TestFrazil->FrazilTTend);
   auto STendH = createHostMirrorCopy(TestFrazil->FrazilSTend);

   if (HTendH(ICell, KBottom0) >= 0.0_Real ||
       TTendH(ICell, KBottom0) <= 0.0_Real ||
       STendH(ICell, KBottom0) >= 0.0_Real) {
      ABORT_ERROR(
          "FrazilTestColumn: bottom cold layer sign check failed (HTend<0, "
          "TTend>0, STend<0 expected)");
   }

   if (HTendH(ICell, KBottom1) >= 0.0_Real ||
       TTendH(ICell, KBottom1) <= 0.0_Real ||
       STendH(ICell, KBottom1) >= 0.0_Real) {
      ABORT_ERROR(
          "FrazilTestColumn: second cold layer sign check failed (HTend<0, "
          "TTend>0, STend<0 expected)");
   }

   if (HTendH(ICell, KWarm) < 0.0_Real || TTendH(ICell, KWarm) > 0.0_Real ||
       STendH(ICell, KWarm) < 0.0_Real) {
      ABORT_ERROR("FrazilTestColumn: warm layer1 sign check failed (HTend>=0, "
                  "TTend<=0, STend>=0 expected)");
   }
   if (HTendH(ICell, KWarm2) < 0.0_Real || TTendH(ICell, KWarm2) > 0.0_Real ||
       STendH(ICell, KWarm) < 0.0_Real) {
      ABORT_ERROR("FrazilTestColumn: warm layer2 sign check failed (HTend>=0, "
                  "TTend<=0, STend>=0 expected)");
   }

   if (HTendH(ICell, KTopCold) >= 0.0_Real ||
       TTendH(ICell, KTopCold) <= 0.0_Real ||
       STendH(ICell, KTopCold) >= 0.0_Real) {
      ABORT_ERROR(
          "FrazilTestColumn: top cold layer sign check failed (HTend<0, "
          "TTend>0, STend<0 expected)");
   }
   if (HTendH(ICell, KTopCold2) >= 0.0_Real ||
       TTendH(ICell, KTopCold2) <= 0.0_Real ||
       STendH(ICell, KTopCold2) >= 0.0_Real) {
      ABORT_ERROR(
          "FrazilTestColumn: top cold layer sign check failed (HTend<0, "
          "TTend>0, STend<0 expected)");
   }

   LOG_INFO(
       "FrazilTestColumn: XTend branch-switching checks passed for ICell={}",
       ICell);
}

// this test exercises the frazil formation and melt functors
// with a depth limit set. Layers deeper than the depth limit
// should have zero frazil tendencies.
void testComputeFrazilDepthLimit() {
   const auto Mesh   = HorzMesh::getDefault();
   const auto VCoord = VertCoord::getDefault();
   auto *TestFrazil  = Frazil::getDefault();

   if (!TestFrazil) {
      ABORT_ERROR("FrazilTestColumn: default frazil object is null");
   }

   const Real RTol               = 1e-12_Real;
   const Real AbsSalinityCold    = 35.0_Real;
   const Real PressureRef        = 100000.0_Real; // computeFrazil() expects Pa
   const Real PseudoThicknessRef = 10.0_Real;
   const Real ConservTempCold    = -2.0_Real;
   const Real ConservTempWarm    = 0.0_Real;
   const Real ConservTempWarm2   = -1.9_Real;

   Array2DReal AbsSalinity("AbsSalinity", Mesh->NCellsSize, NVertLayers);
   Array2DReal ConservTemp("ConservTemp", Mesh->NCellsSize, NVertLayers);
   Array2DReal Pressure("Pressure", Mesh->NCellsSize, NVertLayers);
   Array2DReal PseudoThickness("PseudoThickness", Mesh->NCellsSize,
                               NVertLayers);

   deepCopy(AbsSalinity, AbsSalinityCold);
   deepCopy(ConservTemp, ConservTempWarm);
   deepCopy(Pressure, PressureRef);
   deepCopy(PseudoThickness, PseudoThicknessRef);

   deepCopy(TestFrazil->AccMIce, 0.0_Real);
   deepCopy(TestFrazil->AccMLiq, 0.0_Real);
   deepCopy(TestFrazil->AccMSalt, 0.0_Real);
   deepCopy(TestFrazil->AccELiq, 0.0_Real);
   deepCopy(TestFrazil->AccEIce, 0.0_Real);
   deepCopy(TestFrazil->FrazilHTend, 0.0_Real);
   deepCopy(TestFrazil->FrazilTTend, 0.0_Real);
   deepCopy(TestFrazil->FrazilSTend, 0.0_Real);

   auto MinLayerCellH = createHostMirrorCopy(VCoord->MinLayerCell);
   auto MaxLayerCellH = createHostMirrorCopy(VCoord->MaxLayerCell);

   const I4 ICell = 0;
   const I4 KMin  = MinLayerCellH(ICell);
   const I4 KMax  = MaxLayerCellH(ICell);
   if ((KMax - KMin + 1) < 10) {
      ABORT_ERROR("FrazilTestColumn: cell {} has fewer than 10 active layers",
                  ICell);
   }

   const I4 KBottom0  = KMax;
   const I4 KBottom1  = KMax - 1;
   const I4 KWarm     = KMax - 2;
   const I4 KTopCold  = KMax - 3;
   const I4 KCold2    = KMin + 3;
   const I4 KCold3    = KMin + 2;
   const I4 KWarm2    = KMin + 1;
   const I4 KTopCold2 = KMin;

   auto ConservTempH               = createHostMirrorCopy(ConservTemp);
   ConservTempH(ICell, KBottom0)   = ConservTempCold;
   ConservTempH(ICell, KBottom1)   = ConservTempCold;
   ConservTempH(ICell, KWarm)      = ConservTempWarm;
   ConservTempH(ICell, KTopCold)   = ConservTempCold;
   ConservTempH(ICell, KCold2 + 2) = ConservTempCold;
   ConservTempH(ICell, KCold2 + 1) = ConservTempCold - .5_Real;
   ConservTempH(ICell, KCold2)     = ConservTempCold;
   ConservTempH(ICell, KCold3)     = ConservTempCold;
   ConservTempH(ICell, KWarm2)     = ConservTempWarm2;
   ConservTempH(ICell, KTopCold2)  = ConservTempCold;
   deepCopy(ConservTemp, ConservTempH);

   const bool SavedConservationCheck = TestFrazil->ConservationCheck;
   const Real SavedDepthLimit        = TestFrazil->DepthLimit;
   const Real TestDepthLimit         = 35.0_Real; // this needs to be positive
   // if TestDepthLimit is negative, test will fail:
   // - the code assume depthlimit < 0 mean no limit (i.e. full depth frazil)
   // - the test below will exclude all layers and fail because Tend !=0.

   // Populate GeomZMid explicitly for the test column.
   auto GeomZMidH = createHostMirrorCopy(VCoord->GeomZMid);
   for (I4 K = KMin; K <= KMax; ++K) {
      GeomZMidH(ICell, K) = -10.0_Real * (K - KMin + 1);
   }
   deepCopy(VCoord->GeomZMid, GeomZMidH);

   TestFrazil->ConservationCheck = true;
   TestFrazil->DepthLimit        = TestDepthLimit;
   TestFrazil->computeFrazil(ConservTemp, AbsSalinity, Pressure,
                             PseudoThickness);
   TestFrazil->ConservationCheck = SavedConservationCheck;
   TestFrazil->DepthLimit        = SavedDepthLimit;

   auto HTendH = createHostMirrorCopy(TestFrazil->FrazilHTend);
   auto TTendH = createHostMirrorCopy(TestFrazil->FrazilTTend);
   auto STendH = createHostMirrorCopy(TestFrazil->FrazilSTend);

   bool FoundExcludedLayer = false;
   for (I4 K = KMin; K <= KMax; ++K) {
      const Real Depth    = GeomZMidH(ICell, K);
      const Real AbsDepth = Depth < 0.0_Real ? -Depth : Depth;

      if (AbsDepth > TestDepthLimit) {
         FoundExcludedLayer = true;
         if (!isApprox(HTendH(ICell, K), 0.0_Real, RTol) ||
             !isApprox(TTendH(ICell, K), 0.0_Real, RTol) ||
             !isApprox(STendH(ICell, K), 0.0_Real, RTol)) {
            ABORT_ERROR("FrazilDepthLimitTest: excluded layer K={} has "
                        "non-zero tendencies (HTend={}, TTend={}, STend={})",
                        K, HTendH(ICell, K), TTendH(ICell, K),
                        STendH(ICell, K));
         }
      }
   }
   if (!FoundExcludedLayer) {
      ABORT_ERROR("FrazilDepthLimitTest: no layers were excluded for ICell={} "
                  "with DepthLimit={}",
                  ICell, TestDepthLimit);
   }

   LOG_INFO("FrazilDepthLimitTest: DepthLimit={} exclusion check passed for "
            "ICell={}",
            TestDepthLimit, ICell);
}

void frazilTest(const std::string &MeshFile = "OmegaMesh.nc") {
   initFrazilTest(MeshFile);
   testFrazilFormationCold();
   testFrazilFormationWarm();
   testFrazilFormationMassLimit();
   testFrazilFormationPhi();
   testFrazilMeltMassLimit();
   testFixedPropertyFrazilFormationCold();
   testFixedPropertyFrazilFormationWarm();
   testFixedPropertyFrazilFormationMassLimit();
   testFixedPropertyFrazilMeltMassLimit();
   testComputeFrazilColumn();
   testComputeFrazilDepthLimit();
   finalizeFrazilTest();
}

int main(int argc, char *argv[]) {
   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   frazilTest();

   LOG_INFO("------ Frazil Unit Tests Successful ------");

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   return 0;
}
