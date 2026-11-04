//===-- ocn/Tendencies.cpp - Tendencies ------------------*- C++ -*-===//
//
// The Tendencies class is responsible for managing tendencies of state
// variables and tracers. It contains arrays that store the tendency data
// and provides methods for computing different tendency groups for use
// within the timestepping algorithm. At initialization, it determines which
// tendency terms are enabled.
//
//===----------------------------------------------------------------------===//
#include <iostream>

#include "CustomTendencyTerms.h"
#include "Eos.h"
#include "Error.h"
#include "Field.h"
#include "FillValues.h"
#include "Forcing.h"
#include "HorzOperators.h"
#include "KPPMix.h"
#include "OceanState.h"
#include "PGrad.h"
#include "Pacer.h"
#include "Tendencies.h"
#include "TimeStepper.h"
#include "Tracers.h"
#include "VertAdv.h"
#include "VertMix.h"
#include <string>

namespace OMEGA {

Tendencies *Tendencies::DefaultTendencies = nullptr;
std::map<std::string, std::unique_ptr<Tendencies>> Tendencies::AllTendencies;

//------------------------------------------------------------------------------
// Initialize the tendencies. Assumes that HorzMesh, VertCoord, VertAdv, and
// TimeStepper  has already been initialized.
void Tendencies::init() {
   Error Err; // error code

   HorzMesh *DefHorzMesh = HorzMesh::getDefault();
   OMEGA_REQUIRE(DefHorzMesh,
                 "Null default HorzMesh pointer in Tendencies::init");
   VertCoord *DefVertCoord = VertCoord::getDefault();
   OMEGA_REQUIRE(DefVertCoord,
                 "Null default VertCoord pointer in Tendencies::init");
   VertAdv *DefVertAdv = VertAdv::getDefault();
   OMEGA_REQUIRE(DefVertAdv,
                 "Null default VertAdv pointer in Tendencies::init");
   TimeStepper *DefTimeStepper = TimeStepper::getDefault();
   OMEGA_REQUIRE(DefTimeStepper,
                 "Null default TimeStepper pointer in Tendencies::init");
   Eos *DefEos = Eos::getInstance();
   OMEGA_REQUIRE(DefEos, "Null default Eos pointer in Tendencies::init");
   PressureGrad *DefPGrad = PressureGrad::getDefault();
   OMEGA_REQUIRE(DefPGrad,
                 "Null default PressureGrad pointer in Tendencies::init");
   VertMix *DefVertMix = VertMix::getInstance();
   OMEGA_REQUIRE(DefVertMix,
                 "Null default VertMix pointer in Tendencies::init");

   I4 NTracers = Tracers::getNumTracers();

   // Get TendConfig group
   Config *OmegaConfig = Config::getOmegaConfig();
   OMEGA_REQUIRE(OmegaConfig, "Null OmegaConfig pointer in Tendencies::init");
   Config TendConfig("Tendencies");
   Err += OmegaConfig->get(TendConfig);
   CHECK_ERROR_ABORT(Err, "Tendencies: Tendencies group not found in Config");

   // Check if use the customized tendencies. If it is not found in the
   // config, we assume it is not used (false)
   bool UseCustomTendency = false;
   Err += TendConfig.get("UseCustomTendency", UseCustomTendency);

   /// Instances of custom tendencies - empty by default
   CustomTendencyType CustomThickTend;
   CustomTendencyType CustomVelTend;

   if (UseCustomTendency) {
      // Check if use manufactured tendency terms if it is not found in
      // the config file, we will assume it is not used (false)
      bool ManufacturedTend = false;
      Error ManufacturedTendErr =
          TendConfig.get("ManufacturedSolutionTendency", ManufacturedTend);

      if (ManufacturedTend) {
         ManufacturedSolution ManufacturedSol;
         ManufacturedSol.init();

         CustomThickTend = ManufacturedSol.ManufacturedThickTend;
         CustomVelTend   = ManufacturedSol.ManufacturedVelTend;

      } // if ManufacturedTend

   } // end if UseCustomTendency

   TimeInterval TimeStep = DefTimeStepper->getTimeStep();

   // Ceate default tendencies
   Tendencies::DefaultTendencies =
       create("Default", DefHorzMesh, DefVertCoord, DefVertAdv, DefPGrad,
              DefEos, DefVertMix, NTracers, TimeStep, &TendConfig,
              CustomThickTend, CustomVelTend);

   DefaultTendencies->readConfig(OmegaConfig);

   if (DefaultTendencies->TracerHorzAdv.Enabled)
      DefaultTendencies->TracerHorzAdv.init();
} // end init

//------------------------------------------------------------------------------
// Destroys the tendencies
Tendencies::~Tendencies() {

   // No operations needed, Kokkos arrays removed when no longer in scope

} // end destructor

//------------------------------------------------------------------------------
// Removes all tendencies instances before exit
void Tendencies::clear() {
   AllTendencies.clear();
   DefaultTendencies = nullptr; // prevent dangling pointer
} // end clear

//------------------------------------------------------------------------------
// Removes tendencies from list by name
void Tendencies::erase(const std::string &Name) {

   AllTendencies.erase(Name);

} // end erase

//------------------------------------------------------------------------------
// Get default tendencies
Tendencies *Tendencies::getDefault() {

   return Tendencies::DefaultTendencies;

} // end get default

//------------------------------------------------------------------------------
// Get tendencies by name
Tendencies *Tendencies::get(const std::string &Name ///< [in] Name of tendencies
) {

   auto it = AllTendencies.find(Name);

   if (it != AllTendencies.end()) {
      return it->second.get();
   } else {
      LOG_ERROR(
          "Tendencies::get: Attempt to retrieve non-existent tendencies:");
      LOG_ERROR("{} has not been defined or has been removed", Name);
      return nullptr;
   }

} // end get tendencies

//------------------------------------------------------------------------------
// Check if tendencies exists by name
bool Tendencies::exists(const std::string &Name ///< [in] Name of tendencies
) {

   return AllTendencies.find(Name) != AllTendencies.end();

} // end exists

//------------------------------------------------------------------------------
// read and set config options
void Tendencies::readConfig(Config *OmegaConfig ///< [in] Omega config
) {
   Error Err; // error code

   Config TendConfig("Tendencies");
   Err += OmegaConfig->get(TendConfig);
   CHECK_ERROR_ABORT(Err, "Tendencies: Tendencies group not found in Config");

   Err += TendConfig.get("ThicknessFluxTendencyEnable",
                         this->PseudoThicknessFluxDiv.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: ThicknessFluxTendencyEnable not found in TendConfig");

   Err += TendConfig.get("PVTendencyEnable", this->PotentialVortHAdv.Enabled);
   CHECK_ERROR_ABORT(Err,
                     "Tendencies: PVTendencyEnable not found in TendConfig");

   Err += TendConfig.get("KETendencyEnable", this->KEGrad.Enabled);
   CHECK_ERROR_ABORT(Err,
                     "Tendencies: KETendencyEnable not found in TendConfig");

   Err += TendConfig.get("SSHTendencyEnable", this->SSHGrad.Enabled);
   CHECK_ERROR_ABORT(Err,
                     "Tendencies: SSHTendencyEnable not found in TendConfig");

   Err +=
       TendConfig.get("VelDiffTendencyEnable", this->VelocityDiffusion.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: VelDiffTendencyEnable not found in TendConfig");

   Err += TendConfig.get("VelHyperDiffTendencyEnable",
                         this->VelocityHyperDiff.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: VelHyperDiffTendencyEnable not found in TendConfig");

   if (this->VelocityDiffusion.Enabled) {
      Err += TendConfig.get("ViscDel2", this->VelocityDiffusion.ViscDel2);
      CHECK_ERROR_ABORT(Err, "Tendencies: ViscDel2 not found in TendConfig");
   }

   if (this->VelocityHyperDiff.Enabled) {
      Err += TendConfig.get("ViscDel4", this->VelocityHyperDiff.ViscDel4);
      CHECK_ERROR_ABORT(Err, "Tendencies: ViscDel4 not found in TendConfig");
      Err += TendConfig.get("DivFactor", this->VelocityHyperDiff.DivFactor);
      CHECK_ERROR_ABORT(Err, "Tendencies: DivFactor not found in TendConfig");
   }
   Err += TendConfig.get("TracerHorzAdvTendencyEnable",
                         this->TracerHorzAdv.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: TracerHorzAdvTendencyEnable not found in TendConfig");
   if (this->TracerHorzAdv.Enabled) {
      Config AdvectConfig("Advection");
      Err += OmegaConfig->get(AdvectConfig);
      CHECK_ERROR_ABORT(Err, "Tendencies: Advection group not in Config");

      I4 Order = 0;
      Err += AdvectConfig.get("HorzTracerFluxOrder", Order);
      CHECK_ERROR_ABORT(
          Err, "Tendencies: HorzTracerFluxOrder not found in AdvectConfig");
      OMEGA_REQUIRE(Order >= 2 && Order <= 4,
                    "HorzTracerFluxOrder: Only values are 2, 3, 4, found {}",
                    Order);

      if (Order == 2) {
         this->TracerHorzAdv.ForceLowOrder = true;
         this->TracerHorzAdv.Coef3rdOrder  = 0;
      }
      if (Order == 3) {
         Err += AdvectConfig.get("Coef3rdOrder", TracerHorzAdv.Coef3rdOrder);
         CHECK_ERROR_ABORT(
             Err, "Tendencies: Coef3rdOrder not found in AdvectConfig");
      }
      if (Order == 4) {
         this->TracerHorzAdv.Coef3rdOrder = 0;
      }
      Err += AdvectConfig.get("HorzTracerFluxLimiterEnable", TracerHorzAdv.FCT);
      OMEGA_REQUIRE(
          !TracerHorzAdv.FCT || 2 < Order,
          "HorzTracerFluxOrder: Since HorzTracerFluxLimiterEnable is true, "
          "HorzTracerFluxOrder must be greater than 2. Found Order={}",
          Order);
      CHECK_ERROR_ABORT(
          Err,
          "Tendencies: HorzTracerFluxLimiterEnable not found in AdvectConfig");
      if (TracerHorzAdv.FCT) {
         Err += AdvectConfig.get("HorzTracerFluxLimiterBudgetsEnable",
                                 TracerHorzAdv.ComputeBudgets);
         CHECK_ERROR_ABORT(Err,
                           "Tendencies: HorzTracerFluxLimiterBudgetsEnable not "
                           "found in AdvectConfig");
         Err += AdvectConfig.get("HorzTracerFluxLimiterMonotonicityCheckEnable",
                                 TracerHorzAdv.MonotonicityCheck);
         CHECK_ERROR_ABORT(
             Err, "Tendencies: HorzTracerFluxLimiterMonotonicityCheckEnable "
                  "not found in AdvectConfig");
      }
   }
   Err += TendConfig.get("TracerDiffTendencyEnable",
                         this->TracerDiffusion.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: TracerDiffTendencyEnable not found in TendConfig");

   Err += TendConfig.get("SfcStressForcingTendencyEnable",
                         this->SfcStressForcing.Enabled);
   CHECK_ERROR_ABORT(
       Err,
       "Tendencies: SfcStressForcingTendencyEnable not found in TendConfig");

   Config BottomDragConfig("BottomDragTendency");
   Err += TendConfig.get(BottomDragConfig);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: BottomDragTendency group not found in TendConfig");

   bool BottomDragTendencyEnabled = false;
   Err += BottomDragConfig.get("Enable", BottomDragTendencyEnabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: BottomDragTendency Enable not found in TendConfig");

   std::string BottomDragMode;
   Err += BottomDragConfig.get("Mode", BottomDragMode);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: BottomDragTendency Mode not found in TendConfig");

   std::string BottomDragType;
   Err += BottomDragConfig.get("Type", BottomDragType);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: BottomDragTendency Type not found in TendConfig");

   Real BottomDragCoeff = 0.0_Real;
   Err += BottomDragConfig.get("BottomDragCoeff", BottomDragCoeff);
   CHECK_ERROR_ABORT(Err, "Tendencies: BottomDragTendency BottomDragCoeff not "
                          "found in TendConfig");

   if (BottomDragType != "Constant" && BottomDragType != "constant") {
      ABORT_ERROR("Tendencies: BottomDragTendency Type should be one of "
                  " 'Constant' but got {}:",
                  BottomDragType);
   }

   this->ExplicitBottomDrag.Coeff              = BottomDragCoeff;
   this->VMix->VelVertMixSetup.BottomDragCoeff = BottomDragCoeff;

   if (BottomDragTendencyEnabled) {
      if (BottomDragMode == "Explicit" || BottomDragMode == "explicit") {
         this->ExplicitBottomDrag.Enabled = true;
      } else if (BottomDragMode == "Implicit" || BottomDragMode == "implicit") {
         this->VMix->VelVertMixSetup.ImplicitBottomDragEnabled = true;
      } else {
         ABORT_ERROR("Tendencies: BottomDragTendency Mode must be Explicit or "
                     "Implicit");
      }
   }

   Err += TendConfig.get("SfcThicknessForcingTendencyEnable",
                         this->SfcThicknessForcing.Enabled);
   CHECK_ERROR_ABORT(
       Err,
       "Tendencies: SfcThicknessForcingTendencyEnable not found in TendConfig");

   Err += TendConfig.get("SfcTracerForcingTendencyEnable",
                         this->SfcTracerForcing.Enabled);
   CHECK_ERROR_ABORT(
       Err,
       "Tendencies: SfcTracerForcingTendencyEnable not found in TendConfig");

   Config PenSwConfig("PenetratingShortwaveTendency");
   Err += TendConfig.get(PenSwConfig);
   CHECK_ERROR_ABORT(Err, "Tendencies: PenetratingShortwaveTendency group not "
                          "found in TendConfig");

   Err += PenSwConfig.get("Enable", this->PenetratingShortwave.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: PenetratingShortwaveTendency Enable not found in "
            "PenetratingShortwaveTendency config");
   Err += PenSwConfig.get("NearIrFraction",
                          this->PenetratingShortwave.NearIrFraction);
   CHECK_ERROR_ABORT(Err, "Tendencies: NearIrFraction not found in "
                          "PenetratingShortwaveTendency config");
   Err +=
       PenSwConfig.get("NearIrCoeff", this->PenetratingShortwave.NearIrCoeff);
   CHECK_ERROR_ABORT(Err, "Tendencies: NearIrCoeff not found in "
                          "PenetratingShortwaveTendency config");
   Err +=
       PenSwConfig.get("RedFraction", this->PenetratingShortwave.RedFraction);
   CHECK_ERROR_ABORT(Err, "Tendencies: RedFraction not found in "
                          "PenetratingShortwaveTendency config");
   Err +=
       PenSwConfig.get("BlueFraction", this->PenetratingShortwave.BlueFraction);
   CHECK_ERROR_ABORT(Err, "Tendencies: BlueFraction not found in "
                          "PenetratingShortwaveTendency config");

   const Real TotalFraction = this->PenetratingShortwave.NearIrFraction +
                              this->PenetratingShortwave.RedFraction +
                              this->PenetratingShortwave.BlueFraction;
   if (Kokkos::abs(TotalFraction - 1.0_Real) > 1.0e-5_Real) {
      ABORT_ERROR(
          "Tendencies: The sum of NearIrFraction ({}), RedFraction ({}), "
          "and BlueFraction ({}) must equal 1 (got {})",
          this->PenetratingShortwave.NearIrFraction,
          this->PenetratingShortwave.RedFraction,
          this->PenetratingShortwave.BlueFraction, TotalFraction);
   }

   if (this->PenetratingShortwave.Enabled && !this->SfcTracerForcing.Enabled) {
      ABORT_ERROR(
          "Tendencies: PenetratingShortwaveTendency is enabled, but "
          "SfcTracerForcingTendencyEnable is false. Penetrating shortwave "
          "requires surface tracer forcing to be enabled.");
   }

   this->SfcTracerForcing.IncludeShortWaveHeatFlux =
       !this->PenetratingShortwave.Enabled;

   if (this->TracerDiffusion.Enabled) {
      Err += TendConfig.get("EddyDiff2", this->TracerDiffusion.EddyDiff2);
      CHECK_ERROR_ABORT(Err, "Tendencies: EddyDiff2 not found in TendConfig");
   }

   Err += TendConfig.get("TracerHyperDiffTendencyEnable",
                         this->TracerHyperDiff.Enabled);
   CHECK_ERROR_ABORT(
       Err,
       "Tendencies: TracerHyperDiffTendencyEnable not found in TendConfig");

   if (this->TracerHyperDiff.Enabled) {
      Err += TendConfig.get("EddyDiff4", this->TracerHyperDiff.EddyDiff4);
      CHECK_ERROR_ABORT(Err, "Tendencies: EddyDiff4 not found in TendConfig");
   }

   Err += TendConfig.get("PressureGradTendencyEnable", this->PGrad->Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: PressureGradTendencyEnable not found in TendConfig");

   Err += TendConfig.get("FrazilTendencyEnable", this->FrazilTerm.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: FrazilTendencyEnable not found in TendConfig");

   Err += TendConfig.get("SurfaceTracerRestoringEnable",
                         this->SurfaceTracerRestoring.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: SurfaceTracerRestoringEnable not found in TendConfig");
   if (this->SurfaceTracerRestoring.Enabled) {
      Config SurfRestConfig("SurfaceRestoring");
      Err += OmegaConfig->get(SurfRestConfig);
      Err += SurfRestConfig.get("PistonVelocity",
                                this->SurfaceTracerRestoring.PistonVelocity);
      CHECK_ERROR_ABORT(
          Err,
          "Tendencies: PistonVelocity not found in SurfaceRestoringConfig");

      std::vector<std::string> TracersToRestore;
      SurfRestConfig.get("TracersToRestore", TracersToRestore);

      // Enable restoring for specified individual tracers
      I4 NumInvalidTracers = 0;
      std::vector<I4> TracerIdsToRestoreVec;
      for (const auto &TracerName : TracersToRestore) {
         I4 TracerIndex = -1;
         Tracers::getIndex(TracerIndex, TracerName);
         if (TracerIndex == -1) {
            LOG_ERROR("Tendencies: Tracer {} in TracersToRestore is not "
                      "defined",
                      TracerName);
            NumInvalidTracers++;
         } else {
            bool IsDuplicate = false;
            for (const auto ExistingTracerIndex : TracerIdsToRestoreVec) {
               if (ExistingTracerIndex == TracerIndex) {
                  IsDuplicate = true;
                  break;
               }
            }
            if (!IsDuplicate) {
               TracerIdsToRestoreVec.push_back(TracerIndex);
            }
         }
      }
      if (NumInvalidTracers > 0) {
         ABORT_ERROR("Tendencies: {} invalid tracer(s) in TracersToRestore "
                     "configuration",
                     NumInvalidTracers);
      }

      this->SurfaceTracerRestoring.NTracersToRestore =
          static_cast<I4>(TracerIdsToRestoreVec.size());
      this->SurfaceTracerRestoring.TracerIdsToRestore = Array1DI4(
          "TracerIdsToRestore", this->SurfaceTracerRestoring.NTracersToRestore);
      deepCopy(this->SurfaceTracerRestoring.TracerIdsToRestore,
               HostArray1DI4(TracerIdsToRestoreVec.data(),
                             TracerIdsToRestoreVec.size()));
   }

   // Validate VertMix tendency
   Err += TendConfig.get("VelVertMixTendencyEnable",
                         this->VelVertMixTendencyEnable);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: VelVertMixTendencyEnable not found in TendConfig");
   this->VMix->VelVertMixSetup.Enabled = this->VelVertMixTendencyEnable;

   if (this->VMix->VelVertMixSetup.ImplicitBottomDragEnabled &&
       !this->VMix->VelVertMixSetup.Enabled) {
      ABORT_ERROR("Tendencies: BottomDragTendency Mode Implicit requires "
                  "VelVertMixTendencyEnable to be true");
   }

   if (this->VMix->VelVertMixSetup.RayleighDampingEnabled &&
       !this->VMix->VelVertMixSetup.Enabled) {
      ABORT_ERROR("Tendencies: RayleighDamping Enable requires "
                  "VelVertMixTendencyEnable to be true");
   }

   Err += TendConfig.get("TracerVertMixTendencyEnable",
                         this->TracerVertMixTendencyEnable);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: TracerVertMixTendencyEnable not found in TendConfig");
   this->VMix->TracerVertMixSetup.Enabled = this->TracerVertMixTendencyEnable;

   if (this->VMix->VelVertMixSetup.Enabled ||
       this->VMix->TracerVertMixSetup.Enabled) {

      if (!this->EqState) {
         ABORT_ERROR("Tendencies: Eos must be initialized when"
                     "vertical mixing tendencies are enabled");
      }
      if (!this->VMix) {
         ABORT_ERROR("Tendencies: VertMix must be initialized when"
                     "vertical mixing tendencies are enabled");
      }
      // Optional KPP non-local tracer tendency: no abort if missing
      Error KPPNonLocalTracerErr =
          TendConfig.get("KPPNonLocalTracerFluxTendencyEnable",
                         this->KPPNonLocalTracerFlux.Enabled);
      if (!KPPNonLocalTracerErr.isSuccess()) {
         KPPNonLocalTracerErr.reset();
         this->KPPNonLocalTracerFlux.Enabled = false;
      }

      Error KPPNonLocalTracerDiagErr =
          TendConfig.get("KPPNonLocalTracerDiagnosticsEnable",
                         this->KPPNonLocalTracerDiagnosticsEnable);
      if (!KPPNonLocalTracerDiagErr.isSuccess()) {
         KPPNonLocalTracerDiagErr.reset();
         this->KPPNonLocalTracerDiagnosticsEnable = true;
      }
   }
}

//------------------------------------------------------------------------------
// Configure the velocity tendency for a mode-split time stepper
void Tendencies::setModeSplit(CoriolisTendMode Mode, ///< [in] Coriolis mode
                              Real SplitFactorIn ///< [in] split-explicit factor
) {

   if (SplitFactorIn != 0._Real && this->SSHGrad.Enabled) {
      ABORT_ERROR("Tendencies: SSHTendencyEnable must be false for a "
                  "mode-split time stepper, which uses the barotropic "
                  "pressure anomaly gradient instead");
   }

   this->CoriolisMode = Mode;
   this->SplitFactor  = SplitFactorIn;
}

//------------------------------------------------------------------------------
// Define fields associated with tendencies
void Tendencies::defineFields() {
   std::string PseudoThicknessTendFieldName = "PseudoThicknessTend";
   std::string NormalVelocityTendFieldName  = "NormalVelocityTend";
   std::string TracerTendFieldName          = "TracerTend";
   std::string KPPNonLocalTracerTempDiagFieldName =
       "KPPNonLocalTracerTempTendDiag";
   std::string KPPNonLocalTracerTempColSumFieldName =
       "KPPNonLocalTracerTempColumnSumDiag";
   if (Name != "Default") {
      PseudoThicknessTendFieldName.append(Name);
      NormalVelocityTendFieldName.append(Name);
      TracerTendFieldName.append(Name);
      KPPNonLocalTracerTempDiagFieldName.append(Name);
      KPPNonLocalTracerTempColSumFieldName.append(Name);
   }

   int NDims = 2;
   std::vector<std::string> DimNamesThickness(NDims);
   DimNamesThickness[0]          = "NCells";
   DimNamesThickness[1]          = "NVertLayers";
   auto PseudoThicknessTendField = Field::create(
       PseudoThicknessTendFieldName, "Pseudo-thickness tendency", "m s-1", "",
       -9.99E+10, 9.99E+10, NDims, DimNamesThickness);
   NDims = 3;
   std::vector<std::string> DimNamesTracer(NDims);
   DimNamesTracer[0] = "NTracers";
   DimNamesTracer[1] = "NCells";
   DimNamesTracer[2] = "NVertLayers";
   auto TracerTendField =
       Field::create(TracerTendFieldName, "Tracer tendency", "", "", -9.99E+10,
                     9.99E+10, NDims, DimNamesTracer);
   NDims = 2;
   std::vector<std::string> DimNamesVelocity(NDims);
   DimNamesVelocity[0] = "NEdges";
   DimNamesVelocity[1] = "NVertLayers";
   auto NormalVelocityTendField =
       Field::create(NormalVelocityTendFieldName, "Normal velocity tendency",
                     "m s-2", "", -9.99E+10, 9.99E+10, NDims, DimNamesVelocity);

   NDims = 2;
   std::vector<std::string> DimNamesTempDiag(NDims);
   DimNamesTempDiag[0] = "NCells";
   DimNamesTempDiag[1] = "NVertLayers";
   auto KPPNonLocalTracerTempDiagField =
       Field::create(KPPNonLocalTracerTempDiagFieldName,
                     "Temperature non-local KPP tendency diagnostic", "K m s-1",
                     "", -9.99E+10, 9.99E+10, NDims, DimNamesTempDiag);

   NDims = 1;
   std::vector<std::string> DimNamesCellOnly(NDims);
   DimNamesCellOnly[0]                   = "NCells";
   auto KPPNonLocalTracerTempColSumField = Field::create(
       KPPNonLocalTracerTempColSumFieldName,
       "Temperature non-local tendency column-sum diagnostic", "K m s-1", "",
       -9.99E+10, 9.99E+10, NDims, DimNamesCellOnly);

   std::string TendGroupName = "Tendencies";
   if (Name != "Default") {
      TendGroupName.append(Name);
   }
   auto TendGroup = FieldGroup::create(TendGroupName);

   TendGroup->addField(PseudoThicknessTendFieldName);
   TendGroup->addField(NormalVelocityTendFieldName);
   TendGroup->addField(TracerTendFieldName);
   TendGroup->addField(KPPNonLocalTracerTempDiagFieldName);
   TendGroup->addField(KPPNonLocalTracerTempColSumFieldName);

   PseudoThicknessTendField->attachData<Array2DReal>(PseudoThicknessTend);
   NormalVelocityTendField->attachData<Array2DReal>(NormalVelocityTend);
   TracerTendField->attachData<Array3DReal>(TracerTend);
   KPPNonLocalTracerTempDiagField->attachData<Array2DReal>(
       KPPNonLocalTracerTempTendDiag, false);
   KPPNonLocalTracerTempColSumField->attachData<Array1DReal>(
       KPPNonLocalTracerTempColumnSumDiag, false);
   // Attaching with FillOnAttach=false keeps the zero initialization but does
   // not record the _FillValue that IOStream needs to write these fields
   KPPNonLocalTracerTempDiagField->addMetadata("_FillValue", FillValueReal);
   KPPNonLocalTracerTempColSumField->addMetadata("_FillValue", FillValueReal);

} // end defineFields

//------------------------------------------------------------------------------
// Construct a new group of tendencies
Tendencies::Tendencies(const std::string &Name_, ///< [in] Name for tendencies
                       const HorzMesh *Mesh,     ///< [in] Horizontal mesh
                       VertCoord *VCoord,        ///< [in] Vertical coordinate
                       VertAdv *VAdv_,           ///< [in] Vertical advection
                       PressureGrad *PGrad,      ///< [in] Pressure gradient
                       Eos *EqState,             ///< [in] Equation of state
                       VertMix *VMix,            ///< [in] Vertical mixing
                       int NTracersIn,           ///< [in] Number of tracers
                       TimeInterval TimeStepIn,  ///< [in] Time step
                       Config *Options,          ///< [in] Configuration options
                       CustomTendencyType InCustomThicknessTend,
                       CustomTendencyType InCustomVelocityTend)
    : Mesh(Mesh), VCoord(VCoord), VAdv(VAdv_),
      PseudoThicknessFluxDiv(Mesh, VCoord), PotentialVortHAdv(Mesh, VCoord),
      KEGrad(Mesh, VCoord), SSHGrad(Mesh, VCoord),
      CoriolisAcceleration(Mesh, VCoord), VelocityDiffusion(Mesh, VCoord),
      VelocityHyperDiff(Mesh, VCoord), SfcStressForcing(Mesh, VCoord),
      ExplicitBottomDrag(Mesh, VCoord), SfcThicknessForcing(Mesh, VCoord),
      SfcTracerForcing(Mesh, VCoord, Tracers::IndxTemp, Tracers::IndxSalt,
                       EqState),
      PenetratingShortwave(Mesh, VCoord, Tracers::IndxTemp),
      TracerHorzAdv(Mesh, VCoord, VAdv_), TracerDiffusion(Mesh, VCoord),
      KPPNonLocalTracerFlux(Mesh, VCoord), TracerHyperDiff(Mesh, VCoord),
      SurfaceTracerRestoring(Mesh), FrazilTerm(Mesh, VCoord),
      CustomThicknessTend(InCustomThicknessTend),
      CustomVelocityTend(InCustomVelocityTend), EqState(EqState), PGrad(PGrad),
      VMix(VMix) {

   if (this->VMix) {
      VelVertMixTendencyEnable    = this->VMix->VelVertMixSetup.Enabled;
      TracerVertMixTendencyEnable = this->VMix->TracerVertMixSetup.Enabled;
   }

   // Tendency arrays
   PseudoThicknessTend = Array2DReal("PseudoThicknessTend", Mesh->NCellsSize,
                                     VCoord->NVertLayers);
   NormalVelocityTend =
       Array2DReal("NormalVelocityTend", Mesh->NEdgesSize, VCoord->NVertLayers);
   TracerTend = Array3DReal("TracerTend", NTracersIn, Mesh->NCellsSize,
                            VCoord->NVertLayers);
   // Provisional tracers are only needed for flux-corrected vertical advection
   if (VAdv->VertAdvChoice == VertAdvOption::FCT) {
      ProvTracers = Array3DReal("ProvTracers", NTracersIn, Mesh->NCellsSize,
                                VCoord->NVertLayers);
   }
   KPPNonLocalTracerTempTendDiag = Array2DReal(
       "KPPNonLocalTracerTempTendDiag", Mesh->NCellsSize, VCoord->NVertLayers);
   KPPNonLocalTracerTempColumnSumDiag =
       Array1DReal("KPPNonLocalTracerTempColumnSumDiag", Mesh->NCellsSize);
   deepCopy(KPPNonLocalTracerTempTendDiag, 0.0_Real);
   deepCopy(KPPNonLocalTracerTempColumnSumDiag, 0.0_Real);

   Name = Name_;

   NTracers = NTracersIn;
   TimeStep = TimeStepIn;

   defineFields();

} // end constructor

// Create a non-default group of tendencies
Tendencies *Tendencies::create(const std::string &Name,
                               const HorzMesh *Mesh, ///< [in] Horizontal mesh
                               VertCoord *VCoord, ///< [in] Vertical coordinate
                               VertAdv *VAdv,     ///< [in] Vertical advection
                               PressureGrad *PGrad, ///< [in] Pressure gradient
                               Eos *EqState,        ///< [in] Equation of state
                               VertMix *VMix,       ///< [in] Vertical mixing
                               int NTracersIn,      ///< [in] Number of tracers
                               TimeInterval TimeStep, ///< [in] Time step
                               Config *Options, ///< [in] Configuration options
                               CustomTendencyType CustomThicknessTend,
                               CustomTendencyType CustomVelocityTend) {

   OMEGA_REQUIRE(Mesh,
                 "Null HorzMesh pointer in Tendencies::create with Name = {}",
                 Name);
   OMEGA_REQUIRE(VCoord,
                 "Null VCoord pointer in Tendencies::create with Name = {}",
                 Name);
   OMEGA_REQUIRE(
       VAdv, "Null VertAdv pointer in Tendencies::create with Name = {}", Name);
   OMEGA_REQUIRE(
       PGrad, "Null PressureGrad pointer in Tendencies::create with Name = {}",
       Name);
   OMEGA_REQUIRE(EqState,
                 "Null Eos pointer in Tendencies::create with Name = {}", Name);
   OMEGA_REQUIRE(
       VMix, "Null VertMix pointer in Tendencies::create with Name = {}", Name);
   OMEGA_REQUIRE(Options,
                 "Null Config pointer in Tendencies::create with Name = {}",
                 Name);

   // Check to see if tendencies of the same name already exist and
   // if so, exit with an error
   if (AllTendencies.find(Name) != AllTendencies.end()) {
      LOG_ERROR("Attempted to create Tendencies with name {} but Tendencies of "
                "that name already exists",
                Name);
      return nullptr;
   }

   // create new tendencies on the heap and put it in a map of
   // unique_ptrs, which will manage its lifetime
   auto *NewTendencies = new Tendencies(
       Name, Mesh, VCoord, VAdv, PGrad, EqState, VMix, NTracersIn, TimeStep,
       Options, CustomThicknessTend, CustomVelocityTend);

   AllTendencies.emplace(Name, NewTendencies);

   return get(Name);
}

//------------------------------------------------------------------------------
// Accumulate f times tangential velocity reconstruction for edge-centered 2D
// fields
void Tendencies::computeCoriolisAccelerationOnEdge(
    const Array2DReal &Tend,          ///< [inout] velocity tendency
    const Array2DReal &NormalVelEdge, ///< [in] normal velocity on edges
    const Array1DReal &FEdge          ///< [in] Coriolis parameter on edges
) const {

   computeCoriolisAccelerationOnEdge(Tend, Tend, NormalVelEdge, FEdge);
}

//------------------------------------------------------------------------------
// Write a base tendency plus f times tangential velocity reconstruction into a
// separate output array, for iterative Coriolis solves
void Tendencies::computeCoriolisAccelerationOnEdge(
    const Array2DReal &Tend,          ///< [out] base plus Coriolis tendency
    const Array2DReal &BaseTend,      ///< [in] tendency without Coriolis
    const Array2DReal &NormalVelEdge, ///< [in] normal velocity on edges
    const Array1DReal &FEdge          ///< [in] Coriolis parameter on edges
) const {

   // Coriolis acceleration should be turned off if the PV tendency is disabled,
   // in which case the output is just the base tendency.
   if (!PotentialVortHAdv.Enabled) {
      deepCopy(Tend, BaseTend);
      return;
   }

   OMEGA_SCOPE(LocCoriolisAcceleration, CoriolisAcceleration);

   Pacer::start("Tend:coriolisAccelerationOnEdge2DBase", 2);
   parallelForOuter(
       LaunchConfig({Mesh->NEdgesAll}, TeamScratch<Real>(VCoord->NVertLayers)),
       KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          LocCoriolisAcceleration(Team, Tend, BaseTend, IEdge, NormalVelEdge,
                                  FEdge);
       });
   Pacer::stop("Tend:coriolisAccelerationOnEdge2DBase", 2);
}

//------------------------------------------------------------------------------
// Accumulate f times tangential velocity reconstruction for edge-centered 1D
// fields
void Tendencies::computeCoriolisAccelerationOnEdge(
    const Array1DReal &Tend,          ///< [inout] barotropic velocity tendency
    const Array1DReal &NormalVelEdge, ///< [in] normal velocity on edges
    const Array1DReal &FEdge          ///< [in] Coriolis parameter on edges
) const {

   // Shares PVTendencyEnable with the 2D overload above
   if (!PotentialVortHAdv.Enabled)
      return;

   OMEGA_SCOPE(LocCoriolisAcceleration, CoriolisAcceleration);

   Pacer::start("Tend:coriolisAccelerationOnEdge1D", 2);
   parallelFor(
       {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge) {
          LocCoriolisAcceleration(Tend, IEdge, NormalVelEdge, FEdge);
       });
   Pacer::stop("Tend:coriolisAccelerationOnEdge1D", 2);
}

//------------------------------------------------------------------------------
// Compute tendencies for the pseudo-thickness equation
void Tendencies::computePseudoThicknessTendenciesOnly(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {

   Array2DReal NormalVelEdge = State->getNormalVelocity(VelTimeLevel);

   OMEGA_SCOPE(LocPseudoThicknessTend, PseudoThicknessTend);
   OMEGA_SCOPE(LocThicknessFluxDiv, PseudoThicknessFluxDiv);
   OMEGA_SCOPE(LocSfcThicknessForcing, SfcThicknessForcing);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   Pacer::start("Tend:computePseudoThicknessTendenciesOnly", 1);

   parallelForOuter(
       {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          const int KMin = MinLayerCell(ICell);
          const int KMax = MaxLayerCell(ICell);

          parallelForInner(
              Team, Range{KMin, KMax},
              INNER_LAMBDA(int K) { LocPseudoThicknessTend(ICell, K) = 0; });
       });

   // Compute pseudo-thickness flux divergence
   const Array2DReal &ThickFluxEdge =
       AuxState->PseudoThicknessAux.FluxPseudoThickEdge;
   const auto &NormalTransportVelocity =
       AuxState->TransportAux.NormalTransportVelocity;

   if (LocThicknessFluxDiv.Enabled) {
      Pacer::start("Tend:thicknessFluxDiv", 2);
      parallelForOuter(
          LaunchConfig({Mesh->NCellsAll},
                       TeamScratch<Real>(VCoord->NVertLayers)),
          KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
             LocThicknessFluxDiv(Team, LocPseudoThicknessTend, ICell,
                                 ThickFluxEdge, NormalTransportVelocity);
          });
      Pacer::stop("Tend:thicknessFluxDiv", 2);
   }

   Pacer::start("Tend:computePseudoThicknessVAdvTend", 2);
   // Compute thickness tendency from vertical advection
   VAdv->computePseudoThicknessVAdvTend(PseudoThicknessTend);
   Pacer::stop("Tend:computePseudoThicknessVAdvTend", 2);

   if (LocSfcThicknessForcing.Enabled) {
      Pacer::start("Tend:sfcThicknessForcing", 2);
      const auto *ForcingState = Forcing::getDefault();

      const auto &SnowFlux = ForcingState->TracerForcing.SnowFluxCell;
      const auto &RainFlux = ForcingState->TracerForcing.RainFluxCell;
      const auto &EvaporationFlux =
          ForcingState->TracerForcing.EvaporationFluxCell;
      const auto &SeaIceFreshWaterFlux =
          ForcingState->TracerForcing.SeaIceFreshWaterFluxCell;
      const auto &IceRunoffFlux = ForcingState->TracerForcing.IceRunoffFluxCell;
      const auto &RiverRunoffFlux =
          ForcingState->TracerForcing.RiverRunoffFluxCell;
      const auto &SeaIceSaltFlux =
          ForcingState->TracerForcing.SeaIceSaltFluxCell;

      parallelFor(
          {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell) {
             LocSfcThicknessForcing(LocPseudoThicknessTend, ICell, SnowFlux,
                                    RainFlux, EvaporationFlux,
                                    SeaIceFreshWaterFlux, IceRunoffFlux,
                                    RiverRunoffFlux, SeaIceSaltFlux);
          });
      Pacer::stop("Tend:sfcThicknessForcing", 2);
   }

   if (CustomThicknessTend) {
      Pacer::start("Tend:customThicknessTend", 2);
      CustomThicknessTend(LocPseudoThicknessTend, State, AuxState,
                          ThickTimeLevel, VelTimeLevel, Time);
      Pacer::stop("Tend:customThicknessTend", 2);
   }

   Pacer::stop("Tend:computePseudoThicknessTendenciesOnly", 1);

} // end thickness tendency compute

//------------------------------------------------------------------------------
// Compute tendencies for normal velocity equation
void Tendencies::computeVelocityTendenciesOnly(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    const Array3DReal &TracerArray, ///< [in] Tracer array
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    int TracerTimeLevel,            ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {

   OMEGA_SCOPE(LocNormalVelocityTend, NormalVelocityTend);
   OMEGA_SCOPE(LocPotentialVortHAdv, PotentialVortHAdv);
   OMEGA_SCOPE(LocKEGrad, KEGrad);
   OMEGA_SCOPE(LocSSHGrad, SSHGrad);
   OMEGA_SCOPE(LocVelocityDiffusion, VelocityDiffusion);
   OMEGA_SCOPE(LocVelocityHyperDiff, VelocityHyperDiff);
   OMEGA_SCOPE(LocSfcStressForcing, SfcStressForcing);
   OMEGA_SCOPE(LocExplicitBottomDrag, ExplicitBottomDrag);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);
   OMEGA_SCOPE(LocSshCell, VCoord->SshCell);

   Pacer::start("Tend:computeVelocityTendenciesOnly", 1);

   // Zero NVTend over all layers where at least one neighboring cell is active:
   // [MinLayerEdgeTop, MaxLayerEdgeBot]. This includes boundary edges adjacent
   // to land or bathymetry steps, which have no computed tendency but should
   // show 0, not FillValueReal. Fully inactive layers (outside this range)
   // keep their FillValueReal from attachData().
   VCoord->zeroEdgeField(NormalVelocityTend, Mesh->NEdgesAll);

   // Compute vorticity horizontal advection.
   // With CoriolisTendMode::Separate only the relative vorticity is advected
   // here; the time stepper adds the linear Coriolis acceleration itself so it
   // can be iterated.
   const Array2DReal &FluxPseudoThickEdge =
       AuxState->PseudoThicknessAux.FluxPseudoThickEdge;
   const Array2DReal &NormRVortEdge = AuxState->VorticityAux.NormRelVortEdge;
   const Array2DReal &NormFEdge     = AuxState->VorticityAux.NormPlanetVortEdge;
   Array2DReal NormalVelEdge        = State->getNormalVelocity(VelTimeLevel);
   if (LocPotentialVortHAdv.Enabled) {
      Pacer::start("Tend:PotentialVortHAdv", 2);
      if (CoriolisMode == CoriolisTendMode::Separate) {
         parallelForOuter(
             LaunchConfig({Mesh->NEdgesAll},
                          TeamScratch<Real>(VCoord->NVertLayers)),
             KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
                LocPotentialVortHAdv(Team, LocNormalVelocityTend, IEdge,
                                     NormRVortEdge, FluxPseudoThickEdge,
                                     NormalVelEdge);
             });
      } else {
         parallelForOuter(
             LaunchConfig({Mesh->NEdgesAll},
                          TeamScratch<Real>(VCoord->NVertLayers)),
             KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
                LocPotentialVortHAdv(Team, LocNormalVelocityTend, IEdge,
                                     NormRVortEdge, NormFEdge,
                                     FluxPseudoThickEdge, NormalVelEdge);
             });
      }
      Pacer::stop("Tend:PotentialVortHAdv", 2);
   }

   // Compute kinetic energy gradient
   const Array2DReal &KECell = AuxState->KineticAux.KineticEnergyCell;
   if (LocKEGrad.Enabled) {
      Pacer::start("Tend:KEGrad", 2);
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             LocKEGrad(Team, LocNormalVelocityTend, IEdge, KECell);
          });
      Pacer::stop("Tend:KEGrad", 2);
   }

   // Compute the barotropic pressure gradient for the mode-split time steppers
   const Array1DReal &SSHCell = LocSshCell;
   if (SplitFactor != 0._Real) {
      Pacer::start("Tend:barotropicPressureGradTerm", 2);
      const Array1DReal &BtrPressAnomaly =
          State->getBarotropicPressureAnomaly(VelTimeLevel);
      const Array1DReal &DepthMeanSpecVol = EqState->DepthMeanSpecificVolume;
      PGrad->computeBarotropicPressureGrad(LocNormalVelocityTend,
                                           BtrPressAnomaly, DepthMeanSpecVol);
      Pacer::stop("Tend:barotropicPressureGradTerm", 2);
   } else if (LocSSHGrad.Enabled) {
      Pacer::start("Tend:SSHGrad", 2);
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             LocSSHGrad(Team, LocNormalVelocityTend, IEdge, SSHCell);
          });
      Pacer::stop("Tend:SSHGrad", 2);
   }

   // Compute del2 horizontal diffusion
   const Array2DReal &DivCell     = AuxState->KineticAux.VelocityDivCell;
   const Array2DReal &RVortVertex = AuxState->VorticityAux.RelVortVertex;
   if (LocVelocityDiffusion.Enabled) {
      Pacer::start("Tend:velocityDiffusion", 2);
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             LocVelocityDiffusion(Team, LocNormalVelocityTend, IEdge, DivCell,
                                  RVortVertex);
          });
      Pacer::stop("Tend:velocityDiffusion", 2);
   }

   // Compute del4 horizontal diffusion
   const Array2DReal &Del2DivCell = AuxState->VelocityDel2Aux.Del2DivCell;
   const Array2DReal &Del2RVortVertex =
       AuxState->VelocityDel2Aux.Del2RelVortVertex;
   if (LocVelocityHyperDiff.Enabled) {
      Pacer::start("Tend:velocityHyperDiff", 2);
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             LocVelocityHyperDiff(Team, LocNormalVelocityTend, IEdge,
                                  Del2DivCell, Del2RVortVertex);
          });
      Pacer::stop("Tend:velocityHyperDiff", 2);
   }

   Pacer::start("Tend:computeVelocityVAdvTend", 2);
   // Compute velocity tendency from vertical advection
   VAdv->computeVelocityVAdvTend(NormalVelocityTend, NormalVelEdge,
                                 FluxPseudoThickEdge);
   Pacer::stop("Tend:computeVelocityVAdvTend", 2);

   // Compute surface stress forcing
   const auto *ForcingState = Forcing::getDefault();
   const auto &NormalStressEdge =
       ForcingState->SfcStressForcing.NormalStressEdge;
   const auto &MeanPseudoThickEdge =
       AuxState->PseudoThicknessAux.MeanPseudoThickEdge;
   if (LocSfcStressForcing.Enabled) {
      Pacer::start("Tend:sfcStressForcing", 2);
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             LocSfcStressForcing(Team, LocNormalVelocityTend, IEdge,
                                 NormalStressEdge, MeanPseudoThickEdge);
          });
      Pacer::stop("Tend:sfcStressForcing", 2);
   }

   // Compute explicit bottom drag
   if (LocExplicitBottomDrag.Enabled) {
      Pacer::start("Tend:explicitBottomDrag", 2);
      parallelFor(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge) {
             LocExplicitBottomDrag(LocNormalVelocityTend, IEdge, NormalVelEdge,
                                   KECell, MeanPseudoThickEdge);
          });
      Pacer::stop("Tend:explicitBottomDrag", 2);
   }

   if (CustomVelocityTend) {
      Pacer::start("Tend:customVelocityTend", 2);
      CustomVelocityTend(LocNormalVelocityTend, State, AuxState, ThickTimeLevel,
                         VelTimeLevel, Time);
      Pacer::stop("Tend:customVelocityTend", 2);
   }

   // Compute pressure gradient
   if (PGrad->Enabled) {

      Pacer::start("Tend:pressureGradTerm", 2);
      Array2DReal PseudoThick       = State->getPseudoThickness(ThickTimeLevel);
      const auto &PressureMid       = VCoord->PressureMid;
      const auto &PressureInterface = VCoord->PressureInterface;
      const auto &SpecVol           = EqState->SpecVol;
      const auto &GeomZInterface    = VCoord->GeomZInterface;

      // The FiniteVolume scheme reconstructs the layer-mean temperature and
      // salinity in pressure and uses the specific volume derivatives Eos
      // holds; the Centered scheme ignores both
      I4 ConservTempIdx;
      I4 AbsSalinityIdx;
      Tracers::getIndex(ConservTempIdx, "Temperature");
      Tracers::getIndex(AbsSalinityIdx, "Salinity");
      const Array2DReal ConservTemp = Kokkos::subview(
          TracerArray, ConservTempIdx, Kokkos::ALL, Kokkos::ALL);
      const Array2DReal AbsSalinity = Kokkos::subview(
          TracerArray, AbsSalinityIdx, Kokkos::ALL, Kokkos::ALL);

      PGrad->computePressureGrad(
          LocNormalVelocityTend, PressureMid, PressureInterface, SpecVol,
          GeomZInterface, PseudoThick, ConservTemp, AbsSalinity, EqState);
      Pacer::stop("Tend:pressureGradTerm", 2);
   }

   Pacer::stop("Tend:computeVelocityTendenciesOnly", 1);

} // end velocity tendency compute

void Tendencies::computeTracerTendenciesOnly(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    const Array3DReal &TracerArray, ///< [in] Tracer array
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time,               ///< [in] Time
    Real FinalUpdateWeight          ///< [in] final update weight
) {
   OMEGA_SCOPE(LocTracerArray, TracerArray);
   OMEGA_SCOPE(LocTracerTend, TracerTend);
   OMEGA_SCOPE(LocTracerHorzAdv, TracerHorzAdv);
   OMEGA_SCOPE(LocTracerDiffusion, TracerDiffusion);
   OMEGA_SCOPE(LocKPPNonLocalTracerFlux, KPPNonLocalTracerFlux);
   OMEGA_SCOPE(LocTracerHyperDiff, TracerHyperDiff);
   OMEGA_SCOPE(LocSurfaceTracerRestoring, SurfaceTracerRestoring);
   OMEGA_SCOPE(LocSfcTracerForcing, SfcTracerForcing);
   OMEGA_SCOPE(LocMinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(LocMinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(LocMaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   Array2DReal PseudoThickCell = State->getPseudoThickness(ThickTimeLevel);
   Array2DReal FluxPseudoThickEdge =
       AuxState->PseudoThicknessAux.FluxPseudoThickEdge;

   OMEGA_SCOPE(LocPseudoThickCell, PseudoThickCell);
   OMEGA_SCOPE(LocNormalTransportVelocity,
               AuxState->TransportAux.NormalTransportVelocity);
   OMEGA_SCOPE(LocFluxPseudoThickEdge, FluxPseudoThickEdge);
   OMEGA_SCOPE(LocKPPNonLocalTracerTempTendDiag, KPPNonLocalTracerTempTendDiag);
   OMEGA_SCOPE(LocKPPNonLocalTracerTempColumnSumDiag,
               KPPNonLocalTracerTempColumnSumDiag);
   const bool LocKPPNonLocalTracerDiagnosticsEnable =
       KPPNonLocalTracerDiagnosticsEnable;
   I4 TempTracerIndex = -1;
   const bool LocHasTempTracer =
       (Tracers::getIndex(TempTracerIndex, "Temperature") == 0);
   const I4 LocTempTracerIndex = TempTracerIndex;

   Pacer::start("Tend:computeTracerTendenciesOnly", 1);

   parallelForOuter(
       "Tend:LocTracerTend.init", {NTracers, Mesh->NCellsAll},
       KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
          const int KMin = LocMinLayerCell(ICell);
          const int KMax = LocMaxLayerCell(ICell);
          parallelForInner(
              Team, Range{KMin, KMax},
              INNER_LAMBDA(int K) { LocTracerTend(L, ICell, K) = 0; });
       });

   if (LocKPPNonLocalTracerDiagnosticsEnable) {
      parallelForOuter(
          {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
             const int KMin = LocMinLayerCell(ICell);
             const int KMax = LocMaxLayerCell(ICell);
             parallelForInner(
                 Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                    LocKPPNonLocalTracerTempTendDiag(ICell, K) = 0.0_Real;
                 });
             Kokkos::single(Kokkos::PerTeam(Team), [&]() {
                LocKPPNonLocalTracerTempColumnSumDiag(ICell) = 0.0_Real;
             });
          });
   }

   // compute tracer horizotal advection
   R8 Dt = 0;
   TimeStep.get(Dt, TimeUnits::Seconds);
   if (LocTracerHorzAdv.Enabled) {
      Pacer::start("Tend:tracerHorzAdv", 2);
      if (LocTracerHorzAdv.FCT) {
         parallelForOuter(
             "Tend:FCTProvisionaLayerThicknesses", {Mesh->NCellsHaloH(0)},
             KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
                LocTracerHorzAdv.FCTProvisionaLayerThicknesses(
                    Team, ICell, Dt, LocFluxPseudoThickEdge, LocPseudoThickCell,
                    LocNormalTransportVelocity);
             });
         for (int L = 0; L < NTracers; ++L) {
            parallelForOuter(
                "Tend:FCTTracerCurFill", {Mesh->NCellsAll},
                KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
                   LocTracerHorzAdv.FCTTracerCurFill(Team, L, ICell,
                                                     LocTracerArray);
                });
            parallelForOuter(
                "Tend:FCTTracerMinMax", {Mesh->NCellsAll},
                KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
                   LocTracerHorzAdv.FCTTracerMinMax(Team, ICell);
                });
            parallelForOuter(
                "Tend:FCTHighAndLowOrderFlux", {Mesh->NEdgesHaloH(1)},
                KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
                   LocTracerHorzAdv.FCTHighAndLowOrderFlux(
                       Team, IEdge, LocFluxPseudoThickEdge,
                       LocNormalTransportVelocity);
                });
            parallelForOuter(
                "Tend:FCTFluxInOut", {Mesh->NCellsHaloH(0)},
                KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
                   LocTracerHorzAdv.FCTFluxInOut(Team, ICell, Dt,
                                                 LocPseudoThickCell);
                });
            parallelForOuter(
                "Tend:FCTRescaleHighOrderFlux", {Mesh->NEdgesHaloH(0)},
                KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
                   LocTracerHorzAdv.FCTRescaleHighOrderFlux(Team, IEdge);
                });
            parallelForOuter(
                "Tend:FCTAccumulateHighOrderFlux", {Mesh->NCellsOwned},
                KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
                   const Array2DReal Tend = Kokkos::subview(
                       LocTracerTend, L, Kokkos::ALL, Kokkos::ALL);
                   LocTracerHorzAdv.FCTAccumulateHighOrderFlux(
                       Team, ICell, Dt, Tend, LocPseudoThickCell);
                });
            if (LocTracerHorzAdv.ComputeBudgets) {
               parallelForOuter(
                   "Tend:FCTComputeBudgetAdvectionEdgeFlux",
                   {Mesh->NEdgesHaloH(1)},
                   KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
                      LocTracerHorzAdv.FCTComputeBudgetAdvectionEdgeFlux(
                          Team, L, IEdge);
                   });
               parallelForOuter(
                   "Tend:FCTComputeBudgetAdvectionTendency",
                   {Mesh->NCellsOwned},
                   KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
                      LocTracerHorzAdv.FCTComputeBudgetAdvectionTendency(
                          Team, L, ICell);
                   });
            }
            if (LocTracerHorzAdv.MonotonicityCheck) {
               const I4 NCellsOwned = Mesh->NCellsOwned;
               parallelForOuter(
                   "Tend:FCTMonotonicityCheck", {Mesh->NCellsOwned},
                   KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
                      LocTracerHorzAdv.FCTMonotonicityCheck(Team, ICell);
                   });
            }
         }
      } else {
         const Array2DReal &FluxPseudoThickEdge =
             AuxState->PseudoThicknessAux.FluxPseudoThickEdge;

         parallelForOuter(
             LaunchConfig({NTracers, Mesh->NEdgesAll},
                          TeamScratch<Real>(VCoord->NVertLayers)),
             KOKKOS_LAMBDA(int L, int IEdge, const TeamMember &Team) {
                LocTracerHorzAdv(Team, L, IEdge, TracerArray,
                                 FluxPseudoThickEdge,
                                 LocNormalTransportVelocity);
             });
         parallelForOuter(
             LaunchConfig({NTracers, Mesh->NCellsAll},
                          TeamScratch<Real>(VCoord->NVertLayers)),
             KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
                LocTracerHorzAdv(Team, LocTracerTend, L, ICell);
             });
      }
      Pacer::stop("Tend:tracerHorzAdv", 2);
   }

   // Flux-corrected vertical advection starts from the tracers after
   // horizontal advection alone, as in MPAS-Ocean, so that its limiter bounds
   // the combined horizontal and vertical update. TracerTend holds only the
   // horizontal advection tendency at this point.
   if (VAdv->TracerVertAdvEnabled &&
       VAdv->VertAdvChoice == VertAdvOption::FCT) {
      OMEGA_SCOPE(LocProvTracers, ProvTracers);
      OMEGA_SCOPE(LocProvPseudoThick,
                  AuxState->PseudoThicknessAux.ProvPseudoThickness);
      parallelForOuter(
          "Tend:ProvTracers", {NTracers, Mesh->NCellsOwned},
          KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
             const int KMin = LocMinLayerCell(ICell);
             const int KMax = LocMaxLayerCell(ICell);
             parallelForInner(
                 Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                    LocProvTracers(L, ICell, K) =
                        (LocTracerArray(L, ICell, K) *
                             LocPseudoThickCell(ICell, K) +
                         Dt * LocTracerTend(L, ICell, K)) /
                        LocProvPseudoThick(ICell, K);
                 });
          });
   }

   // compute tracer diffusion
   const Array2DReal &MeanPseudoThickEdge =
       AuxState->PseudoThicknessAux.MeanPseudoThickEdge;
   if (LocTracerDiffusion.Enabled) {
      Pacer::start("Tend:tracerDiffusion", 2);
      parallelForOuter(
          LaunchConfig({NTracers, Mesh->NCellsAll},
                       TeamScratch<Real>(VCoord->NVertLayers)),
          KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
             LocTracerDiffusion(Team, LocTracerTend, L, ICell, TracerArray,
                                MeanPseudoThickEdge);
          });
      Pacer::stop("Tend:tracerDiffusion", 2);
   }

   // compute tracer hyperdiffusion
   const Array3DReal &Del2TracersCell = AuxState->TracerAux.Del2TracersCell;
   if (LocTracerHyperDiff.Enabled) {
      Pacer::start("Tend:tracerHyperDiff", 2);
      parallelForOuter(
          LaunchConfig({NTracers, Mesh->NCellsAll},
                       TeamScratch<Real>(VCoord->NVertLayers)),
          KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
             LocTracerHyperDiff(Team, LocTracerTend, L, ICell, Del2TracersCell);
          });
      Pacer::stop("Tend:tracerHyperDiff", 2);
   }

   // compute tracer tendency from vertical advection
   Pacer::start("Tend:computeTracerVAdvTend", 2);
   Array2DReal ThicknessForVAdv;
   Array3DReal TracersForVAdv;
   if (VAdv->VertAdvChoice == VertAdvOption::Standard) {
      ThicknessForVAdv = State->getPseudoThickness(ThickTimeLevel);
      TracersForVAdv   = LocTracerArray;
   } else if (VAdv->VertAdvChoice == VertAdvOption::FCT) {
      ThicknessForVAdv = AuxState->PseudoThicknessAux.ProvPseudoThickness;
      TracersForVAdv   = ProvTracers;
   }
   VAdv->computeTracerVAdvTend(LocTracerTend, TracersForVAdv, ThicknessForVAdv,
                               TimeStep);
   Pacer::stop("Tend:computeTracerVAdvTend", 2);

   // compute tracer surface restoring
   const Array2DReal &TracersMonthlySurfClimo =
       AuxState->SurfTracerRestAux.TracersMonthlySurfClimoCell;
   const I4 NTracersToRestore = LocSurfaceTracerRestoring.NTracersToRestore;
   const auto &TracerIdsToRestore =
       LocSurfaceTracerRestoring.TracerIdsToRestore;
   if (LocSurfaceTracerRestoring.Enabled && NTracersToRestore > 0) {
      Pacer::start("Tend:surfaceTracerRestoring", 2);
      parallelFor(
          {NTracersToRestore, Mesh->NCellsAll},
          KOKKOS_LAMBDA(int R, int ICell) {
             const int KMin = LocMinLayerCell(ICell);
             const int L    = TracerIdsToRestore(R);
             LocSurfaceTracerRestoring(LocTracerTend, L, ICell, KMin,
                                       TracersMonthlySurfClimo, LocTracerArray);
          });
      Pacer::stop("Tend:surfaceTracerRestoring", 2);
   }

   // compute tracer forcing tendency; this also fills SurfaceTracerFlux,
   // which the KPP non-local flux term below depends on
   const auto *ForcingState = Forcing::getDefault();
   OMEGA_SCOPE(LocExtraEnergy, ForcingState->TracerForcing.ExtraEnergyCell);
   deepCopy(LocExtraEnergy, 0.0_Real);
   if (LocSfcTracerForcing.Enabled) {
      Pacer::start("Tend:sfcTracerForcing", 2);
      const auto &LatentHeatFluxEvap =
          ForcingState->TracerForcing.LatentHeatFluxEvapCell;
      const auto &SensibleHeatFlux =
          ForcingState->TracerForcing.SensibleHeatFluxCell;
      const auto &LongWaveHeatFluxUp =
          ForcingState->TracerForcing.LongWaveHeatFluxUpCell;
      const auto &LongWaveHeatFluxDown =
          ForcingState->TracerForcing.LongWaveHeatFluxDownCell;
      const auto &SeaIceHeatFlux =
          ForcingState->TracerForcing.SeaIceHeatFluxCell;
      const auto &ShortWaveHeatFlux =
          ForcingState->TracerForcing.ShortWaveHeatFluxCell;
      const auto &SnowFlux = ForcingState->TracerForcing.SnowFluxCell;
      const auto &RainFlux = ForcingState->TracerForcing.RainFluxCell;
      const auto &EvaporationFlux =
          ForcingState->TracerForcing.EvaporationFluxCell;
      const auto &IceRunoffFlux = ForcingState->TracerForcing.IceRunoffFluxCell;
      const auto &RiverRunoffFlux =
          ForcingState->TracerForcing.RiverRunoffFluxCell;
      const auto &SeaIceSaltFlux =
          ForcingState->TracerForcing.SeaIceSaltFluxCell;
      const auto &PressureMid = VCoord->PressureMid;
      const auto &SurfaceTracerFlux =
          ForcingState->TracerForcing.SurfaceTracerFluxCell;

      parallelFor(
          {Mesh->NCellsOwned}, KOKKOS_LAMBDA(int ICell) {
             LocSfcTracerForcing(
                 LocTracerTend, SurfaceTracerFlux, LocExtraEnergy, ICell,
                 TracerArray, LatentHeatFluxEvap, SensibleHeatFlux,
                 LongWaveHeatFluxUp, LongWaveHeatFluxDown, SeaIceHeatFlux,
                 ShortWaveHeatFlux, SnowFlux, RainFlux, IceRunoffFlux,
                 RiverRunoffFlux, EvaporationFlux, SeaIceSaltFlux);
          });
      Pacer::stop("Tend:sfcTracerForcing", 2);
   }

   // compute frazil tendency
   if (FrazilTerm.Enabled) {
      Pacer::start("Tend:frazil", 2);
      const auto &PressureMid     = VCoord->PressureMid;
      Array2DReal PseudoThickness = State->getPseudoThickness(ThickTimeLevel);
      FrazilTerm(PseudoThicknessTend, TracerTend, TracerArray, PressureMid,
                 PseudoThickness, TimeStep, FinalUpdateWeight);
      Pacer::stop("Tend:frazil", 2);
   }

   // Compute KPP non-local tracer tendency
   if (TracerVertMixTendencyEnable && LocKPPNonLocalTracerFlux.Enabled &&
       LocSfcTracerForcing.Enabled) {
      KPPMix *KPPInstance = KPPMix::getInstance();
      if (KPPInstance && KPPInstance->Enabled) {
         const auto *ForcingState = Forcing::getDefault();
         OMEGA_REQUIRE(ForcingState,
                       "KPP non-local tracer tendency requires Forcing");
         const auto &SurfaceTracerFlux =
             ForcingState->TracerForcing.SurfaceTracerFluxCell;
         Pacer::start("Tend:KPPNonLocalTracerFlux", 2);
         OMEGA_SCOPE(VertNonLocalFlux, KPPInstance->VertNonLocalFlux);
         parallelForOuter(
             {NTracers, Mesh->NCellsAll},
             KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
                LocKPPNonLocalTracerFlux(Team, LocTracerTend, L, ICell,
                                         SurfaceTracerFlux, VertNonLocalFlux);
             });

         if (LocKPPNonLocalTracerDiagnosticsEnable && LocHasTempTracer) {
            parallelForOuter(
                {Mesh->NCellsAll},
                KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
                   const int KMin = LocMinLayerCell(ICell);
                   const int KMax = LocMaxLayerCell(ICell);
                   parallelForInner(
                       Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                          LocKPPNonLocalTracerTempTendDiag(ICell, K) +=
                              SurfaceTracerFlux(LocTempTracerIndex, ICell) *
                              (VertNonLocalFlux(ICell, K) -
                               VertNonLocalFlux(ICell, K + 1));
                       });
                });

            parallelForOuter(
                {Mesh->NCellsAll},
                KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
                   const int KMin = LocMinLayerCell(ICell);
                   const int KMax = LocMaxLayerCell(ICell);
                   Real Sum       = 0.0_Real;
                   parallelReduceInner(
                       Team, Range{KMin, KMax},
                       [&](int K, Real &LocalSum) {
                          const Real NonLocalTend =
                              SurfaceTracerFlux(LocTempTracerIndex, ICell) *
                              (VertNonLocalFlux(ICell, K) -
                               VertNonLocalFlux(ICell, K + 1));
                          LocalSum += NonLocalTend;
                       },
                       Sum);
                   Kokkos::single(Kokkos::PerTeam(Team), [&]() {
                      LocKPPNonLocalTracerTempColumnSumDiag(ICell) = Sum;
                   });
                });
         }

         Pacer::stop("Tend:KPPNonLocalTracerFlux", 2);
      }
   }
   OMEGA_SCOPE(LocPenetratingShortwave, PenetratingShortwave);
   if (LocPenetratingShortwave.Enabled) {
      Pacer::start("Tend:penetratingShortwave", 2);
      const auto *ForcingState = Forcing::getDefault();
      const auto &ShortWaveHeatFlux =
          ForcingState->TracerForcing.ShortWaveHeatFluxCell;
      const auto &GeomZInterface = VCoord->GeomZInterface;
      const auto &ExtinctionCoeffRed =
          ForcingState->ShortwavePenForcing.ExtinctionCoeffRedCell;
      const auto &ExtinctionCoeffBlue =
          ForcingState->ShortwavePenForcing.ExtinctionCoeffBlueCell;

      parallelForOuter(
          LaunchConfig({Mesh->NCellsOwned},
                       TeamScratch<Real>(VCoord->NVertLayers + 1)),
          KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
             LocPenetratingShortwave(Team, LocTracerTend, ICell, GeomZInterface,
                                     ShortWaveHeatFlux, ExtinctionCoeffRed,
                                     ExtinctionCoeffBlue);
          });
      Pacer::stop("Tend:penetratingShortwave", 2);
   }

   Pacer::stop("Tend:computeTracerTendenciesOnly", 1);
} // end tracer tendency compute

void Tendencies::computePseudoThicknessTendencies(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    const Array3DReal &TracerArray, ///< [in] Tracer array
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time,               ///< [in] Time
    TimeInterval ProjDt ///< [in] Time interval for projection over the
                        ///< current time stepper stage
) {
   // only need PseudoThicknessAux on edge
   Array2DReal PseudoThick   = State->getPseudoThickness(ThickTimeLevel);
   Array2DReal NormalVelEdge = State->getNormalVelocity(VelTimeLevel);
   OMEGA_SCOPE(PseudoThicknessAux, AuxState->PseudoThicknessAux);
   OMEGA_SCOPE(PseudoThickCell, PseudoThick);
   OMEGA_SCOPE(LocNormalVelEdge, NormalVelEdge);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   Pacer::start("Tend:computePseudoThicknessTendencies", 1);

   AuxState->computePseudoThicknessAux(State, TracerArray, ThickTimeLevel,
                                       VelTimeLevel, ProjDt);

   computePseudoThicknessTendenciesOnly(State, AuxState, ThickTimeLevel,
                                        VelTimeLevel, Time);

   Pacer::stop("Tend:computePseudoThicknessTendencies", 1);
}

void Tendencies::computeVelocityTendencies(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    const Array3DReal &TracerArray, ///< [in] Tracer array
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    int TracerTimeLevel,            ///< [in] Time level
    TimeInstant Time,               ///< [in] Time
    TimeInterval ProjDt ///< [in] Time interval for projection over the
                        ///< current time stepper stage
) {
   Pacer::start("Tend:computeVelocityTendencies", 1);

   AuxState->computeMomAux(State, TracerArray, ThickTimeLevel, VelTimeLevel,
                           ProjDt);
   computeVelocityTendenciesOnly(State, AuxState, TracerArray, ThickTimeLevel,
                                 VelTimeLevel, TracerTimeLevel, Time);

   Pacer::stop("Tend:computeVelocityTendencies", 1);
}

void Tendencies::computeTracerTendencies(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    const Array3DReal &TracerArray, ///< [in] Tracer array
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time,               ///< [in] Time
    TimeInterval ProjDt,   ///< [in] Time interval for projection over the
                           ///< current time stepper stage
    Real FinalUpdateWeight ///< [in] final update weight
) {

   Pacer::start("Tend:computeTracerTendencies", 1);

   AuxState->computeTracerAux(State, TracerArray, ThickTimeLevel, VelTimeLevel,
                              ProjDt);

   computeTracerTendenciesOnly(State, AuxState, TracerArray, ThickTimeLevel,
                               VelTimeLevel, Time, FinalUpdateWeight);

   Pacer::stop("Tend:computeTracerTendencies", 1);
}

//------------------------------------------------------------------------------
// Compute both pseudo-thickness and normal velocity tendencies
void Tendencies::computeAllTendencies(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    const Array3DReal &TracerArray, ///< [in] Tracer array
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    int TracerTimeLevel,            ///< [in] Time level
    TimeInstant Time,               ///< [in] Time
    TimeInterval ProjDt,   ///< [in] Time interval for projection over the
                           ///< current time stepper stage
    Real FinalUpdateWeight ///< [in] final update weight
) {
   AuxState->computeAll(State, TracerArray, ThickTimeLevel, VelTimeLevel,
                        ProjDt);

   computePseudoThicknessTendenciesOnly(State, AuxState, ThickTimeLevel,
                                        VelTimeLevel, Time);
   computeVelocityTendenciesOnly(State, AuxState, TracerArray, ThickTimeLevel,
                                 VelTimeLevel, TracerTimeLevel, Time);
   computeTracerTendenciesOnly(State, AuxState, TracerArray, ThickTimeLevel,
                               VelTimeLevel, Time, FinalUpdateWeight);
} // end all tendency compute

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
