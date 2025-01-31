//===-- ocn/CustomTendencyTerms.cpp - Custom tendency terms -----*- C++ -*-===//
//
// The customized tendency terms can be added to the tendency terms based
// based on an option 'UseCustomTendency' in Tendencies Config group.
// This file contains functions for initializing customized tendency terms.
//
//===----------------------------------------------------------------------===//

#include "CustomTendencyTerms.h"
#include "Config.h"
#include "TimeStepper.h"

namespace OMEGA {

//===-----------------------------------------------------------------------===/
// Initialize the manufactured solution tendency terms.
//===-----------------------------------------------------------------------===/
int ManufacturedSolution::init() {
   int Err;

   // Get ManufacturedSolConfig group
   Config *OmegaConfig = Config::getOmegaConfig();
   Config ManufacturedSolConfig("ManufacturedSolution");
   Err = OmegaConfig->get(ManufacturedSolConfig);
   if (Err != 0) {
      LOG_CRITICAL("ManufacturedSolution:: ManufacturedSolution group "
                   "not found in Config");
      return Err;
   }

   // Get TendConfig group
   Config TendConfig("Tendencies");
   Err = OmegaConfig->get(TendConfig);
   if (Err != 0) {
      LOG_CRITICAL("ManufacturedSolution:: Tendencies group "
                   "not found in Config");
      return Err;
   }

   // Get manufactured solution parameters from Config
   R8 WavelengthX;
   R8 WavelengthY;
   R8 Amplitude;

   Err = ManufacturedSolConfig.get("WavelengthX", WavelengthX);
   if (Err != 0) {
      LOG_ERROR("ManufacturedSolution:: WavelengthX not found in "
                "ManufacturedSolConfig");
      return Err;
   }

   Err = ManufacturedSolConfig.get("WavelengthY", WavelengthY);
   if (Err != 0) {
      LOG_ERROR("ManufacturedSolution:: WavelengthY not found in "
                "ManufacturedSolConfig");
      return Err;
   }

   Err = ManufacturedSolConfig.get("Amplitude", Amplitude);
   if (Err != 0) {
      LOG_ERROR("ManufacturedSolution:: Amplitude not found in "
                "ManufacturedSolConfig");
      return Err;
   }

   // Get Tendendices parameters for del2 and del4 source terms
   Err = TendConfig.get("VelDiffTendencyEnable",
                        ManufacturedVelTend.VelDiffTendencyEnable);
   Err += TendConfig.get("VelHyperDiffTendencyEnable",
                         ManufacturedVelTend.VelHyperDiffTendencyEnable);
   Err += TendConfig.get("ViscDel2", ManufacturedVelTend.ViscDel2);
   Err += TendConfig.get("ViscDel4", ManufacturedVelTend.ViscDel4);

   if (Err != 0) {
      LOG_ERROR("ManufacturedSolution::Error reading Tendencies config");
      return Err;
   }

   // Get the reference time to compute the model elapsed time
   /// Get model clock from time stepper
   TimeStepper *DefStepper             = TimeStepper::getDefault();
   Clock *ModelClock                   = DefStepper->getClock();
   ManufacturedThickTend.ReferenceTime = ModelClock->getCurrentTime();
   ManufacturedVelTend.ReferenceTime   = ManufacturedThickTend.ReferenceTime;

   // Get BottomDepth for the resting thickness
   /// This test case assumes that the restingThickness is horizontally uniform
   /// and that only one vertical level is used so only one set of indices is
   /// used here.
   HorzMesh *DefHorzMesh = HorzMesh::getDefault();
   R8 H0                 = DefHorzMesh->BottomDepthH(0);

   // Define and compute common constants
   R8 Grav    = 9.80665_Real;                          // Gravity acceleration
   R8 Pii     = 3.141592653589793_Real;                // Pi
   R8 Kx      = 2.0_Real * Pii / WavelengthX;          // Wave in X-dir
   R8 Ky      = 2.0_Real * Pii / WavelengthY;          // Wave in Y-dir
   R8 AngFreq = sqrt(H0 * Grav * (Kx * Kx + Ky * Ky)); // Angular frequency

   // Assign constants for thickness tendency function
   ManufacturedThickTend.H0      = H0;
   ManufacturedThickTend.Eta0    = Amplitude;
   ManufacturedThickTend.Kx      = Kx;
   ManufacturedThickTend.Ky      = Ky;
   ManufacturedThickTend.AngFreq = AngFreq;

   // Assign constants for velocity tendency function
   ManufacturedVelTend.Grav    = Grav;
   ManufacturedVelTend.Eta0    = Amplitude;
   ManufacturedVelTend.Kx      = Kx;
   ManufacturedVelTend.Ky      = Ky;
   ManufacturedVelTend.AngFreq = AngFreq;

   return Err;

} // end ManufacturedSolution init

//===--------------------------------------------------------------------===/
// Manufactured tendency term for the thickness equation
//===--------------------------------------------------------------------===/
void ManufacturedSolution::ManufacturedThicknessTendency::operator()(
    Array2DReal ThicknessTend, const OceanState *State,
    const AuxiliaryState *AuxState, int ThickTimeLevel, int VelTimeLevel,
    TimeInstant Time) const {

   // Get elapsed time since reference time
   R8 ElapsedTimeSec;
   TimeInterval ElapsedTimeInterval = Time - ReferenceTime;
   ElapsedTimeInterval.get(ElapsedTimeSec, TimeUnits::Seconds);

   auto *Mesh       = HorzMesh::getDefault();
   auto NVertLevels = ThicknessTend.extent_int(1);

   Array1DReal XCell = Mesh->XCell;
   Array1DReal YCell = Mesh->YCell;

   OMEGA_SCOPE(LocH0, H0);
   OMEGA_SCOPE(LocEta0, Eta0);
   OMEGA_SCOPE(LocKx, Kx);
   OMEGA_SCOPE(LocKy, Ky);
   OMEGA_SCOPE(LocAngFreq, AngFreq);

   parallelFor(
       {Mesh->NCellsAll, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KLevel) {
          R8 X     = XCell(ICell);
          R8 Y     = YCell(ICell);
          R8 Phase = LocKx * X + LocKy * Y - LocAngFreq * ElapsedTimeSec;
          ThicknessTend(ICell, KLevel) +=
              LocEta0 *
              (-LocH0 * (LocKx + LocKy) * sin(Phase) - LocAngFreq * cos(Phase) +
               LocEta0 * (LocKx + LocKy) * cos(2.0_Real * Phase));
       });

} // end void ManufacturedThicknessTendency

//===--------------------------------------------------------------------===/
// Manufactured tendency term for the momentum equation
//===--------------------------------------------------------------------===/
void ManufacturedSolution::ManufacturedVelocityTendency::operator()(
    Array2DReal NormalVelTend, const OceanState *State,
    const AuxiliaryState *AuxState, int ThickTimeLevel, int VelTimeLevel,
    TimeInstant Time) const {

   // Get elapsed time since reference time
   R8 ElapsedTimeSec;
   TimeInterval ElapsedTimeInterval = Time - ReferenceTime;
   ElapsedTimeInterval.get(ElapsedTimeSec, TimeUnits::Seconds);

   auto *Mesh       = HorzMesh::getDefault();
   auto NVertLevels = NormalVelTend.extent_int(1);

   Array1DReal FEdge     = Mesh->FEdge;
   Array1DReal XEdge     = Mesh->XEdge;
   Array1DReal YEdge     = Mesh->YEdge;
   Array1DReal AngleEdge = Mesh->AngleEdge;

   OMEGA_SCOPE(LocGrav, Grav);
   OMEGA_SCOPE(LocEta0, Eta0);
   OMEGA_SCOPE(LocKx, Kx);
   OMEGA_SCOPE(LocKy, Ky);
   OMEGA_SCOPE(LocAngFreq, AngFreq);
   OMEGA_SCOPE(LocViscDel2, ViscDel2);
   OMEGA_SCOPE(LocViscDel4, ViscDel4);
   OMEGA_SCOPE(LocVelDiffTendencyEnable, VelDiffTendencyEnable);
   OMEGA_SCOPE(LocVelHyperDiffTendencyEnable, VelHyperDiffTendencyEnable);

   R8 LocKx2 = LocKx * LocKx;
   R8 LocKy2 = LocKy * LocKy;
   R8 LocKx4 = LocKx2 * LocKx2;
   R8 LocKy4 = LocKy2 * LocKy2;

   parallelFor(
       {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
          R8 X = XEdge(IEdge);
          R8 Y = YEdge(IEdge);

          R8 Phase       = LocKx * X + LocKy * Y - LocAngFreq * ElapsedTimeSec;
          R8 SourceTerm0 = LocAngFreq * sin(Phase) - 0.5_Real * LocEta0 *
                                                         (LocKx + LocKy) *
                                                         sin(2.0_Real * Phase);

          R8 U = LocEta0 *
                 ((-FEdge(IEdge) + LocGrav * LocKx) * cos(Phase) + SourceTerm0);
          R8 V = LocEta0 *
                 ((FEdge(IEdge) + LocGrav * LocKy) * cos(Phase) + SourceTerm0);

          // Del2 and del4 source terms
          if (LocVelDiffTendencyEnable) {
             U += LocViscDel2 * LocEta0 * (LocKx2 + LocKy2) * cos(Phase);
             V += LocViscDel2 * LocEta0 * (LocKx2 + LocKy2) * cos(Phase);
          }
          if (LocVelHyperDiffTendencyEnable) {
             U -= LocViscDel4 * LocEta0 *
                  ((LocKx4 + LocKy4 + LocKx2 * LocKy2) * cos(Phase));
             V -= LocViscDel4 * LocEta0 *
                  ((LocKx4 + LocKy4 + LocKx2 * LocKy2) * cos(Phase));
          }

          R8 NormalCompSourceTerm =
              cos(AngleEdge(IEdge)) * U + sin(AngleEdge(IEdge)) * V;
          NormalVelTend(IEdge, KLevel) += NormalCompSourceTerm;
       });

} // end void ManufacturedVelocityTendency

} // end namespace OMEGA

//=-------------------------------------------------------------------------===/
