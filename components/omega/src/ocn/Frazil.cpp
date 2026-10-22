//===-- ocn/Frazil.cpp - Frazil Ice Formation -----------------*- C++ -*-===//
//
// The Frazil class manages frazil tendencies and accumulators.
// The frazil options include a Teos-10 and a FixedProperty implementation.
//
//===----------------------------------------------------------------------===//

#include "Frazil.h"
#include "Eos.h"
#include "Error.h"
#include "Field.h"
#include "Logging.h"
#include "TimeStepper.h"

#include <limits>

namespace OMEGA {

namespace {

KOKKOS_INLINE_FUNCTION bool isApprox(const Real X, const Real Y,
                                     const Real RTol, const Real ATol = 0) {
   if (Kokkos::isnan(X) || Kokkos::isnan(Y) || Kokkos::isinf(X) ||
       Kokkos::isinf(Y)) {
      return false;
   }

   return Kokkos::abs(X - Y) <=
          Kokkos::max(ATol, RTol * Kokkos::max(Kokkos::abs(X), Kokkos::abs(Y)));
}

} // namespace

Frazil *Frazil::DefaultFrazil = nullptr;
std::map<std::string, std::unique_ptr<Frazil>> Frazil::AllFrazil;

/// Constructor for FrazilFormation
FrazilFormation::FrazilFormation() {}

/// Constructor for FrazilMelt
FrazilMelt::FrazilMelt() {}

/// Constructor for FixedPropertyFrazilFormation
FixedPropertyFrazilFormation::FixedPropertyFrazilFormation() {}

/// Constructor for FixedPropertyFrazilMelt
FixedPropertyFrazilMelt::FixedPropertyFrazilMelt() {}

void Frazil::init() {

   if (!HorzMesh::getDefault() or !VertCoord::getDefault()) {
      ABORT_ERROR("Frazil::init: HorzMesh and VertCoord must be initialized");
   }

   // Frazil freezing-temperature calculations depend on EOS configuration.
   Eos::init();

   if (!DefaultFrazil) {
      Error Err;
      bool FrazilTendencyEnable = false;
      Config *OmegaConfig       = Config::getOmegaConfig();
      Config TendConfig("Tendencies");

      Err += OmegaConfig->get(TendConfig);
      CHECK_ERROR_ABORT(Err,
                        "Frazil::init: Tendencies group not found in Config");

      Err += TendConfig.get("FrazilTendencyEnable", FrazilTendencyEnable);
      CHECK_ERROR_ABORT(
          Err, "Frazil::init: FrazilTendencyEnable not found in Tendencies");

      if (!FrazilTendencyEnable) {
         LOG_INFO("Frazil::init: Frazil tendency disabled; skipping default "
                  "frazil object creation");
         LOG_INFO("All frazil is off - frazil parameters will be ignored");
         return;
      }

      TimeStepper *DefTimeStepper = TimeStepper::getDefault();
      if (DefTimeStepper &&
          DefTimeStepper->getType() == TimeStepperType::ForwardBackward) {
         ABORT_ERROR(
             "Frazil is not supported for the Forward-Backward timestepper. "
             "Turn frazil off or use a different timestepper");
      }

      FieldGroup::create("Frazil");
      DefaultFrazil = create("Default");
      DefaultFrazil->registerFields();
   }
}

Frazil::Frazil(const HorzMesh *Mesh, const VertCoord *VCoord)
    : FrazilChoice(FrazilType::TeosFrazil), ComputeFrazilFormation(),
      ComputeFrazilMelt(), NCellsAll(Mesh->NCellsAll),
      NChunks((VCoord->NVertLayers + VecLength - 1) / VecLength), MeshPtr(Mesh),
      VCoordPtr(VCoord) {

   FrazilTTend =
       Array2DReal("FrazilTTend", Mesh->NCellsSize, VCoord->NVertLayers);
   FrazilSTend =
       Array2DReal("FrazilSTend", Mesh->NCellsSize, VCoord->NVertLayers);
   FrazilHTend =
       Array2DReal("FrazilHTend", Mesh->NCellsSize, VCoord->NVertLayers);

   AccMIce          = Array1DReal("AccMIce", Mesh->NCellsSize);
   AccEIce          = Array1DReal("AccEIce", Mesh->NCellsSize);
   AccMLiq          = Array1DReal("AccMLiq", Mesh->NCellsSize);
   AccELiq          = Array1DReal("AccELiq", Mesh->NCellsSize);
   AccMSalt         = Array1DReal("AccMSalt", Mesh->NCellsSize);
   FrazilMassFlux   = Array1DReal("FrazilMassFlux", Mesh->NCellsSize);
   FrazilSaltFlux   = Array1DReal("FrazilSaltFlux", Mesh->NCellsSize);
   FrazilEnergyFlux = Array1DReal("FrazilEnergyFlux", Mesh->NCellsSize);

   deepCopy(FrazilTTend, 0.0_Real);
   deepCopy(FrazilSTend, 0.0_Real);
   deepCopy(FrazilHTend, 0.0_Real);
   deepCopy(AccMIce, 0.0_Real);
   deepCopy(AccEIce, 0.0_Real);
   deepCopy(AccMLiq, 0.0_Real);
   deepCopy(AccELiq, 0.0_Real);
   deepCopy(AccMSalt, 0.0_Real);
   resetOcnStepFluxes();
}

Frazil::~Frazil() { unregisterFields(); }

Frazil *Frazil::create(const std::string &Name) {
   if (AllFrazil.find(Name) != AllFrazil.end()) {
      LOG_ERROR("Attempted to create Frazil {} but it already exists", Name);
      return nullptr;
   }

   auto *NewFrazil =
       new Frazil(HorzMesh::getDefault(), VertCoord::getDefault());
   AllFrazil.emplace(Name, NewFrazil);

   Error Err;
   Config *OmegaConfig = Config::getOmegaConfig();
   Config FrazilConfig("Frazil");
   Err += OmegaConfig->get(FrazilConfig);
   CHECK_ERROR_ABORT(Err, "Frazil::create: Frazil group not found in Config");

   std::string FrazilTypeStr;
   Err += FrazilConfig.get("FrazilType", FrazilTypeStr);
   CHECK_ERROR_ABORT(Err,
                     "Frazil::create: FrazilType not found in Frazil config");

   if ((FrazilTypeStr == "FixedProperty") or
       (FrazilTypeStr == "fixedProperty") or (FrazilTypeStr == "fixed") or
       (FrazilTypeStr == "fixedproperty")) {
      NewFrazil->FrazilChoice = FrazilType::FixedPropertyFrazil;
   } else if ((FrazilTypeStr == "teos") or (FrazilTypeStr == "Teos") or
              (FrazilTypeStr == "TEOS") or (FrazilTypeStr == "Teos10") or
              (FrazilTypeStr == "teos10") or (FrazilTypeStr == "TEOS10")) {
      NewFrazil->FrazilChoice = FrazilType::TeosFrazil;
   } else {
      ABORT_ERROR("Frazil::create: Unknown FrazilType requested");
   }

   Err += FrazilConfig.get("LayerMassFracMax",
                           NewFrazil->ComputeFrazilFormation.LayerMassFracMax);
   Err += FrazilConfig.get("LayerMassFracMax",
                           NewFrazil->ComputeFrazilMelt.LayerMassFracMax);
   CHECK_ERROR_ABORT(
       Err, "Frazil::create: LayerMassFracMax not found in Frazil config");
   if (!(NewFrazil->ComputeFrazilFormation.LayerMassFracMax > 0.0_Real &&
         NewFrazil->ComputeFrazilFormation.LayerMassFracMax < 1.0_Real)) {
      ABORT_ERROR("Frazil::create: LayerMassFracMax must be between 0 and 1 "
                  "(excluded)");
   }

   NewFrazil->ComputeFixedPropertyFrazilFormation.LayerMassFracMax =
       NewFrazil->ComputeFrazilFormation.LayerMassFracMax;
   NewFrazil->ComputeFixedPropertyFrazilMelt.LayerMassFracMax =
       NewFrazil->ComputeFrazilFormation.LayerMassFracMax;

   Err += FrazilConfig.get("Phi", NewFrazil->ComputeFrazilFormation.Phi);
   CHECK_ERROR_ABORT(Err, "Frazil::create: Phi not found in Frazil config");
   if (!(NewFrazil->ComputeFrazilFormation.Phi >= 0.0_Real &&
         NewFrazil->ComputeFrazilFormation.Phi < 1.0_Real)) {
      ABORT_ERROR("Frazil::create: Phi must be between 0 and 1 (1 excluded)");
   }

   Err += FrazilConfig.get("ConservationCheck", NewFrazil->ConservationCheck);
   CHECK_ERROR_ABORT(
       Err, "Frazil::create: ConservationCheck not found in Frazil config");

   Err += FrazilConfig.get("DepthLimit", NewFrazil->DepthLimit);
   CHECK_ERROR_ABORT(Err,
                     "Frazil::create: DepthLimit not found in Frazil config");

   if (Name == "Default") {
      DefaultFrazil = NewFrazil;
   }

   return NewFrazil;
}

Frazil *Frazil::getDefault() { return DefaultFrazil; }

Frazil *Frazil::get(const std::string &Name) {
   auto Iter = AllFrazil.find(Name);
   if (Iter != AllFrazil.end()) {
      return Iter->second.get();
   }

   LOG_ERROR("Frazil::get: Attempted to retrieve non-existent Frazil {}", Name);
   return nullptr;
}

void Frazil::erase(std::string InName) {
   auto *ToErase = get(InName);
   AllFrazil.erase(InName);
   if (ToErase == DefaultFrazil) {
      DefaultFrazil = nullptr;
   }
}

void Frazil::clear() {
   AllFrazil.clear();
   DefaultFrazil = nullptr;
   if (FieldGroup::exists("Frazil")) {
      FieldGroup::destroy("Frazil");
   }
}

void Frazil::resetOcnStepFluxes() {
   deepCopy(FrazilMassFlux, 0.0_Real);
   deepCopy(FrazilSaltFlux, 0.0_Real);
   deepCopy(FrazilEnergyFlux, 0.0_Real);
}

void Frazil::accumulateOcnStepFluxes(const Real FinalUpdateWeight,
                                     const R8 TimeStepSeconds) {
   // This function accumulates the frazil terms produced at each frazil call
   // (multiple per timestep) into the ocean (outer)timestep fluxes.
   OMEGA_SCOPE(LocAccMIce, AccMIce);
   OMEGA_SCOPE(LocAccMLiq, AccMLiq);
   OMEGA_SCOPE(LocAccMSalt, AccMSalt);
   OMEGA_SCOPE(LocAccELiq, AccELiq);
   OMEGA_SCOPE(LocAccEIce, AccEIce);
   OMEGA_SCOPE(LocFrazilMassFlux, FrazilMassFlux);
   OMEGA_SCOPE(LocFrazilSaltFlux, FrazilSaltFlux);
   OMEGA_SCOPE(LocFrazilEnergyFlux, FrazilEnergyFlux);

   // Accumulators are positive into the frazil reservoir; the exported fluxes
   // use the coupler convention of being counted positive into the ocean.
   parallelFor(
       {NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          LocFrazilMassFlux(ICell) -= FinalUpdateWeight *
                                      (LocAccMIce(ICell) + LocAccMLiq(ICell)) /
                                      TimeStepSeconds;
          LocFrazilSaltFlux(ICell) -=
              FinalUpdateWeight * LocAccMSalt(ICell) / TimeStepSeconds;
          LocFrazilEnergyFlux(ICell) -=
              FinalUpdateWeight * (LocAccEIce(ICell) + LocAccELiq(ICell)) /
              TimeStepSeconds;
       });
}

void Frazil::registerFields() {
   constexpr int NDims = 1;
   const std::vector<std::string> DimNames{"NCells"};

   auto FrazilMassField = Field::create(
       FrazilMassFlux.label(),
       "frazil mass flux averaged over ocean timestep, positive into the ocean",
       "kg m^-2 s^-1", "", std::numeric_limits<Real>::lowest(),
       std::numeric_limits<Real>::max(), NDims, DimNames);
   auto FrazilSaltField = Field::create(
       FrazilSaltFlux.label(),
       "frazil salt flux averaged over ocean timestep, positive into the ocean",
       "kg m^-2 s^-1", "", std::numeric_limits<Real>::lowest(),
       std::numeric_limits<Real>::max(), NDims, DimNames);
   auto FrazilEnergyField =
       Field::create(FrazilEnergyFlux.label(),
                     "frazil energy flux averaged over ocean timestep, "
                     "positive into the ocean",
                     "W m^-2", "", std::numeric_limits<Real>::lowest(),
                     std::numeric_limits<Real>::max(), NDims, DimNames);

   FieldGroup::addFieldToGroup(FrazilMassFlux.label(), "Frazil");
   FieldGroup::addFieldToGroup(FrazilSaltFlux.label(), "Frazil");
   FieldGroup::addFieldToGroup(FrazilEnergyFlux.label(), "Frazil");

   FrazilMassField->attachData<Array1DReal>(FrazilMassFlux);
   FrazilSaltField->attachData<Array1DReal>(FrazilSaltFlux);
   FrazilEnergyField->attachData<Array1DReal>(FrazilEnergyFlux);
   FieldsRegistered = true;
}

void Frazil::unregisterFields() {
   if (!FieldsRegistered) {
      return;
   }

   if (Field::exists(FrazilMassFlux.label())) {
      Field::destroy(FrazilMassFlux.label());
   }
   if (Field::exists(FrazilSaltFlux.label())) {
      Field::destroy(FrazilSaltFlux.label());
   }
   if (Field::exists(FrazilEnergyFlux.label())) {
      Field::destroy(FrazilEnergyFlux.label());
   }
   FieldsRegistered = false;
}

void Frazil::checkColumnConservation() const {
   auto MinLayerCellH = createHostMirrorCopy(VCoordPtr->MinLayerCell);
   auto MaxLayerCellH = createHostMirrorCopy(VCoordPtr->MaxLayerCell);
   auto FrazilHTendH  = createHostMirrorCopy(FrazilHTend);
   auto FrazilTTendH  = createHostMirrorCopy(FrazilTTend);
   auto FrazilSTendH  = createHostMirrorCopy(FrazilSTend);
   auto AccMIceH      = createHostMirrorCopy(AccMIce);
   auto AccMLiqH      = createHostMirrorCopy(AccMLiq);
   auto AccMSaltH     = createHostMirrorCopy(AccMSalt);
   auto AccELiqH      = createHostMirrorCopy(AccELiq);
   auto AccEIceH      = createHostMirrorCopy(AccEIce);

   constexpr Real RTol = sizeof(Real) == 4 ? 1.0e-5_Real : 1.0e-10_Real;
   constexpr Real ATol = sizeof(Real) == 4 ? 1.0e-4_Real : 1.0e-8_Real;

   for (I4 ICell = 0; ICell < NCellsAll; ++ICell) {
      const I4 KMin = MinLayerCellH(ICell);
      const I4 KMax = MaxLayerCellH(ICell);

      Real MassTend   = 0.0_Real;
      Real EnergyTend = 0.0_Real;
      Real SaltTend   = 0.0_Real;

      for (I4 K = KMin; K <= KMax; ++K) {
         MassTend += FrazilHTendH(ICell, K);
         EnergyTend += FrazilTTendH(ICell, K);
         SaltTend += FrazilSTendH(ICell, K);
      }

      const Real MassTotal   = AccMIceH(ICell) + AccMLiqH(ICell);
      const Real EnergyTotal = AccELiqH(ICell) + AccEIceH(ICell);
      const Real SaltTotal   = AccMSaltH(ICell);
      if (ICell == 0) {
         LOG_INFO("Frazil column conservation check: cell {} MassTend={} "
                  "MassTotal={} "
                  "EnergyTend={} EnergyTotal={} SaltTend={} SaltTotal={}",
                  ICell, MassTend * RhoSw, MassTotal,
                  EnergyTend * Cp0Sw * RhoSw, EnergyTotal,
                  SaltTend * RhoSw * PPt2Salt, SaltTotal);
         LOG_INFO("Frazil column conservation check: cell {} EpsMass={} "
                  "EpsE={} EpsS={} ",
                  ICell, MassTend * RhoSw + MassTotal,
                  EnergyTend * Cp0Sw * RhoSw + EnergyTotal,
                  SaltTend * RhoSw * PPt2Salt + SaltTotal);
      }

      if (!isApprox(-MassTend * RhoSw, MassTotal, RTol, ATol)) {
         ABORT_ERROR(
             "Frazil column mass check failed: cell {} tendency={} total={}",
             ICell, -MassTend * RhoSw, MassTotal);
      }
      if (!isApprox(-EnergyTend * Cp0Sw * RhoSw, EnergyTotal, RTol, ATol)) {
         ABORT_ERROR(
             "Frazil column energy check failed: cell {} tendency={} total={}",
             ICell, -EnergyTend * Cp0Sw * RhoSw, EnergyTotal);
      }
      if (!isApprox(-SaltTend * RhoSw * PPt2Salt, SaltTotal, RTol, ATol)) {
         ABORT_ERROR(
             "Frazil column salt check failed: cell {} tendency={} total={}",
             ICell, -SaltTend * RhoSw * PPt2Salt, SaltTotal);
      }
   }
}

// Reported once per run; a persistent negative salinity would otherwise warn
// on every tendency evaluation.
void Frazil::warnNegativeSalinity(I4 NClamped) {
   if (NClamped > 0 && !WarnedNegativeSalinity) {
      WarnedNegativeSalinity = true;
      LOG_WARN("Frazil: clamped negative absolute salinity to zero in {} "
               "layer(s); the freezing-point polynomial is undefined for "
               "SA < 0. This warning is issued only once.",
               NClamped);
   }
}

void Frazil::computeFrazilFixedPropertyImpl(
    const Array2DReal &ConservTemp, const Array2DReal &AbsSalinity,
    const Array2DReal &Pressure, const Array2DReal &PseudoThickness) {
   const EosType LocEosChoice = Eos::getInstance()->EosChoice;
   const Real LocDepthLimit   = DepthLimit;

   OMEGA_SCOPE(MinLayerCell, VCoordPtr->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoordPtr->MaxLayerCell);
   OMEGA_SCOPE(LocGeomZMid, VCoordPtr->GeomZMid);

   OMEGA_SCOPE(LocComputeFixedPropertyFrazilFormation,
               ComputeFixedPropertyFrazilFormation);
   OMEGA_SCOPE(LocComputeFixedPropertyFrazilMelt,
               ComputeFixedPropertyFrazilMelt);
   OMEGA_SCOPE(LocFrazilTTend, FrazilTTend);
   OMEGA_SCOPE(LocFrazilSTend, FrazilSTend);
   OMEGA_SCOPE(LocFrazilHTend, FrazilHTend);
   OMEGA_SCOPE(LocAccMIce, AccMIce);
   OMEGA_SCOPE(LocAccEIce, AccEIce);
   OMEGA_SCOPE(LocAccMLiq, AccMLiq);
   OMEGA_SCOPE(LocAccELiq, AccELiq);
   OMEGA_SCOPE(LocAccMSalt, AccMSalt);
   OMEGA_SCOPE(LocIceRefSal, IceRefSal);
   OMEGA_SCOPE(LocLatIce, LatIce);

   I4 NClamped = 0;
   parallelReduce(
       "frazilFixedPropertyImpl", {NCellsAll},
       KOKKOS_LAMBDA(I4 ICell, I4 & Accum) {
          const I4 KMin = MinLayerCell(ICell);
          const I4 KMax = MaxLayerCell(ICell);

          I4 KLim          = KMax;
          bool HasKLim     = true;
          const bool Limit = (LocDepthLimit >= 0.0_Real);

          if (Limit) {
             HasKLim = false;
             for (I4 K = KMax; K >= KMin; --K) {
                if (Kokkos::abs(LocGeomZMid(ICell, K)) <= LocDepthLimit) {
                   KLim    = K;
                   HasKLim = true;
                   break;
                }
             }
          }

          // Explicit accumulation order: bottom layer to top layer.
          for (I4 K = KMax; K >= KMin; --K) {
             if (!HasKLim || K > KLim) {
                LocFrazilHTend(ICell, K) = 0.0_Real;
                LocFrazilTTend(ICell, K) = 0.0_Real;
                LocFrazilSTend(ICell, K) = 0.0_Real;
                continue;
             }

             // sqrt(SA) in the freezing-point polynomial returns NaN for
             // SA < 0, which silently disables frazil since NaN compares false
             const Real AbsSalinityRaw = AbsSalinity(ICell, K);
             if (AbsSalinityRaw < 0.0_Real) {
                ++Accum;
             }
             const Real AbsSalinityIn = Kokkos::max(0.0_Real, AbsSalinityRaw);
             const Real ConservTempIn = ConservTemp(ICell, K);
             const Real PressureIn    = Pressure(ICell, K);
             const Real PressureDb    = PressureIn * Pa2Db;
             const Real PseudoThicknessIn = PseudoThickness(ICell, K);

             const Real CtFreezing = Eos::calcCtFreezing(
                 LocEosChoice, AbsSalinityIn, PressureDb, 0.0_Real);

             Real HTend = 0.0_Real;
             Real TTend = 0.0_Real;
             Real STend = 0.0_Real;

             if (ConservTempIn < CtFreezing) {
                LocComputeFixedPropertyFrazilFormation(
                    AbsSalinityIn, ConservTempIn, PressureDb, PseudoThicknessIn,
                    LocAccMIce(ICell), LocAccMSalt(ICell), LocAccEIce(ICell),
                    HTend, TTend, STend, CtFreezing);
             } else if (LocAccMIce(ICell) > 0.0_Real) {
                LocComputeFixedPropertyFrazilMelt(
                    AbsSalinityIn, ConservTempIn, PressureDb, PseudoThicknessIn,
                    LocAccMIce(ICell), LocAccMSalt(ICell), LocAccEIce(ICell),
                    HTend, TTend, STend, CtFreezing);
             }

             // Per-call increments; FrazilOnCell normalizes these to rates.
             LocFrazilHTend(ICell, K) = HTend;
             LocFrazilTTend(ICell, K) = TTend;
             LocFrazilSTend(ICell, K) = STend;
          } // end of vertical loop

          // in the fractional / conservative treatment of melt
          // there is no excess salt etc. or need to hijack coupling terms

          // Convert to coupler units
          LocAccMIce(ICell)  = LocAccMIce(ICell) * RhoSw;
          LocAccMLiq(ICell)  = LocAccMLiq(ICell) * RhoSw;
          LocAccMSalt(ICell) = LocAccMSalt(ICell) * RhoSw * PPt2Salt;
          LocAccELiq(ICell)  = LocAccELiq(ICell) * RhoSw;
          LocAccEIce(ICell)  = LocAccEIce(ICell) * RhoSw;
       },
       NClamped); // end of NCells loop

   warnNegativeSalinity(NClamped);
}

// TEOS-10 frazil relies on host-only GSW routines, so this implementation
// mirrors all inputs/outputs to the host and runs a plain host loop.
// This is a temporary implementation until a device-callable solution
// is available.
void Frazil::computeFrazilTeosImpl(const Array2DReal &ConservTemp,
                                   const Array2DReal &AbsSalinity,
                                   const Array2DReal &Pressure,
                                   const Array2DReal &PseudoThickness) {
   const EosType LocEosChoice = Eos::getInstance()->EosChoice;
   const Real LocDepthLimit   = DepthLimit;

   const auto MinLayerCell = createHostMirrorCopy(VCoordPtr->MinLayerCell);
   const auto MaxLayerCell = createHostMirrorCopy(VCoordPtr->MaxLayerCell);
   const auto LocGeomZMid  = createHostMirrorCopy(VCoordPtr->GeomZMid);

   const auto AbsSalinityH     = createHostMirrorCopy(AbsSalinity);
   const auto ConservTempH     = createHostMirrorCopy(ConservTemp);
   const auto PressureH        = createHostMirrorCopy(Pressure);
   const auto PseudoThicknessH = createHostMirrorCopy(PseudoThickness);

   auto LocFrazilTTend = createHostMirrorCopy(FrazilTTend);
   auto LocFrazilSTend = createHostMirrorCopy(FrazilSTend);
   auto LocFrazilHTend = createHostMirrorCopy(FrazilHTend);
   auto LocAccMIce     = createHostMirrorCopy(AccMIce);
   auto LocAccEIce     = createHostMirrorCopy(AccEIce);
   auto LocAccMLiq     = createHostMirrorCopy(AccMLiq);
   auto LocAccELiq     = createHostMirrorCopy(AccELiq);
   auto LocAccMSalt    = createHostMirrorCopy(AccMSalt);

   // Copies of the functors so the lambda below stays a plain (non-device)
   // lambda instead of a KOKKOS_LAMBDA, keeping this loop host-only.
   const auto LocComputeFrazilFormation = ComputeFrazilFormation;
   const auto LocComputeFrazilMelt      = ComputeFrazilMelt;

   I4 NClamped = 0;
   Kokkos::parallel_reduce(
       "frazilTeosImpl", Kokkos::RangePolicy<HostExecSpace>(0, NCellsAll),
       [=](const I4 ICell, I4 &Accum) {
          const I4 KMin = MinLayerCell(ICell);
          const I4 KMax = MaxLayerCell(ICell);

          I4 KLim          = KMax;
          bool HasKLim     = true;
          const bool Limit = (LocDepthLimit >= 0.0_Real);

          // calculates the depth limit based on geometric height
          if (Limit) {
             HasKLim = false;
             for (I4 K = KMax; K >= KMin; --K) {
                if (Kokkos::abs(LocGeomZMid(ICell, K)) <= LocDepthLimit) {
                   KLim    = K;
                   HasKLim = true;
                   break;
                }
             }
          }

          // Explicit accumulation order: bottom layer to top layer.
          for (I4 K = KMax; K >= KMin; --K) {
             if (!HasKLim || K > KLim) {
                LocFrazilHTend(ICell, K) = 0.0_Real;
                LocFrazilTTend(ICell, K) = 0.0_Real;
                LocFrazilSTend(ICell, K) = 0.0_Real;
                continue;
             }

             // sqrt(SA) in the freezing-point polynomial returns NaN for
             // SA < 0, which silently disables frazil since NaN compares false
             const Real AbsSalinityRaw = AbsSalinityH(ICell, K);
             if (AbsSalinityRaw < 0.0_Real) {
                ++Accum;
             }
             const Real AbsSalinityIn = Kokkos::max(0.0_Real, AbsSalinityRaw);
             const Real ConservTempIn = ConservTempH(ICell, K);
             const Real PressureIn    = PressureH(ICell, K);
             const Real PressureDb    = PressureIn * Pa2Db;
             const Real PseudoThicknessIn = PseudoThicknessH(ICell, K);

             const Real CtFreezing = Eos::calcCtFreezing(
                 LocEosChoice, AbsSalinityIn, PressureDb, 0.0_Real);

             Real HTend = 0.0_Real;
             Real TTend = 0.0_Real;
             Real STend = 0.0_Real;

             if (ConservTempIn < CtFreezing) {
                LocComputeFrazilFormation(
                    AbsSalinityIn, ConservTempIn, PressureDb, PseudoThicknessIn,
                    LocAccMIce(ICell), LocAccMLiq(ICell), LocAccMSalt(ICell),
                    LocAccELiq(ICell), LocAccEIce(ICell), HTend, TTend, STend);
             } else if (LocAccMIce(ICell) > 0.0_Real) {
                LocComputeFrazilMelt(AbsSalinityIn, ConservTempIn, PressureDb,
                                     PseudoThicknessIn, LocAccMIce(ICell),
                                     LocAccMLiq(ICell), LocAccMSalt(ICell),
                                     LocAccELiq(ICell), LocAccEIce(ICell),
                                     HTend, TTend, STend, CtFreezing);
             }

             // Per-call increments; FrazilOnCell normalizes these to rates.
             LocFrazilHTend(ICell, K) = HTend;
             LocFrazilTTend(ICell, K) = TTend;
             LocFrazilSTend(ICell, K) = STend;
          } // end of vertical loop

          // Convert to coupler units
          LocAccMIce(ICell)  = LocAccMIce(ICell) * RhoSw;
          LocAccMLiq(ICell)  = LocAccMLiq(ICell) * RhoSw;
          LocAccMSalt(ICell) = LocAccMSalt(ICell) * RhoSw * PPt2Salt;
          LocAccELiq(ICell)  = LocAccELiq(ICell) * RhoSw;
          LocAccEIce(ICell)  = LocAccEIce(ICell) * RhoSw;
       },
       NClamped); // end of NCells loop

   warnNegativeSalinity(NClamped);

   deepCopy(FrazilTTend, LocFrazilTTend);
   deepCopy(FrazilSTend, LocFrazilSTend);
   deepCopy(FrazilHTend, LocFrazilHTend);
   deepCopy(AccMIce, LocAccMIce);
   deepCopy(AccEIce, LocAccEIce);
   deepCopy(AccMLiq, LocAccMLiq);
   deepCopy(AccELiq, LocAccELiq);
   deepCopy(AccMSalt, LocAccMSalt);
}

void Frazil::computeFrazil(const Array2DReal &ConservTemp,
                           const Array2DReal &AbsSalinity,
                           const Array2DReal &Pressure,
                           const Array2DReal &PseudoThickness) {
   Eos *DefEos = Eos::getInstance();
   if (!DefEos) {
      ABORT_ERROR("Frazil::computeFrazil: Eos must be initialized before "
                  "computeFrazil");
   }

   switch (FrazilChoice) {
   case FrazilType::FixedPropertyFrazil:
      computeFrazilFixedPropertyImpl(ConservTemp, AbsSalinity, Pressure,
                                     PseudoThickness);
      break;
   case FrazilType::TeosFrazil:
      computeFrazilTeosImpl(ConservTemp, AbsSalinity, Pressure,
                            PseudoThickness);
      break;
   default:
      ABORT_ERROR("Frazil::computeFrazil: Unknown FrazilChoice");
      break;
   }

   if (ConservationCheck) {
      checkColumnConservation();
   }
} // end of computeFrazil

} // namespace OMEGA
