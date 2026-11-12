#ifndef OMEGA_FRAZIL_H
#define OMEGA_FRAZIL_H
//===-- ocn/Frazil.h - Frazil Ice Formation -------------------*- C++ -*-===//
//
// The Frazil class manages frazil tendencies and accumulators.
// This initial implementation only has a teos-10 configuration.
// but carries scaffolding for other implementations.
//
//===----------------------------------------------------------------------===//

#include "Config.h"
// #include "DataTypes.h"
#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"

#include <map>
#include <memory>
#include <string>

#include <gswteos-10.h>

namespace OMEGA {

enum class FrazilType {
   FixedPropertyFrazil, ///< fixed-property frazil option
   TeosFrazil           ///< TEOS frazil option
};

class FixedPropertyFrazilFormation {
 public:
   FixedPropertyFrazilFormation();

   Real LayerMassFracMax; ///< layer mass fraction limit (set in config)
   Real FrazilIceSalinity = IceRefSal; // Global constant
   Real LatFrazil         = LatIce;    // Global constant

   Real FrazilPorosity =
       1.0_Real; // Internal for now; can move to config later.

   KOKKOS_FUNCTION void operator()(const Real AbsSalinity,
                                   const Real ConservTemp,
                                   const Real PressureDb,
                                   const Real PseudoThickness,
                                   Real &SumIceThickness, Real &SumSalt,
                                   Real &SumEnergy, Real &HTend, Real &TTend,
                                   Real &STend, const Real CtFreezing) const {

      const Real Potential =
          PseudoThickness * Cp0Sw * RhoSw * (ConservTemp - CtFreezing);
      const Real FreezingEnergy = Kokkos::max(0.0_Real, -Potential);

      HTend = 0.0_Real;
      TTend = 0.0_Real;
      STend = 0.0_Real;

      Real NewFrzThickness =
          FreezingEnergy /
          (LatFrazil * RhoSw); // frazil (ice) mass in pseudo-thickness terms

      NewFrzThickness =
          Kokkos::min(NewFrzThickness, PseudoThickness * LayerMassFracMax);
      Real NewFrzEnergy =
          NewFrzThickness *
          (-LatFrazil +
           Cp0Sw * CtFreezing); // (<0; enthalpy of frazil, i.e. phase
                                // change and enthalpy of melted equ)

      // MANUAL TOGGLE: uncomment line below to use porosity
      // const Real FrazilIceSalinity = FrazilPorosity * AbsSalinity;

      const Real FrazilSalinity = Kokkos::min(FrazilIceSalinity, AbsSalinity);
      const Real NewSaltContent =
          NewFrzThickness * FrazilSalinity; // in m.(g/kg)

      // HTend does not include the salt mass contribution to align with mpas-o.
      HTend = -NewFrzThickness;
      // TTend includes the enthalpy associated with the mass flux to be
      // conservative with the melt impl. (this differs from old mpas-o)
      TTend =
          -(NewFrzEnergy) / (Cp0Sw); // (E< 0 so TTend>0) // scaled to h.CT tend
      STend = -NewSaltContent;

      SumIceThickness += NewFrzThickness;
      SumSalt += NewSaltContent;
      SumEnergy += NewFrzEnergy;
   }
};

class FixedPropertyFrazilMelt {
 public:
   FixedPropertyFrazilMelt();

   Real LayerMassFracMax;   ///< layer mass fraction limit (set in config)
   Real LatFrazil = LatIce; // Global constant

   KOKKOS_FUNCTION void operator()(const Real AbsSalinity,
                                   const Real ConservTemp,
                                   const Real PressureDb,
                                   const Real PseudoThickness,
                                   Real &SumIceThickness, Real &SumSalt,
                                   Real &SumEnergy, Real &HTend, Real &TTend,
                                   Real &STend, const Real CtFreezing) const {
      constexpr Real Eps = 1.0e-12_Real;

      if (SumIceThickness <=
          Eps) { // skipping melt if noise-level Ice is present
         HTend = 0.0_Real;
         TTend = 0.0_Real;
         STend = 0.0_Real;
         return;
      }
      const Real Potential =
          PseudoThickness * Cp0Sw * RhoSw * (ConservTemp - CtFreezing);
      const Real AvailableEnergy = Kokkos::max(0.0_Real, Potential);

      HTend = 0.0_Real;
      TTend = 0.0_Real;
      STend = 0.0_Real;

      Real MeltThickness =
          AvailableEnergy /
          (LatFrazil * RhoSw); // mass in pseudo-thickness units
      MeltThickness = Kokkos::min(MeltThickness, SumIceThickness);
      MeltThickness = Kokkos::min(
          MeltThickness,
          PseudoThickness * LayerMassFracMax); // also 0.1h lim on added mass
      const Real FrazilFractionMelted =
          MeltThickness / SumIceThickness; // mass fraction melted

      HTend = FrazilFractionMelted * (SumIceThickness); // (>0 so HTend>0)
      TTend = FrazilFractionMelted * SumEnergy /
              Cp0Sw; // (SumE <0 thus TTend < 0 when melting for phase change)
      STend = +FrazilFractionMelted * SumSalt; // (STend > 0 when melting)

      const Real FrazilFractionLeft = 1.0_Real - FrazilFractionMelted;
      SumIceThickness               = FrazilFractionLeft * SumIceThickness;
      SumSalt                       = FrazilFractionLeft * SumSalt;
      SumEnergy =
          FrazilFractionLeft * SumEnergy; // conservative by construction
   }
};

class FrazilMelt {
 public:
   /// constructor declaration
   FrazilMelt();

   // layer mass fraction limit parameter (set in config)
   Real LayerMassFracMax;

   //   The functor for FrazilMelt takes as inputs:
   //   the local ocean layer state (AbsSalinity, ConservTemp, PressureDb,
   //   PseudoThickness),
   //   the accumulated frazil solid and liquid mass, energy, and salt
   //   and outputs the frazil tendencies (HTend, TTend, STend) and updated
   //   accumulators.
   //   Host-only: relies on GSW TEOS-10 routines that are not device-callable.
   //   This is a temporary implementation until a device-callable solution is
   //   available.
   void operator()(const Real AbsSalinity, const Real ConservTemp,
                   const Real PressureDb, const Real PseudoThickness,
                   Real &AccMIce, Real &AccMLiq, Real &AccMSalt, Real &AccELiq,
                   Real &AccEIce, Real &HTend, Real &TTend, Real &STend,
                   const Real CtFreezing) const {

      constexpr Real Eps = 1.0e-12_Real;

      // this check on AccMIce and E should be done in the calling function, but
      // is here for safety we can do a better implementation of the checks

      if (AccMIce <= 0.0_Real || AccMLiq <= 0.0_Real ||
          AccEIce >= 0.0_Real) { // below calculations assume sign
         ABORT_ERROR("FrazilMelt: Invalid accumulator signs: AccMIce={}, "
                     "AccMLiq={}, AccELiq={}, AccEIce={}",
                     AccMIce, AccMLiq, AccELiq, AccEIce);
      }

      if (AccMIce <= Eps ||
          AccMLiq <= Eps) { // skipping melt calculation for noise-level AccMIce
         HTend = 0.0_Real;
         TTend = 0.0_Real;
         STend = 0.0_Real;
         return;
      }
      // 1. we start by adding the solid ice to the ocean layer (no brine yet)
      const Real PotEnthalpyIce = AccEIce / AccMIce;
      const Real NewLayerMass =
          PseudoThickness +
          AccMIce; // AccMIce > Eps and PseudoThickness passed ocean_validate()
      const Real NewLayerIceFraction = AccMIce / NewLayerMass;

      // 2. we calculate the (mass- and energy-conserving) ocean layer evolution
      // ... how much (pure) ice does this layer melt?

      // typecasting for now but will be simplified once gsw functions are
      // ported
      const double SA_d     = static_cast<double>(AbsSalinity);
      const double CT_d     = static_cast<double>(ConservTemp);
      const double P_d      = static_cast<double>(PressureDb);
      const double WIhIn_d  = static_cast<double>(NewLayerIceFraction);
      const double Pt0Ice_d = gsw_pt_from_pot_enthalpy_ice_poly(
          static_cast<double>(PotEnthalpyIce));
      const double TIce_d = gsw_t_from_pt0_ice(Pt0Ice_d, P_d);
      double SANew_d      = SA_d;
      double CTNew_d      = CT_d;
      double WIhOut_d     = WIhIn_d;

      gsw_melting_ice_into_seawater(SA_d, CT_d, P_d, WIhIn_d, TIce_d, &SANew_d,
                                    &CTNew_d, &WIhOut_d);

      // GSW marks all three of its failure exits with the same sentinel. The
      // ct < ctf exit is reachable in normal operation because the caller
      // gates melt on the polynomial freezing point while GSW uses the exact
      // one; at that point there is no melt energy available anyway.
      if (WIhOut_d > GSW_ERROR_LIMIT) {
         constexpr Real FreezingTol = 1.0e-3_Real;
         if (ConservTemp - CtFreezing < FreezingTol) {
            HTend = 0.0_Real;
            TTend = 0.0_Real;
            STend = 0.0_Real;
            return;
         }
         ABORT_ERROR("FrazilMelt: GSW returned invalid values for "
                     "AbsSalinity={}, ConservTemp={}, PressureDb={}, "
                     "PseudoThickness={}, CtFreezing={}",
                     AbsSalinity, ConservTemp, PressureDb, PseudoThickness,
                     CtFreezing);
      }

      const Real WIhOut = static_cast<Real>(
          WIhOut_d); // by def 0<= WIhOut <= 1 ; above checksfor invalid values

      // 3. we calculate the mass fraction of the frazil (pure) ice that was
      // melted - limited by a total mass limit of 0.1h
      const Real SolidMassMelted = Kokkos::max(
          0.0_Real,
          AccMIce - WIhOut * NewLayerMass); // original - left-over solid ice,
      const Real FrazilFractionMelted =
          Kokkos::min(SolidMassMelted / AccMIce,
                      PseudoThickness * LayerMassFracMax /
                          (AccMIce + AccMLiq)); // added mass < 0.1h

      // the frazil fraction based on the solid ice also sets the (proportional)
      // contributions from the frazil brine
      HTend = +(FrazilFractionMelted * (AccMLiq + AccMIce));
      TTend = +(FrazilFractionMelted * (AccELiq + AccEIce)) / Cp0Sw;
      STend = +(FrazilFractionMelted * AccMSalt);

      const Real FrazilFractionLeft = 1.0_Real - FrazilFractionMelted;

      AccMIce  = FrazilFractionLeft * AccMIce;
      AccMLiq  = FrazilFractionLeft * AccMLiq;
      AccMSalt = FrazilFractionLeft * AccMSalt;
      AccELiq  = FrazilFractionLeft * AccELiq;
      AccEIce  = FrazilFractionLeft * AccEIce;
   }
};

class FrazilFormation {
 public:
   /// constructor declaration
   FrazilFormation();

   /// Parameters for FrazilFormation (set by yaml file)
   Real Phi;              ///< liquid mass fraction of new frazil (0 < Phi < 1)
   Real LayerMassFracMax; ///< layer mass fraction limit for thickness tendency

   //   The functor for FrazilFormation takes as inputs:
   //   the local ocean layer state (AbsSalinity, ConservTemp, PressureDb,
   //   PseudoThickness),
   //   the accumulated frazil solid and liquid mass, energy, and salt
   //   and outputs the frazil tendencies (HTend, TTend, STend) and updated
   //   accumulators.
   //   Host-only: relies on GSW TEOS-10 routines that are not device-callable.
   //   This is a temporary implementation until a device-callable solution is
   //   available.
   void operator()(const Real AbsSalinity, const Real ConservTemp,
                   const Real PressureDb, const Real PseudoThickness,
                   Real &AccMIce, Real &AccMLiq, Real &AccMSalt, Real &AccELiq,
                   Real &AccEIce, Real &HTend, Real &TTend, Real &STend) const {

      Real ConservTempNew;
      Real AbsSalinityNew;
      Real WIh            = 0.0_Real;
      Real SolidMass      = 0.0_Real;
      Real LiquidMass     = 0.0_Real;
      Real SolidEnthalpy  = 0.0_Real;
      Real LiquidEnthalpy = 0.0_Real;

      // type casting for now because gsw expects double
      // and Real can be either float or double depending on build configuration
      double SANew_d = 0.0;
      double CTNew_d = 0.0;
      double WIh_d   = 0.0;
      // double PTNew_d = 0.0;

      gsw_frazil_properties_potential_poly(
          static_cast<double>(AbsSalinity),
          static_cast<double>(Cp0Sw * ConservTemp),
          static_cast<double>(PressureDb), &SANew_d, &CTNew_d, &WIh_d);

      // GSW flags out-of-domain input (w_Ih > 0.9) by setting all three outputs
      // to GSW_INVALID_VALUE; the LayerMassFracMax clamp below would otherwise
      // hide it.
      if (WIh_d > GSW_ERROR_LIMIT) {
         ABORT_ERROR("FrazilFormation: GSW returned invalid values for "
                     "AbsSalinity={}, ConservTemp={}, PressureDb={}, "
                     "PseudoThickness={}",
                     AbsSalinity, ConservTemp, PressureDb, PseudoThickness);
      }

      double PTNew_d =
          gsw_pt_from_ct(SANew_d, CTNew_d); // convert to potential temperature

      AbsSalinityNew = static_cast<Real>(SANew_d);
      ConservTempNew = static_cast<Real>(CTNew_d);
      WIh            = static_cast<Real>(WIh_d);

      const Real OneMinusPhi = Kokkos::max(1.0e-12_Real, 1.0_Real - Phi);
      // anything called mass below is in pseudo-thickness units (m) and needs
      // to be scaled by RhoSw for coupling
      SolidMass =
          PseudoThickness * Kokkos::min(WIh, OneMinusPhi * LayerMassFracMax);
      LiquidMass     = (Phi / OneMinusPhi) * SolidMass;
      SolidEnthalpy  = SolidMass * gsw_pot_enthalpy_from_pt_ice_poly(PTNew_d);
      LiquidEnthalpy = LiquidMass * Cp0Sw * ConservTempNew;
      // per timestep (not scaled by dt here)
      HTend = -(SolidMass +
                LiquidMass); // because Phi is a *mass* fraction, LiquidMass
                             // includes the salt contribution to mass.
      TTend = -(LiquidEnthalpy + SolidEnthalpy) / Cp0Sw;
      STend = -(LiquidMass * AbsSalinityNew);

      // Local unit of mass is pseudo thickness (m)
      // these all need a RhoSw factor before coupling
      AccMIce += SolidMass;                    // m
      AccMLiq += LiquidMass;                   // m
      AccMSalt += LiquidMass * AbsSalinityNew; // (m)(g/kg)
      AccELiq += LiquidEnthalpy;               // (m)(J/kg)
      AccEIce += SolidEnthalpy;                // (m)(J/kg)
   }
};

class Frazil {
 public:
   static void init();
   /// Creates a new frazil object and stores it in the AllFrazil map.
   static Frazil *create(const std::string &Name);

   /// Retrieve frazil object by name.
   static Frazil *get(const std::string &Name);

   /// Retrieve default frazil object.
   static Frazil *getDefault();

   /// Destructor
   ~Frazil();

   /// Deallocates arrays
   static void clear();

   /// Remove frazil object by name.
   static void erase(std::string InName); ///< [in] name to remove

   Array2DReal FrazilTTend;
   Array2DReal FrazilSTend;
   Array2DReal FrazilHTend;
   Array1DReal AccMIce;
   Array1DReal AccEIce;
   Array1DReal AccMLiq;
   Array1DReal AccELiq;
   Array1DReal AccMSalt;
   Array1DReal FrazilMassFlux;
   Array1DReal FrazilSaltFlux;
   Array1DReal FrazilEnergyFlux;

   void computeFrazil(const Array2DReal &ConservTemp,
                      const Array2DReal &AbsSalinity,
                      const Array2DReal &Pressure,
                      const Array2DReal &PseudoThickness);
   void resetOcnStepFluxes();
   void accumulateOcnStepFluxes(Real FinalUpdateWeight, R8 TimeStepSeconds);
   void registerFields();
   void unregisterFields();
   void computeFrazilFixedPropertyImpl(const Array2DReal &ConservTemp,
                                       const Array2DReal &AbsSalinity,
                                       const Array2DReal &Pressure,
                                       const Array2DReal &PseudoThickness);
   void computeFrazilTeosImpl(const Array2DReal &ConservTemp,
                              const Array2DReal &AbsSalinity,
                              const Array2DReal &Pressure,
                              const Array2DReal &PseudoThickness);
   bool ConservationCheck = false;
   Real DepthLimit        = -1.0_Real;

 private:
   static Frazil *DefaultFrazil;
   static std::map<std::string, std::unique_ptr<Frazil>> AllFrazil;

   Frazil(const HorzMesh *Mesh, const VertCoord *VCoord);

   // Forbid copy and move construction/assignment.
   Frazil(const Frazil &)            = delete;
   Frazil &operator=(const Frazil &) = delete;
   Frazil(Frazil &&)                 = delete;
   Frazil &operator=(Frazil &&)      = delete;

   FrazilType FrazilChoice;
   FixedPropertyFrazilFormation ComputeFixedPropertyFrazilFormation;
   FixedPropertyFrazilMelt ComputeFixedPropertyFrazilMelt;
   FrazilFormation ComputeFrazilFormation;
   FrazilMelt ComputeFrazilMelt;
   I4 NCellsAll;
   I4 NChunks;

   const HorzMesh *MeshPtr;
   const VertCoord *VCoordPtr;
   bool FieldsRegistered       = false;
   bool WarnedNegativeSalinity = false;

   void warnNegativeSalinity(I4 NClamped);

   void checkColumnConservation() const;
};

} // namespace OMEGA

#endif
