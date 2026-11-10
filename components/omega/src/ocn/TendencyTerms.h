#ifndef OMEGA_TENDENCYTERMS_H
#define OMEGA_TENDENCYTERMS_H
//===-- ocn/TendencyTerms.h - Tendency Terms --------------------*- C++ -*-===//
//
/// \file
/// \brief Contains functors for calculating tendency terms
///
/// This header defines functors to be called by the time-stepping scheme
/// to calculate tendencies used to update state variables.
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "Eos.h"
#include "Frazil.h"
#include "GlobalConstants.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "SurfaceFlux.h"
#include "VertCoord.h"

#include <cmath> // for std::copysign
#include <iomanip>

namespace OMEGA {

/// Divergence of pseudo-thickness flux at cell centers, for updating
/// pseudo-thickness arrays
class PseudoThicknessFluxDivOnCell {
 public:
   bool Enabled = false;

   /// constructor declaration
   PseudoThicknessFluxDivOnCell(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes cell index, vertical chunk index, and pseudo-thickness
   /// flux array as inputs, outputs the tendency array
   KOKKOS_FUNCTION void operator()(const TeamMember &Team,
                                   const Array2DReal &Tend, I4 ICell,
                                   const Array2DReal &PseudoThicknessFlux,
                                   const Array2DReal &NormalVelEdge) const {

      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      ScratchArray1DReal DivTmp(teamScratch(Team), NVertLayers);
      parallelForInner(
          Team, NVertLayers, INNER_LAMBDA(int K) { DivTmp(K) = 0; });

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge = EdgesOnCell(ICell, J);

         const int MinLyrEdgeBot = MinLayerEdgeBot(JEdge);
         const int MaxLyrEdgeTop = MaxLayerEdgeTop(JEdge);

         parallelForInner(
             Team, Range{MinLyrEdgeBot, MaxLyrEdgeTop}, INNER_LAMBDA(int K) {
                DivTmp(K) -= DvEdge(JEdge) * EdgeSignOnCell(ICell, J) *
                             PseudoThicknessFlux(JEdge, K) *
                             NormalVelEdge(JEdge, K) * InvAreaCell;
             });
      }

      const int MinLyrCell = MinLayerCell(ICell);
      const int MaxLyrCell = MaxLayerCell(ICell);

      parallelForInner(
          Team, Range{MinLyrCell, MaxLyrCell},
          INNER_LAMBDA(int K) { Tend(ICell, K) -= DivTmp(K); });
   }

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array1DReal DvEdge;
   Array1DReal AreaCell;
   Array2DReal EdgeSignOnCell;
   I4 NVertLayers;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Horizontal advection of potential vorticity defined on edges, for
/// momentum equation
class PotentialVortHAdvOnEdge {
 public:
   bool Enabled = false;

   /// constructor declaration
   PotentialVortHAdvOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes edge index, vertical chunk index, and arrays for
   /// normalized relative vorticity, normalized planetary vorticity, layer
   /// thickness on edges, and normal velocity on edges as inputs,
   /// outputs the tendency array
   KOKKOS_FUNCTION void operator()(const TeamMember &Team,
                                   const Array2DReal &Tend, I4 IEdge,
                                   const Array2DReal &NormRVortEdge,
                                   const Array2DReal &NormFEdge,
                                   const Array2DReal &FluxPseudoThickEdge,
                                   const Array2DReal &NormVelEdge) const {

      ScratchArray1DReal VortTmp(teamScratch(Team), NVertLayers);

      const int KMin = MinLayerEdgeBot(IEdge);
      const int KMax = MaxLayerEdgeTop(IEdge);

      parallelForInner(
          Team, NVertLayers, INNER_LAMBDA(int K) { VortTmp(K) = 0; });

      for (int J = 0; J < NEdgesOnEdge(IEdge); ++J) {
         I4 JEdge = EdgesOnEdge(IEdge, J);
         parallelForInner(
             Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                Real NormVort =
                    (NormRVortEdge(IEdge, K) + NormFEdge(IEdge, K) +
                     NormRVortEdge(JEdge, K) + NormFEdge(JEdge, K)) *
                    0.5_Real;

                VortTmp(K) += WeightsOnEdge(IEdge, J) *
                              FluxPseudoThickEdge(JEdge, K) *
                              NormVelEdge(JEdge, K) * NormVort;
             });
      }

      parallelForInner(
          Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
             Tend(IEdge, K) += EdgeMask(IEdge, K) * VortTmp(K);
          });
   }

   /// Relative vorticity horizontal advection without Coriolis. Used for the
   /// split-explicit baroclinic velocity forcing term.
   KOKKOS_FUNCTION void operator()(const TeamMember &Team,
                                   const Array2DReal &Tend, I4 IEdge,
                                   const Array2DReal &NormRVortEdge,
                                   const Array2DReal &FluxPseudoThickEdge,
                                   const Array2DReal &NormVelEdge) const {

      ScratchArray1DReal VortTmp(teamScratch(Team), NVertLayers);

      const int KMin = MinLayerEdgeBot(IEdge);
      const int KMax = MaxLayerEdgeTop(IEdge);

      parallelForInner(
          Team, NVertLayers, INNER_LAMBDA(int K) { VortTmp(K) = 0; });

      for (int J = 0; J < NEdgesOnEdge(IEdge); ++J) {
         I4 JEdge = EdgesOnEdge(IEdge, J);
         parallelForInner(
             Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                const Real NormVort =
                    (NormRVortEdge(IEdge, K) + NormRVortEdge(JEdge, K)) *
                    0.5_Real;

                VortTmp(K) += WeightsOnEdge(IEdge, J) *
                              FluxPseudoThickEdge(JEdge, K) *
                              NormVelEdge(JEdge, K) * NormVort;
             });
      }

      parallelForInner(
          Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
             Tend(IEdge, K) += EdgeMask(IEdge, K) * VortTmp(K);
          });
   }

 private:
   Array1DI4 NEdgesOnEdge;
   Array2DI4 EdgesOnEdge;
   Array2DReal WeightsOnEdge;
   I4 NVertLayers;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Coriolis acceleration on edges, f times tangential velocity reconstruction
class CoriolisAccelerationOnEdge {
 public:
   /// constructor declaration
   CoriolisAccelerationOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes thread team, edge index, base tendency, velocity on
   /// edges, and Coriolis parameter on edges as inputs, updates the tendency
   /// array with base tendency plus Coriolis acceleration
   KOKKOS_FUNCTION void operator()(const TeamMember &Team,
                                   const Array2DReal &Tend,
                                   const Array2DReal &BaseTend, I4 IEdge,
                                   const Array2DReal &NormalVelEdge,
                                   const Array1DReal &FEdge) const {

      ScratchArray1DReal AccelTmp(teamScratch(Team), NVertLayers);

      const int KMin = MinLayerEdgeBot(IEdge);
      const int KMax = MaxLayerEdgeTop(IEdge);

      parallelForInner(
          Team, NVertLayers, INNER_LAMBDA(int K) { AccelTmp(K) = 0; });

      for (int J = 0; J < NEdgesOnEdge(IEdge); ++J) {
         const I4 JEdge = EdgesOnEdge(IEdge, J);
         parallelForInner(
             Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                AccelTmp(K) += WeightsOnEdge(IEdge, J) *
                               NormalVelEdge(JEdge, K) * FEdge(JEdge);
             });
      }

      parallelForInner(
          Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
             Tend(IEdge, K) = BaseTend(IEdge, K) + AccelTmp(K);
          });
   }

   /// The functor takes edge index, barotropic velocity on edges, and Coriolis
   /// parameter on edges as inputs, updates the barotropic tendency array
   KOKKOS_FUNCTION void operator()(const Array1DReal &Tend, I4 IEdge,
                                   const Array1DReal &NormalVelEdge,
                                   const Array1DReal &FEdge) const {

      Real AccelTmp = 0._Real;

      for (int J = 0; J < NEdgesOnEdge(IEdge); ++J) {
         const I4 JEdge = EdgesOnEdge(IEdge, J);
         AccelTmp +=
             WeightsOnEdge(IEdge, J) * NormalVelEdge(JEdge) * FEdge(JEdge);
      }

      Tend(IEdge) += AccelTmp;
   }

 private:
   Array1DI4 NEdgesOnEdge;
   Array2DI4 EdgesOnEdge;
   Array2DReal WeightsOnEdge;
   I4 NVertLayers;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Gradient of kinetic energy defined on edges, for momentum equation
class KEGradOnEdge {
 public:
   bool Enabled = false;

   /// constructor declaration
   KEGradOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes edge index, vertical chunk index, and kinetic energy
   /// array as inputs, outputs the tendency array
   KOKKOS_FUNCTION void operator()(const TeamMember &Team,
                                   const Array2DReal &Tend, I4 IEdge,
                                   const Array2DReal &KECell) const {

      const I4 JCell0 = CellsOnEdge(IEdge, 0);
      const I4 JCell1 = CellsOnEdge(IEdge, 1);

      const Real InvDcEdge = 1._Real / DcEdge(IEdge);

      const int KMin = MinLayerEdgeBot(IEdge);
      const int KMax = MaxLayerEdgeTop(IEdge);

      parallelForInner(
          Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
             Tend(IEdge, K) -= EdgeMask(IEdge, K) *
                               (KECell(JCell1, K) - KECell(JCell0, K)) *
                               InvDcEdge;
          });
   }

 private:
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Gradient of sea surface height defined on edges multipled by gravitational
/// acceleration, for momentum equation
class SSHGradOnEdge {
 public:
   bool Enabled = false;

   /// constructor declaration
   SSHGradOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes edge index, vertical chunk index, and array of
   /// pseudo-thickness/SSH, outputs tendency array
   KOKKOS_FUNCTION void operator()(const TeamMember &Team,
                                   const Array2DReal &Tend, I4 IEdge,
                                   const Array1DReal &SshCell) const {

      const I4 ICell0 = CellsOnEdge(IEdge, 0);
      const I4 ICell1 = CellsOnEdge(IEdge, 1);

      const Real InvDcEdge = 1._Real / DcEdge(IEdge);

      const int KMin = MinLayerEdgeBot(IEdge);
      const int KMax = MaxLayerEdgeTop(IEdge);

      parallelForInner(
          Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
             Tend(IEdge, K) -= EdgeMask(IEdge, K) * Gravity *
                               (SshCell(ICell1) - SshCell(ICell0)) * InvDcEdge;
          });
   }

 private:
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Laplacian horizontal mixing, for momentum equation
class VelocityDiffusionOnEdge {
 public:
   bool Enabled = false;

   Real ViscDel2;

   /// constructor declaration
   VelocityDiffusionOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes edge index, vertical chunk index, and arrays for
   /// divergence of horizontal velocity (defined at cell centers) and relative
   /// vorticity (defined at vertices), outputs tendency array
   KOKKOS_FUNCTION void operator()(const TeamMember &Team,
                                   const Array2DReal &Tend, I4 IEdge,
                                   const Array2DReal &DivCell,
                                   const Array2DReal &RVortVertex) const {

      const I4 ICell0 = CellsOnEdge(IEdge, 0);
      const I4 ICell1 = CellsOnEdge(IEdge, 1);

      const I4 IVertex0 = VerticesOnEdge(IEdge, 0);
      const I4 IVertex1 = VerticesOnEdge(IEdge, 1);

      const Real DcEdgeInv = 1._Real / DcEdge(IEdge);
      const Real DvEdgeInv = 1._Real / DvEdge(IEdge);

      const int KMin = MinLayerEdgeBot(IEdge);
      const int KMax = MaxLayerEdgeTop(IEdge);

      parallelForInner(
          Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
             const Real Del2U =
                 ((DivCell(ICell1, K) - DivCell(ICell0, K)) * DcEdgeInv -
                  (RVortVertex(IVertex1, K) - RVortVertex(IVertex0, K)) *
                      DvEdgeInv);

             Tend(IEdge, K) +=
                 EdgeMask(IEdge, K) * ViscDel2 * MeshScalingDel2(IEdge) * Del2U;
          });
   }

 private:
   Array2DI4 CellsOnEdge;
   Array2DI4 VerticesOnEdge;
   Array1DReal DcEdge;
   Array1DReal DvEdge;
   Array1DReal MeshScalingDel2;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Biharmonic horizontal mixing, for momentum equation
class VelocityHyperDiffOnEdge {
 public:
   bool Enabled = false;

   Real ViscDel4;
   Real DivFactor;

   /// Constructor declaration
   VelocityHyperDiffOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes the edge index, vertical chunk index, and arrays for
   /// the laplacian of divergence of horizontal velocity and the laplacian of
   /// the relative vorticity, outputs tendency array
   KOKKOS_FUNCTION void operator()(const TeamMember &Team,
                                   const Array2DReal &Tend, I4 IEdge,
                                   const Array2DReal &Del2DivCell,
                                   const Array2DReal &Del2RVortVertex) const {

      const I4 ICell0 = CellsOnEdge(IEdge, 0);
      const I4 ICell1 = CellsOnEdge(IEdge, 1);

      const I4 IVertex0 = VerticesOnEdge(IEdge, 0);
      const I4 IVertex1 = VerticesOnEdge(IEdge, 1);

      const Real DcEdgeInv = 1._Real / DcEdge(IEdge);
      const Real DvEdgeInv = 1._Real / DvEdge(IEdge);

      const int KMin = MinLayerEdgeBot(IEdge);
      const int KMax = MaxLayerEdgeTop(IEdge);

      parallelForInner(
          Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
             const Real Del2U =
                 (DivFactor *
                      (Del2DivCell(ICell1, K) - Del2DivCell(ICell0, K)) *
                      DcEdgeInv -
                  (Del2RVortVertex(IVertex1, K) -
                   Del2RVortVertex(IVertex0, K)) *
                      DvEdgeInv);

             Tend(IEdge, K) -=
                 EdgeMask(IEdge, K) * ViscDel4 * MeshScalingDel4(IEdge) * Del2U;
          });
   }

 private:
   Array2DI4 CellsOnEdge;
   Array2DI4 VerticesOnEdge;
   Array1DReal DcEdge;
   Array1DReal DvEdge;
   Array1DReal MeshScalingDel4;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Surface stress forcing (eg. wind)
class SfcStressForcingOnEdge {
 public:
   bool Enabled = false;

   /// constructor declaration
   SfcStressForcingOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes the edge index, vertical chunk index, and arrays for
   /// normal surface stress and edge pseudo-thickness, outputs tendency array
   KOKKOS_FUNCTION void operator()(const TeamMember &Team,
                                   const Array2DReal &Tend, I4 IEdge,
                                   const Array1DReal &NormalStressEdge,
                                   const Array2DReal &PseudoThickEdge) const {
      const I4 KMin = MinLayerEdgeBot(IEdge);
      const I4 KMax = MaxLayerEdgeTop(IEdge);

      parallelForInner(
          Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
             if (K == KMin) {
                const Real InvThickEdge = 1._Real / PseudoThickEdge(IEdge, K);
                Tend(IEdge, K) += EdgeMask(IEdge, K) * InvThickEdge *
                                  NormalStressEdge(IEdge) / RhoSw;
             }
          });
   }

 private:
   Array2DReal EdgeMask;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// Bottom drag
class BottomDragOnEdge {
 public:
   bool Enabled = false;
   Real Coeff;

   /// constructor declaration
   BottomDragOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord);

   /// The functor takes the edge index and arrays for
   /// horizontal velocity, kinetic energy,
   /// and edge pseudo-thickness, outputs tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge,
                                   const Array2DReal &NormalVelEdge,
                                   const Array2DReal &KECell,
                                   const Array2DReal &PseudoThickEdge) const {
      const I4 KBot = MaxLayerEdgeTop(IEdge);

      // Land edges and the outermost edges of the halo have no active layer
      // on both sides, and MaxLayerEdgeTop is set to -1 for them. Unlike the
      // other tendency terms, which are chunked over [MinLayerEdgeBot,
      // MaxLayerEdgeTop] and so skip such edges automatically, bottom drag
      // indexes the bottom layer directly and must exclude them explicitly.
      if (KBot < 0)
         return;

      const I4 JCell0 = CellsOnEdge(IEdge, 0);
      const I4 JCell1 = CellsOnEdge(IEdge, 1);

      const Real VelNormEdge =
          Kokkos::sqrt(KECell(JCell0, KBot) + KECell(JCell1, KBot));

      const Real InvThickEdge = 1._Real / PseudoThickEdge(IEdge, KBot);
      Tend(IEdge, KBot) -= EdgeMask(IEdge, KBot) * Coeff * VelNormEdge *
                           InvThickEdge * NormalVelEdge(IEdge, KBot);
   }

 private:
   I4 NVertLayers;
   Array2DI4 CellsOnEdge;
   Array2DReal EdgeMask;
   Array1DI4 MaxLayerEdgeTop;
};

/// Coupled freshwater flux forcing for thickness equation.
class SfcThicknessForcingOnCell {
 public:
   bool Enabled = false;

   SfcThicknessForcingOnCell(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 ICell,
                                   const Array1DReal &SnowFlux,
                                   const Array1DReal &RainFlux,
                                   const Array1DReal &EvaporationFlux,
                                   const Array1DReal &SeaIceFreshWaterFlux,
                                   const Array1DReal &IceRunoffFlux,
                                   const Array1DReal &RiverRunoffFlux,
                                   const Array1DReal &SeaIceSaltFlux) const {

      const I4 KTop = MinLayerCell(ICell);
      if (KTop > MaxLayerCell(ICell)) {
         return;
      }

      const Real FreshWaterFlux = sfcFreshWaterFlux(
          ICell, SnowFlux, RainFlux, EvaporationFlux, SeaIceFreshWaterFlux,
          IceRunoffFlux, RiverRunoffFlux);

      Tend(ICell, KTop) += (FreshWaterFlux + SeaIceSaltFlux(ICell)) / RhoSw;
   }

 private:
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

/// Coupled surface flux forcing for active tracers.
class SfcTracerForcingOnCell {
 public:
   bool Enabled                  = false;
   bool IncludeShortWaveHeatFlux = true;

   SfcTracerForcingOnCell(const HorzMesh *Mesh, const VertCoord *VCoord,
                          I4 TempTracerIndex, I4 SaltTracerIndex,
                          const Eos *EosInst);

   KOKKOS_FUNCTION void operator()(
       const Array3DReal &Tend, const Array2DReal &SurfaceTracerFlux,
       const Array1DReal &ExtraEnergy, I4 ICell, const Array3DReal &TracerCell,
       const Array1DReal &LatentHeatFluxEvap,
       const Array1DReal &SensibleHeatFlux,
       const Array1DReal &LongWaveHeatFluxUp,
       const Array1DReal &LongWaveHeatFluxDown,
       const Array1DReal &SeaIceHeatFlux, const Array1DReal &ShortWaveHeatFlux,
       const Array1DReal &SnowFlux, const Array1DReal &RainFlux,
       const Array1DReal &IceRunoffFlux, const Array1DReal &RiverRunoffFlux,
       const Array1DReal &EvaporationFlux,
       const Array1DReal &SeaIceSaltFlux) const {

      const I4 KTop = MinLayerCell(ICell);
      if (KTop > MaxLayerCell(ICell)) {
         return;
      }

      if (TempIndex >= 0) {
         const Real CtTop = TracerCell(TempIndex, ICell, KTop);

         // CT tendencies are due to direct heat fluxes + pot enthalpy fluxes
         // Each mass flux has an associated potential enthalpy flux.
         // Levels of simplification can be done here. For now:
         // - We approximate PotEnthalpyIce(Tinsitu, P=0) ~ -LatIce (constant);
         // Altrntively, we could use cnst PotEnthalpyIce(0, 0) ​= −333360
         // J/kg a 0.1% / 340 J/kg difference with the LatIce value from pcd.
         // The full expression is gsw_pot_enthalpy_ice(T, P).
         const Real PotEnthalpyIce = -LatIce;
         // - We assume dry snow and use the same PotEnthalpyIce for snow.
         // - We assume liquid water comes in at the specific enthalpy as the
         // top ocean layer, (i.e. mass flux is CT-neutral). The enthalpy is
         // capped by a lower bound of pot enthalpy of freshwater at 0.0 C.
         // Another choice would be to add it at the same in situ temp as ocean
         // i.e. CT(Sa=0, max(0, T)).
         //- We assume evaporation removes the same specific enthalpy as the
         // top ocean layer; not capped, to keep the mass flux CT-neutral.
         const Real CtLim =
             (EosChoice == EosType::Teos10Eos) ? Ct0Fw : 0.0_Real;
         const Real PotEnthalpyFwIn  = Cp0Sw * Kokkos::max(CtLim, CtTop);
         const Real PotEnthalpyFwOut = Cp0Sw * CtTop;
         const Real ShortWaveFlux =
             IncludeShortWaveHeatFlux ? ShortWaveHeatFlux(ICell) : 0.0_Real;
         const Real HeatFlux =
             LongWaveHeatFluxUp(ICell) + LongWaveHeatFluxDown(ICell) +
             ShortWaveFlux + SensibleHeatFlux(ICell) +
             SeaIceHeatFlux(ICell) + // includes enthalpy of meltwater already
             (RainFlux(ICell) + RiverRunoffFlux(ICell)) * PotEnthalpyFwIn +
             LatentHeatFluxEvap(ICell) +
             EvaporationFlux(ICell) * PotEnthalpyFwOut +
             (SnowFlux(ICell) + IceRunoffFlux(ICell)) * PotEnthalpyIce;

         const Real TempFlux = HeatFlux * HFluxFac;
         Tend(TempIndex, ICell, KTop) += TempFlux;
         SurfaceTracerFlux(TempIndex, ICell) = TempFlux;
         // Accounting for extra energy/enthalpy (to be taken from atm)
         // - All direct heat fluxes are accounted for (via coupler)
         // - SeaIceHeatFlux comes entirely from the sea ice
         // - Phase Change for evaporation and frozen should be taken out of
         // ocean Extra terms are the enthalpy of liq water wrt. mass fluxes
         // Formulation allows for non-zero term if PotEnthalpyIce is updated
         ExtraEnergy(ICell) =
             (RainFlux(ICell) + RiverRunoffFlux(ICell)) * PotEnthalpyFwIn +
             EvaporationFlux(ICell) * PotEnthalpyFwOut +
             (SnowFlux(ICell) + IceRunoffFlux(ICell)) *
                 (PotEnthalpyIce + LatIce);
      }

      if (SaltIndex >= 0) {
         const Real SaltFlux = SeaIceSaltFlux(ICell) * SFluxFac;
         Tend(SaltIndex, ICell, KTop) += SaltFlux;
         SurfaceTracerFlux(SaltIndex, ICell) = SaltFlux;
      }
   }

 private:
   I4 TempIndex;
   I4 SaltIndex;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   EosType EosChoice;
};

// Tracer horizontal advection term
class TracerHorzAdvOnCell {
 public:
   bool Enabled           = true;
   bool ForceLowOrder     = false;
   bool FCT               = false;
   bool ComputeBudgets    = false;
   bool MonotonicityCheck = false;
   // coefficient for blending high-order terms
   Real Coef3rdOrder = 0.25;
   TracerHorzAdvOnCell(const HorzMesh *Mesh, const VertCoord *VCoord,
                       const VertAdv *VAdv);
   KOKKOS_FUNCTION ~TracerHorzAdvOnCell() {}
   void init();
   KOKKOS_FUNCTION void operator()(const TeamMember &Team, const I4 L,
                                   const I4 IEdge,
                                   const Array3DReal &TracerCell,
                                   const Array2DReal &FluxPseudoThickEdge,
                                   const Array2DReal &NormVelEdge) const {

      ScratchArray1DReal FlxTmp(teamScratch(Team), NVertLayers);

      const auto LTracerCell =
          subviewUnmanaged(TracerCell, L, Kokkos::ALL, Kokkos::ALL);

      const int KMin = MinLayerEdgeBot(IEdge);
      const int KMax = MaxLayerEdgeTop(IEdge);

      parallelForInner(
          Team, NVertLayers, INNER_LAMBDA(int K) { FlxTmp(K) = 0; });

      // Stay at low order at boundaries
      parallelForInner(
          Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
             const I4 JCell0 = CellsOnEdge(IEdge, 0);
             const I4 JCell1 = CellsOnEdge(IEdge, 1);
             const Real NormalThicknessFlux =
                 FluxPseudoThickEdge(IEdge, K) * NormVelEdge(IEdge, K);
             const Real TracerWgt =
                 DvEdge(IEdge) * 0.5_Real * NormalThicknessFlux;
             FlxTmp(K) += TracerWgt * (1._Real - AdvMaskHighOrder(IEdge, K)) *
                          (LTracerCell(JCell1, K) + LTracerCell(JCell0, K));
          });

      // High order (3rd or 4th) fluxes elsewhere when requested
      //    - If HorzTracerFluxOrder = 2, NAdvCellsForEdge = 0 and
      //      this loop is skipped.
      for (int I = 0; I < NAdvCellsForEdge(IEdge); ++I) {
         const I4 ICell = AdvCellsForEdge(IEdge, I);
         parallelForInner(
             Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                const Real NormalThicknessFlux =
                    FluxPseudoThickEdge(IEdge, K) * NormVelEdge(IEdge, K);
                const Real TracerWgt =
                    (AdvCoefs(I, IEdge) +
                     Coef3rdOrder *
                         std::copysign(1._Real, NormalThicknessFlux) *
                         AdvCoefs3rd(I, IEdge)) *
                    NormalThicknessFlux;
                FlxTmp(K) += TracerWgt * LTracerCell(ICell, K) *
                             AdvMaskHighOrder(IEdge, K);
             });
      }

      const auto LHighOrderFlxHorz =
          subviewUnmanaged(HighOrderFlxHorz, L, Kokkos::ALL, Kokkos::ALL);

      parallelForInner(
          Team, Range{KMin, KMax},
          INNER_LAMBDA(int K) { LHighOrderFlxHorz(IEdge, K) = FlxTmp(K); });
   }

   KOKKOS_FUNCTION void operator()(const TeamMember &Team,
                                   const Array3DReal &Tend, const I4 L,
                                   const I4 ICell) const {

      const auto LTend = subviewUnmanaged(Tend, L, Kokkos::ALL, Kokkos::ALL);
      const auto LHighOrderFlxHorz =
          subviewUnmanaged(HighOrderFlxHorz, L, Kokkos::ALL, Kokkos::ALL);

      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      ScratchArray1DReal TendTmp(teamScratch(Team), NVertLayers);
      parallelForInner(
          Team, NVertLayers, INNER_LAMBDA(int K) { TendTmp(K) = 0; });

      for (int I = 0; I < NEdgesOnCell(ICell); ++I) {
         const I4 IEdge = EdgesOnCell(ICell, I);

         const int MinLyrEdgeBot = MinLayerEdgeBot(IEdge);
         const int MaxLyrEdgeTop = MaxLayerEdgeTop(IEdge);

         parallelForInner(
             Team, Range{MinLyrEdgeBot, MaxLyrEdgeTop}, INNER_LAMBDA(int K) {
                TendTmp(K) += EdgeSignOnCell(ICell, I) *
                              LHighOrderFlxHorz(IEdge, K) * InvAreaCell;
             });
      }

      const int MinLyrCell = MinLayerCell(ICell);
      const int MaxLyrCell = MaxLayerCell(ICell);

      parallelForInner(
          Team, Range{MinLyrCell, MaxLyrCell},
          INNER_LAMBDA(int K) { LTend(ICell, K) += TendTmp(K); });
   }

   KOKKOS_FUNCTION void
   FCTProvisionaLayerThicknesses(const TeamMember &Team, const I4 ICell,
                                 const Real Dt,
                                 const Array2DReal &FluxPseudoThickEdge,
                                 const Array2DReal &LayerThickness,
                                 const Array2DReal &NormVelEdge) const {
      const int KMin = MinLayerCell(ICell);
      const int KMax = MaxLayerCell(ICell);
      parallelForInner(
          Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
             FCTProvisionaLayerThicknesses(ICell, K, Dt, FluxPseudoThickEdge,
                                           LayerThickness, NormVelEdge);
          });
   }

   KOKKOS_FUNCTION void
   FCTProvisionaLayerThicknesses(const I4 ICell, const I4 K, const Real Dt,
                                 const Array2DReal &FluxPseudoThickEdge,
                                 const Array2DReal &LayerThickness,
                                 const Array2DReal &NormVelEdge) const {

      const Real InvAreaCell = Dt / AreaCell(ICell);
      HProv(ICell, K)        = LayerThickness(ICell, K);
      for (I4 I = 0; I < NEdgesOnCell(ICell); ++I) {
         const I4 IEdge = EdgesOnCell(ICell, I);
         const Real SignedFactor =
             InvAreaCell * DvEdge(IEdge) * EdgeSignOnCell(ICell, I);
         // Provisional layer thickness is after horizontal
         // thickness flux only
         const Real NormalThicknessFlux =
             FluxPseudoThickEdge(IEdge, K) * NormVelEdge(IEdge, K);
         HProv(ICell, K) += SignedFactor * NormalThicknessFlux;
      }
      // New layer thickness is after horizontal and vertical
      // thickness flux
      HProvInv(ICell, K) = 1.0_Real / HProv(ICell, K);
      HNewInv(ICell, K) =
          1.0_Real / (HProv(ICell, K) -
                      Dt * TotalVerticalTransportPseudoVelocity(ICell, K) +
                      Dt * TotalVerticalTransportPseudoVelocity(ICell, K + 1));
   }

   KOKKOS_FUNCTION void FCTTracerCurFill(const TeamMember &Team, const I4 L,
                                         const I4 ICell,
                                         const Array3DReal &TracerArray) const {
      const int KMin = MinLayerCell(ICell);
      const int KMax = MaxLayerCell(ICell);
      parallelForInner(
          Team, Range(KMin, KMax),
          INNER_LAMBDA(int K) { FCTTracerCurFill(L, ICell, K, TracerArray); });
   }

   KOKKOS_FUNCTION void FCTTracerCurFill(const I4 L, const I4 ICell, const I4 K,
                                         const Array3DReal &TracerArray) const {
      TracerCur(ICell, K) = TracerArray(L, ICell, K);
      TracerMin(ICell, K) = TracerArray(L, ICell, K);
      TracerMax(ICell, K) = TracerArray(L, ICell, K);
   }
   KOKKOS_FUNCTION void FCTTracerMinMax(const TeamMember &Team,
                                        const I4 ICell) const {
      const I4 KMin   = MinLayerCell(ICell);
      const I4 KMax   = MaxLayerCell(ICell);
      const I4 NEdges = NEdgesOnCell(ICell);
      for (I4 I = 0; I < NEdges; ++I) {
         const I4 ICell2 = CellsOnCell(ICell, I);
         const I4 KMin1  = Kokkos::max(KMin, MinLayerCell(ICell2));
         const I4 KMax1  = Kokkos::min(KMax, MaxLayerCell(ICell2));
         parallelForInner(
             Team, Range(KMin1, KMax1), INNER_LAMBDA(int K) {
                TracerMax(ICell, K) =
                    Kokkos::max(TracerMax(ICell, K), TracerCur(ICell2, K));
                TracerMin(ICell, K) =
                    Kokkos::min(TracerMin(ICell, K), TracerCur(ICell2, K));
             });
      }
   }

   KOKKOS_FUNCTION void
   FCTHighAndLowOrderFlux(const TeamMember &Team, const I4 IEdge,
                          const Array2DReal &FluxPseudoThickEdge,
                          const Array2DReal &NormVelEdge) const {

      const I4 ICell1 = CellsOnEdge(IEdge, 0);
      const I4 ICell2 = CellsOnEdge(IEdge, 1);

      const I4 KMin = MinLayerEdgeBot(IEdge);
      const I4 KMax = MaxLayerEdgeTop(IEdge);
      parallelForInner(
          Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
             HighOrderFlx(IEdge, K) = 0;
             LowOrderFlx(IEdge, K)  = 0;
          });
      // Compute 3rd or 4th fluxes where requested.
      for (int I = 0; I < NAdvCellsForEdge(IEdge); ++I) {
         const I4 ICell   = AdvCellsForEdge(IEdge, I);
         const Real Coef1 = AdvCoefs(I, IEdge);
         const Real Coef3 = AdvCoefs3rd(I, IEdge) * Coef3rdOrder;
         const I4 KMin    = MinLayerCell(ICell);
         const I4 KMax    = MaxLayerCell(ICell);
         parallelForInner(
             Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
                const Real NormalThicknessFlux =
                    FluxPseudoThickEdge(IEdge, K) * NormVelEdge(IEdge, K);
                const Real TracerWgt =
                    NormalThicknessFlux *
                    (Coef1 +
                     Coef3 * std::copysign(1.0_Real, NormalThicknessFlux));
                HighOrderFlx(IEdge, K) += TracerWgt * TracerCur(ICell, K) *
                                          AdvMaskHighOrder(IEdge, K);
             });
      }
      // Compute 2nd order fluxes where needed.
      // Also compute low order upwind horizontal flux (monotonic)
      // Remove low order flux from the high order flux
      // Store left over high order flux in highOrderFlx array
      parallelForInner(
          Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
             const Real NormalThicknessFlux =
                 FluxPseudoThickEdge(IEdge, K) * NormVelEdge(IEdge, K);
             const Real TracerWeight = (1.0_Real - AdvMaskHighOrder(IEdge, K)) *
                                       (DvEdge(IEdge) * 0.5_Real) *
                                       NormalThicknessFlux;
             LowOrderFlx(IEdge, K) =
                 DvEdge(IEdge) * (Kokkos::max(0.0_Real, NormalThicknessFlux) *
                                      TracerCur(ICell1, K) +
                                  Kokkos::min(0.0_Real, NormalThicknessFlux) *
                                      TracerCur(ICell2, K));
             HighOrderFlx(IEdge, K) +=
                 TracerWeight * (TracerCur(ICell1, K) + TracerCur(ICell2, K));
             HighOrderFlx(IEdge, K) -= LowOrderFlx(IEdge, K);
          });
   }

   KOKKOS_FUNCTION void FCTFluxInOut(const TeamMember &Team, const I4 ICell,
                                     const Real Dt,
                                     const Array2DReal &LayerThickness) const {
      const I4 KMin          = MinLayerCell(ICell);
      const I4 KMax          = MaxLayerCell(ICell);
      const Real InvAreaCell = 1._Real / AreaCell(ICell);
      parallelForInner(
          Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
             WorkTend(ICell, K) = 0;
             FlxIn(ICell, K)    = 0;
             FlxOut(ICell, K)   = 0;
          });
      // Finish computing the low order horizontal fluxes
      // Upwind fluxes are accumulated in workTend
      for (I4 I = 0; I < NEdgesOnCell(ICell); ++I) {
         const I4 IEdge          = EdgesOnCell(ICell, I);
         const Real SignedFactor = EdgeSignOnCell(ICell, I) * InvAreaCell;
         const I4 KMin           = MinLayerEdgeBot(IEdge);
         const I4 KMax           = MaxLayerEdgeTop(IEdge);
         parallelForInner(
             Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
                // Here workTend is the advection tendency due to the
                // upwind (low order) fluxes.
                WorkTend(ICell, K) += SignedFactor * LowOrderFlx(IEdge, K);
                // Accumulate remaining high order fluxes
                FlxOut(ICell, K) += Kokkos::min(
                    0.0_Real, SignedFactor * HighOrderFlx(IEdge, K));
                FlxIn(ICell, K) += Kokkos::max(
                    0.0_Real, SignedFactor * HighOrderFlx(IEdge, K));
             });
      }
      // Build the factors for the FCT
      // Computed using the bounds that were computed previously,
      // and the bounds on the newly updated value
      // Factors are placed in the flxIn and flxOut arrays
      parallelForInner(
          Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
             // Here workTend is the upwind tendency
             const Real TracerUpwindNew =
                 (TracerCur(ICell, K) * LayerThickness(ICell, K) +
                  Dt * WorkTend(ICell, K)) *
                 HProvInv(ICell, K);
             const Real TracerMinNew =
                 TracerUpwindNew + Dt * FlxOut(ICell, K) * HProvInv(ICell, K);
             const Real TracerMaxNew =
                 TracerUpwindNew + Dt * FlxIn(ICell, K) * HProvInv(ICell, K);
             const Real ScaleFactorIn =
                 (TracerMax(ICell, K) - TracerUpwindNew) /
                 (TracerMaxNew - TracerUpwindNew + Eps);
             FlxIn(ICell, K) =
                 Kokkos::min(1.0_Real, Kokkos::max(0.0_Real, ScaleFactorIn));
             const Real ScaleFactorOut =
                 (TracerUpwindNew - TracerMin(ICell, K)) /
                 (TracerUpwindNew - TracerMinNew + Eps);
             FlxOut(ICell, K) =
                 Kokkos::min(1.0_Real, Kokkos::max(0.0_Real, ScaleFactorOut));
          });
   }

   KOKKOS_FUNCTION void FCTRescaleHighOrderFlux(const TeamMember &Team,
                                                const I4 IEdge) const {
      const int KMin = MinLayerEdgeBot(IEdge);
      const int KMax = MaxLayerEdgeTop(IEdge);
      parallelForInner(
          Team, Range(KMin, KMax),
          INNER_LAMBDA(int K) { FCTRescaleHighOrderFlux(IEdge, K); });
   }
   KOKKOS_FUNCTION void FCTRescaleHighOrderFlux(const I4 IEdge,
                                                const I4 K) const {
      const I4 ICell1 = CellsOnEdge(IEdge, 0);
      const I4 ICell2 = CellsOnEdge(IEdge, 1);
      HighOrderFlx(IEdge, K) =
          Kokkos::max(0.0_Real, HighOrderFlx(IEdge, K)) *
              Kokkos::min(FlxOut(ICell1, K), FlxIn(ICell2, K)) +
          Kokkos::min(0.0_Real, HighOrderFlx(IEdge, K)) *
              Kokkos::min(FlxIn(ICell1, K), FlxOut(ICell2, K));
   }

   KOKKOS_FUNCTION void
   FCTAccumulateHighOrderFlux(const TeamMember &Team, const I4 ICell,
                              const Real Dt, const Array2DReal &TracerArray,
                              const Array2DReal &LayerThickness) const {

      // Accumulate the scaled high order vertical tendencies
      // and the upwind tendencies
      const Real InvAreaCell1 = 1.0_Real / AreaCell(ICell);
      // Accumulate the scaled high order horizontal tendencies
      for (I4 I = 0; I < NEdgesOnCell(ICell); ++I) {
         const I4 IEdge          = EdgesOnCell(ICell, I);
         const Real SignedFactor = InvAreaCell1 * EdgeSignOnCell(ICell, I);
         const I4 KMin           = MinLayerEdgeBot(IEdge);
         const I4 KMax           = MaxLayerEdgeTop(IEdge);
         parallelForInner(
             Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
                // WorkTend on RHS is upwind tendency
                // WorkTend on LHS is total horiz advect tendency
                WorkTend(ICell, K) += SignedFactor * HighOrderFlx(IEdge, K);
             });
      }
      const int KMin = MinLayerCell(ICell);
      const int KMax = MaxLayerCell(ICell);
      parallelForInner(
          Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
             // workTend  on RHS is total horiz advection tendency
             // TracerCur on LHS is provisional tracer after
             //                     horizontal fluxes only.
             TracerCur(ICell, K) =
                 (TracerCur(ICell, K) * LayerThickness(ICell, K) +
                  Dt * WorkTend(ICell, K)) *
                 HProvInv(ICell, K);
             TracerArray(ICell, K) += WorkTend(ICell, K);
          });
   }

   KOKKOS_FUNCTION void
   FCTComputeBudgetAdvectionEdgeFlux(const TeamMember &Team, const I4 L,
                                     const I4 IEdge) const {
      const int KMin = MinLayerEdgeBot(IEdge);
      const int KMax = MaxLayerEdgeTop(IEdge);
      parallelForInner(
          Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
             FCTComputeBudgetAdvectionEdgeFlux(L, IEdge, K);
          });
   }
   KOKKOS_FUNCTION void FCTComputeBudgetAdvectionEdgeFlux(const I4 L,
                                                          const I4 IEdge,
                                                          const I4 K) const {
      // Save u*h*T flux on edge for analysis. This variable will be
      // divided by h at the end of the time step.
      ActiveTracerHorizontalAdvectionEdgeFlux(L, IEdge, K) =
          (LowOrderFlx(IEdge, K) + HighOrderFlx(IEdge, K)) / DvEdge(IEdge);
   }

   KOKKOS_FUNCTION void
   FCTComputeBudgetAdvectionTendency(const TeamMember &Team, const I4 L,
                                     const I4 ICell) const {
      const int KMin = MinLayerCell(ICell);
      const int KMax = MaxLayerCell(ICell);
      parallelForInner(
          Team, Range(KMin, KMax), INNER_LAMBDA(int K) {
             FCTComputeBudgetAdvectionTendency(L, ICell, K);
          });
   }
   KOKKOS_FUNCTION void FCTComputeBudgetAdvectionTendency(const I4 L,
                                                          const I4 ICell,
                                                          const I4 K) const {
      ActiveTracerHorizontalAdvectionTendency(L, ICell, K) = WorkTend(ICell, K);
   }

   KOKKOS_FUNCTION void FCTMonotonicityCheck(const TeamMember &Team,
                                             const I4 ICell) const {
      const int KMin = MinLayerCell(ICell);
      const int KMax = MaxLayerCell(ICell);
      parallelForInner(
          Team, Range(KMin, KMax),
          INNER_LAMBDA(int K) { FCTMonotonicityCheck(ICell, K); });
   }
   KOKKOS_FUNCTION void FCTMonotonicityCheck(const I4 ICell, const I4 K) const {
      // Check tracer values against local min,max to detect
      // non-monotone values and write warning if found
      if (TracerCur(ICell, K) < TracerMin(ICell, K) - Eps) {
         Kokkos::printf(
             "Horizontal minimum out of bounds on cell %d, level %d "
             "new tracer value %lg is smaller than previous tracer minimum "
             "%lg\n",
             ICell, K, TracerCur(ICell, K), TracerMin(ICell, K));
      }
      if (TracerCur(ICell, K) > TracerMax(ICell, K) + Eps) {
         Kokkos::printf(
             "Horizontal maximum out of bounds on cell %d, level %d "
             "new tracer value %lg is larger than previous tracer maximum "
             "%lg\n",
             ICell, K, TracerCur(ICell, K), TracerMax(ICell, K));
      }
   }

 protected:
   const Real Eps = 1.e-10_Real;
   const HorzMesh *HorzontalMesh;
   const VertCoord *VerticalCoord;
   const I4 NVertLayers;

   Array1DI4 NAdvCellsForEdge;
   Array2DI4 AdvCellsForEdge;
   Array2DI4 AdvMaskHighOrder;
   Array2DI4 CellsOnCell;
   Array1DI4 NEdgesOnCell;
   Array2DReal AdvCoefs;
   Array2DReal AdvCoefs3rd;
   Array3DReal HighOrderFlxHorz;
   Array2DReal TracerCur;

   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnEdge;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
   Array2DReal EdgeSignOnCell;
   Array1DReal DvEdge;
   Array1DReal AreaCell;

   Array2DReal TotalVerticalTransportPseudoVelocity;
   Array2DReal HProvInv;
   Array2DReal HNewInv;
   Array2DReal HProv;
   Array2DReal TracerMax;
   Array2DReal TracerMin;
   Array2DReal HighOrderFlx;
   Array2DReal LowOrderFlx;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array2DReal WorkTend;
   Array2DReal FlxIn;
   Array2DReal FlxOut;
   Array3DReal ActiveTracerHorizontalAdvectionEdgeFlux;
   Array3DReal ActiveTracerHorizontalAdvectionTendency;
};

// Tracer horizontal diffusion term
class TracerDiffOnCell {
 public:
   bool Enabled = false;

   Real EddyDiff2;

   TracerDiffOnCell(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void
   operator()(const TeamMember &Team, const Array3DReal &Tend, I4 L, I4 ICell,
              const Array3DReal &TracerCell,
              const Array2DReal &MeanPseudoThickEdge) const {

      const auto LTend = subviewUnmanaged(Tend, L, Kokkos::ALL, Kokkos::ALL);
      const auto LTracerCell =
          subviewUnmanaged(TracerCell, L, Kokkos::ALL, Kokkos::ALL);

      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      ScratchArray1DReal DiffTmp(teamScratch(Team), NVertLayers);
      parallelForInner(
          Team, NVertLayers, INNER_LAMBDA(int K) { DiffTmp(K) = 0; });

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge = EdgesOnCell(ICell, J);

         const I4 JCell0 = CellsOnEdge(JEdge, 0);
         const I4 JCell1 = CellsOnEdge(JEdge, 1);

         const Real RTemp =
             MeshScalingDel2(JEdge) * DvEdge(JEdge) / DcEdge(JEdge);

         const int MinLyrEdgeBot = MinLayerEdgeBot(JEdge);
         const int MaxLyrEdgeTop = MaxLayerEdgeTop(JEdge);

         parallelForInner(
             Team, Range{MinLyrEdgeBot, MaxLyrEdgeTop}, INNER_LAMBDA(int K) {
                const Real TracerGrad =
                    (LTracerCell(JCell1, K) - LTracerCell(JCell0, K));

                DiffTmp(K) -= EdgeMask(JEdge, K) * EdgeSignOnCell(ICell, J) *
                              RTemp * MeanPseudoThickEdge(JEdge, K) *
                              TracerGrad;
             });
      }
      const int MinLyrCell = MinLayerCell(ICell);
      const int MaxLyrCell = MaxLayerCell(ICell);

      parallelForInner(
          Team, Range{MinLyrCell, MaxLyrCell}, INNER_LAMBDA(int K) {
             LTend(ICell, K) += EddyDiff2 * DiffTmp(K) * InvAreaCell;
          });
   }

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnEdge;
   Array2DReal EdgeSignOnCell;
   Array1DReal DvEdge;
   Array1DReal DcEdge;
   Array1DReal AreaCell;
   Array1DReal MeshScalingDel2;
   I4 NVertLayers;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};

/// KPP vertical non-local tracer flux tendency
class KPPNonLocalTracerFluxOnCell {
 public:
   bool Enabled = false;

   KPPNonLocalTracerFluxOnCell(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void operator()(const TeamMember &Team,
                                   const Array3DReal &Tend, I4 L, I4 ICell,
                                   const Array2DReal &SurfaceTracerFlux,
                                   const Array2DReal &VertNonLocalFlux) const {
      const int KMin = MinLayerCell(ICell);
      const int KMax = MaxLayerCell(ICell);

      parallelForInner(
          Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
             Tend(L, ICell, K) +=
                 SurfaceTracerFlux(L, ICell) *
                 (VertNonLocalFlux(ICell, K) - VertNonLocalFlux(ICell, K + 1));
          });
   }

 private:
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

// Tracer biharmonic horizontal mixing term
class TracerHyperDiffOnCell {
 public:
   bool Enabled = false;

   Real EddyDiff4;

   TracerHyperDiffOnCell(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void operator()(const TeamMember &Team,
                                   const Array3DReal &Tend, I4 L, I4 ICell,
                                   const Array3DReal &TrDel2Cell) const {

      const auto LTend = subviewUnmanaged(Tend, L, Kokkos::ALL, Kokkos::ALL);
      const auto LTrDel2Cell =
          subviewUnmanaged(TrDel2Cell, L, Kokkos::ALL, Kokkos::ALL);

      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      ScratchArray1DReal HypTmp(teamScratch(Team), NVertLayers);
      parallelForInner(
          Team, NVertLayers, INNER_LAMBDA(int K) { HypTmp(K) = 0; });

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge = EdgesOnCell(ICell, J);

         const I4 JCell0 = CellsOnEdge(JEdge, 0);
         const I4 JCell1 = CellsOnEdge(JEdge, 1);

         const Real RTemp =
             MeshScalingDel4(JEdge) * DvEdge(JEdge) / DcEdge(JEdge);

         const int MinLyrEdgeBot = MinLayerEdgeBot(JEdge);
         const int MaxLyrEdgeTop = MaxLayerEdgeTop(JEdge);

         parallelForInner(
             Team, Range{MinLyrEdgeBot, MaxLyrEdgeTop}, INNER_LAMBDA(int K) {
                const Real Del2TrGrad =
                    (LTrDel2Cell(JCell1, K) - LTrDel2Cell(JCell0, K));

                HypTmp(K) -= EdgeMask(JEdge, K) * EdgeSignOnCell(ICell, J) *
                             RTemp * Del2TrGrad;
             });
      }
      const int MinLyrCell = MinLayerCell(ICell);
      const int MaxLyrCell = MaxLayerCell(ICell);

      parallelForInner(
          Team, Range{MinLyrCell, MaxLyrCell}, INNER_LAMBDA(int K) {
             LTend(ICell, K) -= EddyDiff4 * HypTmp(K) * InvAreaCell;
          });
   }

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnEdge;
   Array2DReal EdgeSignOnCell;
   Array1DReal DvEdge;
   Array1DReal DcEdge;
   Array1DReal AreaCell;
   Array1DReal MeshScalingDel4;
   I4 NVertLayers;
   Array2DReal EdgeMask;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
};
/// Penetrating shortwave radiation forcing for conservative temperature.
class PenetratingShortwaveOnCell {
 public:
   bool Enabled        = false;
   Real NearIrFraction = 0.58_Real;
   Real NearIrCoeff    = 2.86_Real;
   Real RedFraction    = 0.21_Real;
   Real BlueFraction   = 0.21_Real;

   PenetratingShortwaveOnCell(const HorzMesh *Mesh, const VertCoord *VCoord,
                              I4 TempTracerIndex);

   KOKKOS_FUNCTION void
   operator()(const TeamMember &Team, const Array3DReal &Tend, I4 ICell,
              const Array2DReal &GeomZInterface,
              const Array1DReal &ShortWaveHeatFlux,
              const Array1DReal &ExtinctionCoeffRed,
              const Array1DReal &ExtinctionCoeffBlue) const {

      const I4 KTop = MinLayerCell(ICell);
      const I4 KBot = MaxLayerCell(ICell);
      if (KTop > KBot || TempIndex < 0) {
         return;
      }

      const Real SurfaceFlux = ShortWaveHeatFlux(ICell);
      const Real Kr          = ExtinctionCoeffRed(ICell);
      const Real Kb          = ExtinctionCoeffBlue(ICell);
      const Real ZSurface    = GeomZInterface(ICell, KTop);

      // First compute the shortwave flux crossing each layer interface and
      // store in scratch. The top interface receives the full surface flux and
      // the bottom interface is set to zero.
      ScratchArray1DReal FluxAtInterface(teamScratch(Team), NVertLayers + 1);
      parallelForInner(
          Team, Range{KTop, KBot + 1}, INNER_LAMBDA(int K) {
             if (K == KTop) {
                FluxAtInterface(K) = SurfaceFlux;
             } else if (K == KBot + 1) {
                FluxAtInterface(K) = 0.0_Real;
             } else {
                const Real Depth =
                    Kokkos::abs(GeomZInterface(ICell, K) - ZSurface);
                FluxAtInterface(K) =
                    SurfaceFlux *
                    (NearIrFraction * Kokkos::exp(-NearIrCoeff * Depth) +
                     RedFraction * Kokkos::exp(-Kr * Depth) +
                     BlueFraction * Kokkos::exp(-Kb * Depth));
             }
          });

      teamBarrier(Team);

      // Deposit the flux divergence across each layer as the heating tendency.
      parallelForInner(
          Team, Range{KTop, KBot}, INNER_LAMBDA(int K) {
             Tend(TempIndex, ICell, K) +=
                 (FluxAtInterface(K) - FluxAtInterface(K + 1)) * HFluxFac;
          });
   }

 private:
   I4 TempIndex;
   I4 NVertLayers;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

/// Surface tracer restoring term
class SurfaceTracerRestoringOnCell {
 public:
   bool Enabled;
   Real PistonVelocity  = 1.585e-5; ///< piston velocity
   I4 NTracersToRestore = 0;        ///< number of tracers to restore
   Array1DI4 TracerIdsToRestore;    ///< tracer IDs to restore
   /// Need to add under sea ice restoring option when that is available

   /// constructor declaration
   SurfaceTracerRestoringOnCell(const HorzMesh *Mesh);

   /// The functor takes the cell index and the array for the tracer surface
   /// restoring values, outputs tendency array
   KOKKOS_FUNCTION void
   operator()(const Array3DReal &Tend, I4 L, I4 ICell, I4 KMin,
              const Array2DReal &TracersMonthlySurfClimoCell,
              const Array3DReal &TracerCell) const {

      Tend(L, ICell, KMin) +=
          PistonVelocity *
          (TracersMonthlySurfClimoCell(L, ICell) - TracerCell(L, ICell, KMin));
   }
};

/// Frazil tendency hook-up term
class FrazilOnCell {
 public:
   bool Enabled = false;

   FrazilOnCell(const HorzMesh *Mesh, const VertCoord *VCoord);

   void operator()(const Array2DReal &PseudoThicknessTend,
                   const Array3DReal &TracerTend,
                   const Array3DReal &TracerArray,
                   const Array2DReal &PressureMid,
                   const Array2DReal &PseudoThickness,
                   const TimeInterval &TimeStep,
                   Real FinalUpdateWeight = 1.0_Real) const;

 private:
   I4 NCellsAll;
   I4 TempTracerIndex;
   I4 SaltTracerIndex;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

} // namespace OMEGA
#endif
