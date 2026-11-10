#ifndef OMEGA_TENDENCIES_H
#define OMEGA_TENDENCIES_H
//===-- ocn/Tendencies.h - Tendencies --------------------*- C++ -*-===//
//
/// \file
/// \brief Manages the tendencies for state variables and tracers
///
/// The Tendencies class contains the tendency data for state variables and
/// tracers and provides methods for computing different tendency groups.
/// Tendencies are configured in the input configuration file using:
/// \ConfigInput
/// # Sample tendencies input configuration (for Default config)
/// Tendencies:
///    ThicknessFluxTendencyEnable: true
///    PVTendencyEnable: true
///    KETendencyEnable: true
///    SSHTendencyEnable: true
///    VelDiffTendencyEnable: true
///    ViscDel2: 1.0e3
///    VelHyperDiffTendencyEnable: true
///    ViscDel4: 1.2e11
///    DivFactor: 1.0
///    BottomDragTendency:
///       Enable: true
///       Mode: Implicit
///       Type: Constant
///       BottomDragCoeff: 1.0e-3
///    TracerHorzAdvTendencyEnable: true
///    TracerDiffTendencyEnable: true
///    EddyDiff2: 10.0
///    TracerHyperDiffTendencyEnable: true
///    EddyDiff4: 0.0
///    UseCustomTendency: false
///    ManufacturedSolutionTendency: false
/// \EndConfigInput
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "Config.h"
#include "Eos.h"
#include "HorzMesh.h"
#include "OceanState.h"
#include "PGrad.h"
#include "TendencyTerms.h"
#include "TimeMgr.h"
#include "VertAdv.h"
#include "VertCoord.h"
#include "VertMix.h"

#include <functional>
#include <memory>

namespace OMEGA {

/// Describes how the Coriolis force enters the normal velocity tendency.
enum class CoriolisTendMode {
   PVFlux,  ///< planetary vorticity carried inside the PV flux term (default)
   Separate ///< relative vorticity only, Coriolis applied by the time stepper
};

/// A class that can be used to calculate the thickness,
/// velocity, and tracer tendencies within the timestepping algorithm.
class Tendencies {
 public:
   using CustomTendencyType =
       std::function<void(Array2DReal, const OceanState *,
                          const AuxiliaryState *, int, int, TimeInstant)>;
   // Arrays for accumulating tendencies
   Array2DReal PseudoThicknessTend;
   Array2DReal NormalVelocityTend;
   Array3DReal TracerTend;

   // Instances of tendency terms
   PseudoThicknessFluxDivOnCell PseudoThicknessFluxDiv;
   PotentialVortHAdvOnEdge PotentialVortHAdv;
   CoriolisAccelerationOnEdge CoriolisAcceleration;
   KEGradOnEdge KEGrad;
   SSHGradOnEdge SSHGrad;
   VelocityDiffusionOnEdge VelocityDiffusion;
   VelocityHyperDiffOnEdge VelocityHyperDiff;
   SfcStressForcingOnEdge SfcStressForcing;
   BottomDragOnEdge ExplicitBottomDrag;
   SfcThicknessForcingOnCell SfcThicknessForcing;
   SfcTracerForcingOnCell SfcTracerForcing;
   PenetratingShortwaveOnCell PenetratingShortwave;
   TracerHorzAdvOnCell TracerHorzAdv;
   TracerDiffOnCell TracerDiffusion;
   KPPNonLocalTracerFluxOnCell KPPNonLocalTracerFlux;
   TracerHyperDiffOnCell TracerHyperDiff;
   SurfaceTracerRestoringOnCell SurfaceTracerRestoring;
   FrazilOnCell FrazilTerm;

   /// Mode-split configuration of the velocity tendency
   ///  - Coriolis treatment in the vorticity flux term
   CoriolisTendMode CoriolisMode = CoriolisTendMode::PVFlux;
   //   - The split factor for the barotropic pressure anomaly gradient
   Real SplitFactor = 0._Real;

   // Diagnostics for temperature forcing pathways used in KPP comparison.
   // These are raw contributions added to TracerTend before tracer update.
   Array2DReal KPPNonLocalTracerTempTendDiag;
   Array1DReal KPPNonLocalTracerTempColumnSumDiag;

   // Enable diagnostics that isolate temperature non-local terms.
   bool KPPNonLocalTracerDiagnosticsEnable = true;

   // Enable vertical mixing tendencies
   bool VelVertMixTendencyEnable    = false;
   bool TracerVertMixTendencyEnable = false;

   std::string Name;

   /// Configure the velocity tendency for a mode-split time stepper
   void setModeSplit(CoriolisTendMode Mode, Real SplitFactorIn);

   // Methods to compute tendency groups
   void computePseudoThicknessTendencies(const OceanState *State,
                                         const AuxiliaryState *AuxState,
                                         const Array3DReal &TracerArray,
                                         int ThickTimeLevel, int VelTimeLevel,
                                         TimeInstant Time, TimeInterval ProjDt);
   void computeVelocityTendencies(const OceanState *State,
                                  const AuxiliaryState *AuxState,
                                  const Array3DReal &TracerArray,
                                  int ThickTimeLevel, int VelTimeLevel,
                                  int TracerTimeLevel, TimeInstant Time,
                                  TimeInterval ProjDt);
   void computeTracerTendencies(const OceanState *State,
                                const AuxiliaryState *AuxState,
                                const Array3DReal &TracerArray,
                                int ThickTimeLevel, int VelTimeLevel,
                                TimeInstant Time, TimeInterval ProjDt,
                                Real FinalUpdateWeight = 1.0_Real);
   void computeAllTendencies(const OceanState *State,
                             const AuxiliaryState *AuxState,
                             const Array3DReal &TracerArray, int ThickTimeLevel,
                             int VelTimeLevel, int TracerTimeLevel,
                             TimeInstant Time, TimeInterval ProjDt,
                             Real FinalUpdateWeight = 1.0_Real);

   void computeAllTendenciesOnly(const OceanState *State,
                                 const AuxiliaryState *AuxState,
                                 const Array3DReal &TracerArray,
                                 int ThickTimeLevel, int VelTimeLevel,
                                 int TracerTimeLevel, TimeInstant Time,
                                 Real FinalUpdateWeight = 1.0_Real);

   void computePseudoThicknessTendenciesOnly(const OceanState *State,
                                             const AuxiliaryState *AuxState,
                                             int ThickTimeLevel,
                                             int VelTimeLevel,
                                             TimeInstant Time);
   void computeVelocityTendenciesOnly(const OceanState *State,
                                      const AuxiliaryState *AuxState,
                                      const Array3DReal &TracerArray,
                                      int ThickTimeLevel, int VelTimeLevel,
                                      int TracerTimeLevel, TimeInstant Time);
   void computeTracerTendenciesOnly(const OceanState *State,
                                    const AuxiliaryState *AuxState,
                                    const Array3DReal &TracerArray,
                                    int ThickTimeLevel, int VelTimeLevel,
                                    TimeInstant Time,
                                    Real FinalUpdateWeight = 1.0_Real);
   void computeCoriolisAccelerationOnEdge(
       const Array2DReal &Tend,          ///< [inout] velocity tendency
       const Array2DReal &NormalVelEdge, ///< [in] normal velocity on edges
       const Array1DReal &FEdge          ///< [in] Coriolis parameter on edges
   ) const;
   void computeCoriolisAccelerationOnEdge(
       const Array2DReal &Tend,          ///< [out] base plus Coriolis tendency
       const Array2DReal &BaseTend,      ///< [in] tendency without Coriolis
       const Array2DReal &NormalVelEdge, ///< [in] normal velocity on edges
       const Array1DReal &FEdge          ///< [in] Coriolis parameter on edges
   ) const;
   void computeCoriolisAccelerationOnEdge(
       const Array1DReal &Tend, ///< [inout] barotropic velocity tendency
       const Array1DReal &NormalVelEdge, ///< [in] normal velocity on edges
       const Array1DReal &FEdge          ///< [in] Coriolis parameter on edges
   ) const;

   // Create a non-default group of tendencies
   static Tendencies *
   create(const std::string &Name, ///< [in] Name for tendencies
          const HorzMesh *Mesh,    ///< [in] Horizontal mesh
          VertCoord *VCoord,       ///< [in] Vertical coordinate
          VertAdv *VAdv,           ///< [in] Vertical advection
          PressureGrad *PGrad,     ///< [in] Pressure gradient
          Eos *EqState,            ///< [in] Equation of state
          VertMix *VMix,           ///< [in] Vertical mixing
          int NTracersIn,          ///< [in] Number of tracers
          TimeInterval TimeStep,   ///< [in] Time step
          Config *Options,         ///< [in] Configuration options
          CustomTendencyType CustomThicknessTend = CustomTendencyType{},
          CustomTendencyType CustomVelocityTend  = CustomTendencyType{});

   // Destructor
   ~Tendencies();

   // Initialize Omega tendencies
   static void init();

   // Deallocates arrays
   static void clear();

   // Remove tendencies object by name
   static void erase(const std::string &Name ///< [in]
   );

   // get default tendencies
   static Tendencies *getDefault();

   // get tendencies by name
   static Tendencies *get(const std::string &Name ///< [in]
   );

   // check if tendencies exists by name
   static bool exists(const std::string &Name ///< [in]
   );

   // read and set config options
   void readConfig(Config *OmegaConfig);

 private:
   // Construct a new tendency object
   Tendencies(const std::string &Name, ///< [in] Name for tendencies
              const HorzMesh *Mesh,    ///< [in] Horizontal mesh
              VertCoord *VCoord,       ///< [in] Vertical coordinate
              VertAdv *VAdv,           ///< [in] Vertical advection
              PressureGrad *PGrad,     ///< [in] Pressure gradient
              Eos *EqState,            ///< [in] Equation of state
              VertMix *VMix,           ///< [in] Vertical mixing
              int NTracersIn,          ///< [in] Number of tracers
              TimeInterval TimeStep,   ///< [in] Time step
              Config *Options,         ///< [in] Configuration options
              CustomTendencyType InCustomThicknessTend,
              CustomTendencyType InCustomVelocityTend);

   void defineFields();

   // forbid copy and move construction
   Tendencies(const Tendencies &) = delete;
   Tendencies(Tendencies &&)      = delete;

   const HorzMesh *Mesh; ///< Pointer to horizontal mesh
   VertCoord *VCoord;    ///< Pointer to vertical coordinate
   VertAdv *VAdv;        ///< Pointer to vertical advection
   CustomTendencyType CustomThicknessTend;
   CustomTendencyType CustomVelocityTend;
   Eos *EqState;          ///< Pointer to equation of state
   PressureGrad *PGrad;   ///< Pointer to pressure gradient
   VertMix *VMix;         ///< Pointer to vertical mixing
   I4 NTracers;           ///< Number of tracers
   TimeInterval TimeStep; ///< Time step

   /// Tracers after horizontal advection alone, the starting point for
   /// flux-corrected vertical advection
   Array3DReal ProvTracers;

   // Pointer to default tendencies
   static Tendencies *DefaultTendencies;

   // Map of all tendency objects
   static std::map<std::string, std::unique_ptr<Tendencies>> AllTendencies;

}; // end class Tendencies

} // namespace OMEGA
#endif
