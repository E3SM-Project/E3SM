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

#include <functional>
#include <memory>

namespace OMEGA {

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
   KEGradOnEdge KEGrad;
   SSHGradOnEdge SSHGrad;
   VelocityDiffusionOnEdge VelocityDiffusion;
   VelocityHyperDiffOnEdge VelocityHyperDiff;
   WindForcingOnEdge WindForcing;
   BottomDragOnEdge BottomDrag;
   TracerHorzAdvOnCell TracerHorzAdv;
   TracerDiffOnCell TracerDiffusion;
   TracerHyperDiffOnCell TracerHyperDiff;
   SurfaceTracerRestoringOnCell SurfaceTracerRestoring;

   std::string Name;

   // Methods to compute tendency groups
   void computeThicknessTendencies(const OceanState *State,
                                   const AuxiliaryState *AuxState,
                                   int ThickTimeLevel, int VelTimeLevel,
                                   TimeInstant Time);
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
                                TimeInstant Time);
   void computeAllTendencies(const OceanState *State,
                             const AuxiliaryState *AuxState,
                             const Array3DReal &TracerArray, int ThickTimeLevel,
                             int VelTimeLevel, int TracerTimeLevel,
                             TimeInstant Time, TimeInterval ProjDt);
   void computeThicknessTendenciesOnly(const OceanState *State,
                                       const AuxiliaryState *AuxState,
                                       int ThickTimeLevel, int VelTimeLevel,
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
                                    TimeInstant Time);

   // Create a non-default group of tendencies
   template <class... ArgTypes>
   static Tendencies *create(const std::string &Name, ArgTypes &&...Args) {
      // Check to see if tendencies of the same name already exist and
      // if so, exit with an error
      if (AllTendencies.find(Name) != AllTendencies.end()) {
         LOG_ERROR(
             "Attempted to create Tendencies with name {} but Tendencies of "
             "that name already exists",
             Name);
         return nullptr;
      }

      // create new tendencies on the heap and put it in a map of
      // unique_ptrs, which will manage its lifetime
      auto *NewTendencies =
          new Tendencies(Name, std::forward<ArgTypes>(Args)...);
      AllTendencies.emplace(Name, NewTendencies);

      return get(Name);
   }

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
              int NTracersIn,          ///< [in] Number of tracers
              TimeInterval TimeStep,   ///< [in] Time step
              Config *Options,         ///< [in] Configuration options
              CustomTendencyType InCustomThicknessTend,
              CustomTendencyType InCustomVelocityTend);

   Tendencies(const std::string &Name, ///< [in] Name for tendencies
              const HorzMesh *Mesh,    ///< [in] Horizontal mesh
              VertCoord *VCoord,       ///< [in] Vertical coordinate
              VertAdv *VAdv,           ///< [in] Vertical advection
              PressureGrad *PGrad,     ///< [in] Pressure gradient
              Eos *EqState,            ///< [in] Equation of state
              int NTracersIn,          ///< [in] Number of tracers
              TimeInterval TimeStep,   ///< [in] Time step
              Config *Options          ///< [in] Configuration options
   );

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
   I4 NTracers;           ///< Number of tracers
   TimeInterval TimeStep; ///< Time step

   // Pointer to default tendencies
   static Tendencies *DefaultTendencies;

   // Map of all tendency objects
   static std::map<std::string, std::unique_ptr<Tendencies>> AllTendencies;

}; // end class Tendencies

} // namespace OMEGA
#endif
