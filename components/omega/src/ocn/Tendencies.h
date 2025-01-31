#ifndef OMEGA_TENDENCIES_H
#define OMEGA_TENDENCIES_H
//===-- ocn/Tendencies.h - Tendencies --------------------*- C++ -*-===//
//
/// \file
/// \brief Manages the tendencies for state variables and tracers
///
/// The Tendencies class contains the tendency data for state variables and
/// tracers and provides methods for computing different tendency groups.
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "Config.h"
#include "HorzMesh.h"
#include "OceanState.h"
#include "TendencyTerms.h"
#include "TimeMgr.h"

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
   Array2DReal LayerThicknessTend;
   Array2DReal NormalVelocityTend;
   Array3DReal TracerTend;

   // Instances of tendency terms
   ThicknessFluxDivOnCell ThicknessFluxDiv;
   PotentialVortHAdvOnEdge PotientialVortHAdv;
   KEGradOnEdge KEGrad;
   SSHGradOnEdge SSHGrad;
   VelocityDiffusionOnEdge VelocityDiffusion;
   VelocityHyperDiffOnEdge VelocityHyperDiff;
   TracerHorzAdvOnCell TracerHorzAdv;
   TracerDiffOnCell TracerDiffusion;
   TracerHyperDiffOnCell TracerHyperDiff;

   // Methods to compute tendency groups
   void computeThicknessTendencies(const OceanState *State,
                                   const AuxiliaryState *AuxState,
                                   int ThickTimeLevel, int VelTimeLevel,
                                   TimeInstant Time);
   void computeVelocityTendencies(const OceanState *State,
                                  const AuxiliaryState *AuxState,
                                  int ThickTimeLevel, int VelTimeLevel,
                                  TimeInstant Time);
   void computeTracerTendencies(const OceanState *State,
                                const AuxiliaryState *AuxState,
                                const Array3DReal &TracerArray,
                                int ThickTimeLevel, int VelTimeLevel,
                                TimeInstant Time);
   void computeAllTendencies(const OceanState *State,
                             const AuxiliaryState *AuxState,
                             const Array3DReal &TracerArray, int ThickTimeLevel,
                             int VelTimeLevel, TimeInstant Time);
   void computeThicknessTendenciesOnly(const OceanState *State,
                                       const AuxiliaryState *AuxState,
                                       int ThickTimeLevel, int VelTimeLevel,
                                       TimeInstant Time);
   void computeVelocityTendenciesOnly(const OceanState *State,
                                      const AuxiliaryState *AuxState,
                                      int ThickTimeLevel, int VelTimeLevel,
                                      TimeInstant Time);
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
   static int init();

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
   int readTendConfig(Config *TendConfig);

 private:
   // Construct a new tendency object
   Tendencies(const std::string &Name, ///< [in] Name for tendencies
              const HorzMesh *Mesh,    ///< [in] Horizontal mesh
              int NVertLevels,         ///< [in] Number of vertical levels
              int NTracersIn,          ///< [in] Number of tracers
              Config *Options,         ///< [in] Configuration options
              CustomTendencyType InCustomThicknessTend,
              CustomTendencyType InCustomVelocityTend);

   Tendencies(const std::string &Name, ///< [in] Name for tendencies
              const HorzMesh *Mesh,    ///< [in] Horizontal mesh
              int NVertLevels,         ///< [in] Number of vertical levels
              int NTracersIn,          ///< [in] Number of tracers
              Config *Options          ///< [in] Configuration options
   );

   // forbid copy and move construction
   Tendencies(const Tendencies &) = delete;
   Tendencies(Tendencies &&)      = delete;

   // Mesh sizes
   I4 NCellsAll; ///< Number of cells including full halo
   I4 NEdgesAll; ///< Number of edges including full halo
   I4 NTracers;  ///< Number of tracers
   I4 NChunks;   ///< Number of vertical level chunks

   // Pointer to default tendencies
   static Tendencies *DefaultTendencies;

   // Map of all tendency objects
   static std::map<std::string, std::unique_ptr<Tendencies>> AllTendencies;

   CustomTendencyType CustomThicknessTend;
   CustomTendencyType CustomVelocityTend;

}; // end class Tendencies

} // namespace OMEGA
#endif
