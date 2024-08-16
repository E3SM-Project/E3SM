#ifndef OMEGA_OCEANSTATE_H
#define OMEGA_OCEANSTATE_H
//===-- ocn/OceanState.h - ocean state --------------------*- C++ -*-===//
//
/// \file
/// \brief Contains the state variables for an OMEGA sub-domain
///
/// The OceanState class contains the non-tracer prognostic variable data for a
/// sub-domain of the global horizontal mesh.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "Decomp.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "MachEnv.h"

#include <string>

namespace OMEGA {

/// A class for the ocean prognostic variable information

/// The OceanState class provides a container for the layer thickness,
/// and normal velocity variables. It contains methods which handle
/// IO and time level updates.
class OceanState {

 private:
   void initParallelIO(Decomp *MeshDecomp);

   void finalizeParallelIO();

   void read();

   void defineFields();

   I4 CellDecompR8;
   I4 EdgeDecompR8;
   Halo *MeshHalo;

   static OceanState *DefaultOceanState;

   static std::map<std::string, std::unique_ptr<OceanState>> AllOceanStates;

   /// Construct a new local state for a given decomposition
   OceanState(const std::string &Name, ///< [in] Name for mesh
              HorzMesh *Mesh,          ///< [in] Horizontal mesh
              Decomp *MeshDecomp,      ///< [in] Decomp for Mesh
              Halo *MeshHalo_,         ///< [in] Halo for Mesh
              const int NVertLevels_,  ///< [in] Number of vertical levels
              const int NTimeLevels_   ///< [in] Number of time levels
   );

   // Forbid copy and move construction
   OceanState(const OceanState &) = delete;
   OceanState(OceanState &&)      = delete;

 public:
   // Variables
   // Since these are used frequently, we make them public to reduce the
   // number of retrievals required.

   std::string StateFileName;
   std::string Name;
   int StateFileID;

   // Sizes and global IDs
   // Note that all sizes are actual counts (1-based) so that loop extents
   // should always use the 0:NCellsXX-1 form.

   I4 NCellsOwned; ///< Number of cells owned by this task
   I4 NCellsAll;   ///< Total number of local cells (owned + all halo)
   I4 NCellsSize;  ///< Array size (incl padding, bndy cell) for cell arrays

   I4 NEdgesOwned; ///< Number of edges owned by this task
   I4 NEdgesAll;   ///< Total number (owned+halo) of local edges
   I4 NEdgesSize;  ///< Array length (incl padding, bndy) for edge dim

   static const I4 MaxTimeLevels = 5; ///< Maximum number of time levels

   I4 NTimeLevels; ///< Number of time levels in state variable arrays
   I4 NVertLevels; ///< Number of vertical levels in state variable arrays

   I4 CurLevel; ///< Time dimension index for current level
   I4 NewLevel; ///< Time dimension index for new level

   // Prognostic variables

   Kokkos::Array<Array2DReal, MaxTimeLevels>
       LayerThickness; ///< Device LayerThickness array
   Kokkos::Array<HostArray2DReal, MaxTimeLevels>
       LayerThicknessH; ///< Host LayerThickness array

   Kokkos::Array<Array2DReal, MaxTimeLevels>
       NormalVelocity; ///< Device NormalVelocity array
   Kokkos::Array<HostArray2DReal, MaxTimeLevels>
       NormalVelocityH; ///< Host NormalVelocity array

   // Field names
   // These are appended with the State name for non-Default state instances
   std::string LayerThicknessFldName; ///< Field name for LayerThickness
   std::string NormalVelocityFldName; ///< Field name for NormalVelocity
   std::string StateGroupName;

   // Methods

   /// Initialize Omega local state
   static int init();

   /// Create a new state by calling the constructor and put it in the
   /// AllOceanStates map
   static OceanState *
   create(const std::string &Name, ///< [in] Name for mesh
          HorzMesh *Mesh,          ///< [in] Horizontal mesh
          Decomp *MeshDecomp,      ///< [in] Decomp for Mesh
          Halo *MeshHalo,          ///< [in] Halo for Mesh
          const int NVertLevels,   ///< [in] Number of vertical levels
          const int NTimeLevels    ///< [in] Number of time levels
   );

   /// Swap time levels to update state arrays
   void updateTimeLevels();

   /// Copy state variables from host to device
   void copyToDevice(int TimeLevel);

   /// Copy state variables from device to host
   void copyToHost(int TimeLevel);

   /// Destructor - deallocates all memory and deletes an OceanState
   ~OceanState();

   /// Deallocates arrays
   static void clear();

   /// Remove state by name
   static void erase(std::string InName ///< [in] name of state to remove
   );

   static OceanState *getDefault();

   static OceanState *get(std::string name);

}; // end class OceanState

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_OCEANSTATE_H
