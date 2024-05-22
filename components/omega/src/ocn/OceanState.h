#ifndef OMEGA_OCEANSTATE_H
#define OMEGA_OCEANSTATE_H
//===-- ocn/OceanState.h - ocean state --------------------*- C++ -*-===//
//
/// \file
/// \brief Contains the state variables for an OMEGA sub-domain
///
/// The OceanState class contains the prognostic variable data for a
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

   I4 CellDecompR8;
   I4 EdgeDecompR8;
   Halo *MeshHalo;

   static OceanState *DefaultOceanState;

   static std::map<std::string, OceanState> AllOceanStates;

 public:
   // Variables
   // Since these are used frequently, we make them public to reduce the
   // number of retrievals required.

   std::string StateFileName;
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

   I4 NTimeLevels;     ///< Number of time levels in state variable arrays
   I4 NVerticalLevels; ///< Number of vertical levels in state variable arrays

   // Prognostic variables

   std::vector<Array2DR8> LayerThickness;      ///< Device LayerThickness array
   std::vector<HostArray2DR8> LayerThicknessH; ///< Host LayerThickness array

   std::vector<Array2DR8> NormalVelocity;      ///< Device LayerThickness array
   std::vector<HostArray2DR8> NormalVelocityH; ///< Host LayterThickness array

   // Methods

   /// Initialize Omega local mesh
   static int init();

   /// Construct a new local mesh for a given decomposition
   OceanState(const std::string &Name,    ///< [in] Name for mesh
              HorzMesh *Mesh,             ///< [in] Horizontal mesh
              Decomp *MeshDecomp,         ///< [in] Decomp for Mesh
              Halo *MeshHalo_,            ///< [in] Halo for Mesh
              const int NVerticalLevels_, ///< [in] Number of vertical levels
              const int NTimeLevels_      ///< [in] Number of time levels
   );

   /// Swap time levels to update state arrays
   void swapTimeLevels(int FromLevel, int ToLevel);

   void copyToDevice(int TimeLevel);

   void copyToHost(int TimeLevel);

   /// Destructor - deallocates all memory and deletes a HorzMesh
   ~OceanState();

   /// Deallocates arrays
   static void clear();

   /// Remove mesh by name
   static void erase(std::string InName ///< [in] name of mesh to remove
   );

   static OceanState *getDefault();

   static OceanState *get(std::string name);

}; // end class OceanState

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_OCEANSTATE_H
