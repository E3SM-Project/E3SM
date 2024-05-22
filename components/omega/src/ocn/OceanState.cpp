//===-- ocn/OceanState.cpp - ocean state methods -------------*- C++ -*-===//
//
// The OceanState class initializes the prognostic variables in OMEGA.
// It contains a method to update the time levels for each variable.
// It is meant to provide a container for passing prognostic variables
// throughout the OMEGA tendency computation routines.
//
//===----------------------------------------------------------------------===//

#include "OceanState.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Halo.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"

namespace OMEGA {

// create the static class members
OceanState *OceanState::DefaultOceanState = nullptr;
std::map<std::string, OceanState> OceanState::AllOceanStates;

//------------------------------------------------------------------------------
// Initialize the state. Assumes that Decomp has already been initialized.

int OceanState::init() {

   int Err = 0; // default successful return code

   // Retrieve the default decomposition and mesh
   Decomp *DefDecomp     = Decomp::getDefault();
   HorzMesh *DefHorzMesh = HorzMesh::getDefault();
   Halo *DefHalo         = Halo::getDefault();

   // Create the default state
   OceanState DefOceanState("Default", DefHorzMesh, DefDecomp, DefHalo, 60, 2);

   // Retrieve this mesh and set pointer to DefaultOceanState
   OceanState::DefaultOceanState = OceanState::get("Default");
   return Err;
}

//------------------------------------------------------------------------------
// Construct a new local mesh given a decomposition

OceanState::OceanState(
    const std::string &Name,    //< [in] Name for new state
    HorzMesh *Mesh,             //< [in] HorzMesh for state
    Decomp *MeshDecomp,         //< [in] Decomp for Mesh
    Halo *MeshHalo_,            //< [in] Halo for Mesh
    const int NVerticalLevels_, //< [in] number of vertical levels
    const int NTimeLevels_      //< [in] number of time levels
    )
    : LayerThicknessH(NTimeLevels_), NormalVelocityH(NTimeLevels_),
      LayerThickness(NTimeLevels_), NormalVelocity(NTimeLevels_) {

   // Retrieve mesh cell/edge/vertex totals from Decomp
   NCellsOwned = Mesh->NCellsOwned;
   NCellsAll   = Mesh->NCellsAll;
   NCellsSize  = Mesh->NCellsSize;

   NEdgesOwned = Mesh->NEdgesOwned;
   NEdgesAll   = Mesh->NEdgesAll;
   NEdgesSize  = Mesh->NEdgesSize;

   NVerticalLevels = NVerticalLevels_;
   NTimeLevels     = NTimeLevels_;

   MeshHalo = MeshHalo_;

   StateFileName = Mesh->MeshFileName;

   // Allocate state arrays
   for (int I = 0; I < NTimeLevels; I++) {
      LayerThicknessH[I] = HostArray2DR8("LayerThickness" + std::to_string(I),
                                         NCellsSize, NVerticalLevels);
      NormalVelocityH[I] = HostArray2DR8("NormalVelocity" + std::to_string(I),
                                         NEdgesSize, NVerticalLevels);
   }

   // Create device arrays
   for (int I = 0; I < NTimeLevels; I++) {
      LayerThickness[I] = Array2DR8("LayerThickness" + std::to_string(I),
                                    NCellsSize, NVerticalLevels);
      NormalVelocity[I] = Array2DR8("NormalVelocity" + std::to_string(I),
                                    NEdgesSize, NVerticalLevels);
   }

   // Open the mesh file for reading (assume IO has already been initialized)
   I4 Err;
   Err = OMEGA::IO::openFile(StateFileID, StateFileName, IO::ModeRead);
   if (Err != 0)
      LOG_CRITICAL("OceanState: error opening state file");

   // Create the parallel IO decompositions required to read in mesh variables
   initParallelIO(MeshDecomp);

   // Read x/y/z and lon/lat coordinates for cells, edges, and vertices
   read();

   // Destroy the parallel IO decompositions
   finalizeParallelIO();

   // Copy host data to device
   copyToDevice(0);

   // Associate this instance with a name
   AllOceanStates.emplace(Name, *this);

} // end horizontal mesh constructor

//------------------------------------------------------------------------------
// Destroys a local mesh and deallocates all arrays
OceanState::~OceanState() {

   // No operations needed, Kokkos arrays removed when no longer in scope

} // end deconstructor

//------------------------------------------------------------------------------
// Removes a mesh from list by name
void OceanState::erase(std::string InName // [in] name of mesh to remove
) {

   AllOceanStates.erase(InName); // remove the mesh from the list and in
                                 // the process, calls the destructor

} // end mesh erase
//------------------------------------------------------------------------------
// Removes all meshes to clean up before exit
void OceanState::clear() {

   AllOceanStates.clear(); // removes all meshes from the list and in
                           // the porcess, calls the destructors for each

} // end clear

//------------------------------------------------------------------------------
// Initialize the parallel IO decompositions for the mesh variables
void OceanState::initParallelIO(Decomp *MeshDecomp) {

   I4 Err;
   I4 NDims             = 3;
   IO::Rearranger Rearr = IO::RearrBox;

   // Create the IO decomp for arrays with (NCells) dimensions
   std::vector<I4> CellDims{1, MeshDecomp->NCellsGlobal, NVerticalLevels};
   std::vector<I4> CellID(NCellsAll * NVerticalLevels, -1);
   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      for (int Level = 0; Level < NVerticalLevels; ++Level) {
         I4 GlobalID =
             (MeshDecomp->CellIDH(Cell) - 1) * NVerticalLevels + Level;
         CellID[Cell * NVerticalLevels + Level] = GlobalID;
      }
   }

   Err = IO::createDecomp(CellDecompR8, IO::IOTypeR8, NDims, CellDims,
                          NCellsAll * NVerticalLevels, CellID, Rearr);
   if (Err != 0)
      LOG_CRITICAL("OceanState: error creating cell IO decomposition");

   // Create the IO decomp for arrays with (NEdges) dimensions
   std::vector<I4> EdgeDims{1, MeshDecomp->NEdgesGlobal, NVerticalLevels};
   std::vector<I4> EdgeID(NEdgesAll * NVerticalLevels, -1);
   for (int Edge = 0; Edge < NEdgesAll; ++Edge) {
      for (int Level = 0; Level < NVerticalLevels; ++Level) {
         I4 GlobalID =
             (MeshDecomp->EdgeIDH(Edge) - 1) * NVerticalLevels + Level;
         EdgeID[Edge * NVerticalLevels + Level] = GlobalID;
      }
   }

   Err = IO::createDecomp(EdgeDecompR8, IO::IOTypeR8, NDims, EdgeDims,
                          NEdgesAll * NVerticalLevels, EdgeID, Rearr);
   if (Err != 0)
      LOG_CRITICAL("OceanStateh: error creating edge IO decomposition");

} // end initParallelIO

//------------------------------------------------------------------------------
// Destroy parallel decompositions
void OceanState::finalizeParallelIO() {

   int Err = 0; // default return code

   // Destroy the IO decomp for arrays with (NCells) dimensions
   Err = IO::destroyDecomp(CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("OceanState: error destroying cell IO decomposition");

   // Destroy the IO decomp for arrays with (NEdges) dimensions
   Err = IO::destroyDecomp(EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("OceanState: error destroying edge IO decomposition");

} // end finalizeParallelIO

//------------------------------------------------------------------------------
// Read Ocean State
void OceanState::read() {

   I4 Err;

   // Read LayerThickness
   int LayerThicknessID;
   Err = IO::readArray(LayerThicknessH[0].data(), NCellsAll, "layerThickness",
                       StateFileID, CellDecompR8, LayerThicknessID);
   if (Err != 0)
      LOG_CRITICAL("OceanState: error reading layerThickness");

   // Read NormalVelocity
   int NormalVelocityID;
   Err = IO::readArray(NormalVelocityH[0].data(), NEdgesAll, "normalVelocity",
                       StateFileID, EdgeDecompR8, NormalVelocityID);
   if (Err != 0)
      LOG_CRITICAL("OceanState: error reading normalVelocity");

} // end read

//------------------------------------------------------------------------------
// Perform copy to device for state variables
void OceanState::copyToDevice(int TimeLevel) {

   deepCopy(LayerThickness[TimeLevel], LayerThicknessH[TimeLevel]);
   deepCopy(NormalVelocity[TimeLevel], NormalVelocityH[TimeLevel]);

} // end copyToDevice

//------------------------------------------------------------------------------
// Perform copy to host for state variables
void OceanState::copyToHost(int TimeLevel) {

   deepCopy(LayerThicknessH[TimeLevel], LayerThickness[TimeLevel]);
   deepCopy(NormalVelocityH[TimeLevel], NormalVelocity[TimeLevel]);

} // end copyToHost

//------------------------------------------------------------------------------
// Perform time level swap
void OceanState::swapTimeLevels(int FromLevel, int ToLevel) {

   copyToHost(FromLevel);
   MeshHalo->exchangeFullArrayHalo(LayerThicknessH[FromLevel], OMEGA::OnCell);
   copyToDevice(FromLevel);

   Array2DR8 Temp            = LayerThickness[ToLevel];
   LayerThickness[ToLevel]   = LayerThickness[FromLevel];
   LayerThickness[FromLevel] = Temp;

   copyToHost(FromLevel);
   MeshHalo->exchangeFullArrayHalo(NormalVelocityH[FromLevel], OMEGA::OnEdge);
   copyToDevice(FromLevel);

   Temp                      = NormalVelocity[ToLevel];
   NormalVelocity[ToLevel]   = NormalVelocity[FromLevel];
   NormalVelocity[FromLevel] = Temp;

} // end swapTimeLevels

//------------------------------------------------------------------------------
// Get default state
OceanState *OceanState::getDefault() { return OceanState::DefaultOceanState; }

//------------------------------------------------------------------------------
// Get state by name
OceanState *OceanState::get(const std::string Name ///< [in] Name of state
) {

   // look for an instance of this name
   auto it = AllOceanStates.find(Name);

   // if found, return the mesh pointer
   if (it != AllOceanStates.end()) {
      return &(it->second);

      // otherwise print error and return null pointer
   } else {
      LOG_ERROR("OceanState::get: Attempt to retrieve non-existent state:");
      LOG_ERROR("{} has not been defined or has been removed", Name);
      return nullptr;
   }
} // end get state

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
