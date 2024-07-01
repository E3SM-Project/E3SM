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
   Decomp *DefDecomp = Decomp::getDefault();
   HorzMesh *DefHorzMesh = HorzMesh::getDefault();
   Halo *DefHalo = Halo::getDefault();

   // Create the default state 
   OceanState DefOceanState("Default", DefHorzMesh, DefDecomp, DefHalo, 60, 2);

   // Retrieve this mesh and set pointer to DefaultOceanState
   OceanState::DefaultOceanState = OceanState::get("Default");
   return Err;
}

//------------------------------------------------------------------------------
// Construct a new local mesh given a decomposition

OceanState::OceanState(const std::string &Name,   //< [in] Name for new state
                       HorzMesh *Mesh,            //< [in] HorzMesh for state
                       Decomp *MeshDecomp,        //< [in] Decomp for Mesh
                       Halo *MeshHalo_,           //< [in] Halo for Mesh
                       const int VerticalLevels_, //< [in] number of vertical levels
                       const int TimeLevels_      //< [in] number of time levels
) {

   // Retrieve mesh cell/edge/vertex totals from Decomp
   NCellsOwned = Mesh->NCellsOwned;
   NCellsAll   = Mesh->NCellsAll;
   NCellsSize  = Mesh->NCellsSize;

   NEdgesOwned    = Mesh->NEdgesOwned;
   NEdgesAll      = Mesh->NEdgesAll;
   NEdgesSize     = Mesh->NEdgesSize;

   VerticalLevels = VerticalLevels_;
   TimeLevels = TimeLevels_;
 
   MeshHalo = MeshHalo_;

   StateFileName = Mesh->MeshFileName;

   // Allocate state arrays
   LayerThicknessH = HostArray3DR8("LayerThickness", TimeLevels, NCellsSize, VerticalLevels);
   NormalVelocityH = HostArray3DR8("NormalVelocity", TimeLevels, NEdgesSize, VerticalLevels);

   // Create device arrays
   LayerThickness = createDeviceMirrorCopy(LayerThicknessH);
   NormalVelocity = createDeviceMirrorCopy(NormalVelocityH);

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
   copyToDevice();

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
   std::vector<I4> CellDims{1, MeshDecomp->NCellsGlobal, VerticalLevels};
   std::vector<I4> CellID(NCellsAll * VerticalLevels, -1);
   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      for(int Level = 0; Level < VerticalLevels; ++Level) {
         I4 GlobalID = (MeshDecomp->CellIDH(Cell) - 1) * VerticalLevels + Level;
         CellID[Cell * VerticalLevels + Level] = GlobalID;
      }
   }

   Err = IO::createDecomp(CellDecompR8, IO::IOTypeR8, NDims, CellDims,
                          NCellsAll * VerticalLevels, CellID, Rearr);
   if (Err != 0)
      LOG_CRITICAL("OceanState: error creating cell IO decomposition");

   // Create the IO decomp for arrays with (NEdges) dimensions
   std::vector<I4> EdgeDims{1, MeshDecomp->NEdgesGlobal, VerticalLevels};
   std::vector<I4> EdgeID(NEdgesAll * VerticalLevels, -1);
   for (int Edge = 0; Edge < NEdgesAll; ++Edge) {
      for (int Level = 0; Level < VerticalLevels; ++Level) {
         I4 GlobalID = (MeshDecomp->EdgeIDH(Edge) - 1) * VerticalLevels + Level;
         EdgeID[Edge * VerticalLevels + Level] = GlobalID;
      }
   }

   Err = IO::createDecomp(EdgeDecompR8, IO::IOTypeR8, NDims, EdgeDims,
                          NEdgesAll * VerticalLevels, EdgeID, Rearr);
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
   auto LayerThicknessSubH = Kokkos::subview(LayerThicknessH, 0, Kokkos::ALL, Kokkos::ALL);
   Err    = IO::readArray(LayerThicknessSubH.data(), NCellsAll, "layerThickness", StateFileID,
                          CellDecompR8, LayerThicknessID);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading xCell");

   // Read NormalVelocity
   int NormalVelocityID;
   auto NormalVelocitySubH = Kokkos::subview(NormalVelocityH, 0, Kokkos::ALL, Kokkos::ALL);
   Err    = IO::readArray(NormalVelocitySubH.data(), NEdgesAll, "normalVelocity", StateFileID,
                          EdgeDecompR8, NormalVelocityID);

} // end read


//------------------------------------------------------------------------------
// Perform copy to device for state variables
void OceanState::copyToDevice() {

   deepCopy(LayerThickness, LayerThicknessH);
   deepCopy(NormalVelocity, NormalVelocityH);

} // end copyToDevice

//------------------------------------------------------------------------------
// Perform copy to host for state variables
void OceanState::copyToHost() {

   deepCopy(LayerThicknessH, LayerThickness);
   deepCopy(NormalVelocityH, NormalVelocity);

} // end copyToHost

//------------------------------------------------------------------------------
// Perform time level swap
void OceanState::swapTimeLevels(int FromLevel, int ToLevel) {
  
   copyToHost(); 
   MeshHalo->exchangeFullArrayHalo(LayerThicknessH, OMEGA::OnCell);
   copyToDevice();

   auto LayerThicknessSubTo = Kokkos::subview(LayerThickness, ToLevel, Kokkos::ALL, Kokkos::ALL);
   auto LayerThicknessTemp = Kokkos::create_mirror_view(LayerThicknessSubTo);

   auto LayerThicknessSubFrom = Kokkos::subview(LayerThickness, FromLevel, Kokkos::ALL, Kokkos::ALL);

   Kokkos::deep_copy(LayerThicknessSubTo, LayerThicknessSubFrom);
   Kokkos::deep_copy(LayerThicknessSubFrom, LayerThicknessTemp);

   copyToHost();
   MeshHalo->exchangeFullArrayHalo(NormalVelocityH, OMEGA::OnEdge); 
   copyToDevice();

   auto NormalVelocitySubTo = Kokkos::subview(NormalVelocity, ToLevel, Kokkos::ALL, Kokkos::ALL);
   auto NormalVelocityTemp = Kokkos::create_mirror_view(NormalVelocitySubTo);

   auto NormalVelocitySubFrom = Kokkos::subview(NormalVelocity, FromLevel, Kokkos::ALL, Kokkos::ALL);

   Kokkos::deep_copy(NormalVelocitySubTo, NormalVelocitySubFrom);
   Kokkos::deep_copy(NormalVelocitySubFrom, NormalVelocityTemp);

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
