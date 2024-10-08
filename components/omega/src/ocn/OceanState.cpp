//===-- ocn/OceanState.cpp - ocean state methods -------------*- C++ -*-===//
//
// The OceanState class initializes the prognostic variables in OMEGA.
// It contains a method to update the time levels for each variable.
// It is meant to provide a container for passing (non-tracer) prognostic
// variables throughout the OMEGA tendency computation routines.
//
//===----------------------------------------------------------------------===//

#include "OceanState.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Field.h"
#include "Halo.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "TimeStepper.h"

namespace OMEGA {

// create the static class members
OceanState *OceanState::DefaultOceanState = nullptr;
std::map<std::string, std::unique_ptr<OceanState>> OceanState::AllOceanStates;

//------------------------------------------------------------------------------
// Initialize the state. Assumes that Decomp has already been initialized.

int OceanState::init() {

   int Err = 0; // default successful return code

   // Retrieve the default decomposition and mesh
   Decomp *DefDecomp     = Decomp::getDefault();
   HorzMesh *DefHorzMesh = HorzMesh::getDefault();
   Halo *DefHalo         = Halo::getDefault();

   int NVertLevels = DefHorzMesh->NVertLevels;

   auto *DefTimeStepper = TimeStepper::getDefault();
   if (!DefTimeStepper) {
      LOG_ERROR("TimeStepper needs to be initialized before OceanState");
   }
   int NTimeLevels = DefTimeStepper->getNTimeLevels();

   // Create the default state and set pointer to it
   OceanState::DefaultOceanState =
       create("Default", DefHorzMesh, DefHalo, NVertLevels, NTimeLevels);

   DefaultOceanState->loadStateFromFile(DefHorzMesh->MeshFileName, DefDecomp);

   return Err;
}

//------------------------------------------------------------------------------
// Construct a new local state

OceanState::OceanState(
    const std::string &Name_, //< [in] Name for new state
    HorzMesh *Mesh,           //< [in] HorzMesh for state
    Halo *MeshHalo_,          //< [in] Halo for Mesh
    const int NVertLevels_,   //< [in] number of vertical levels
    const int NTimeLevels_    //< [in] number of time levels
) {

   // Retrieve mesh cell/edge/vertex totals from Decomp
   NCellsOwned = Mesh->NCellsOwned;
   NCellsAll   = Mesh->NCellsAll;
   NCellsSize  = Mesh->NCellsSize;

   NEdgesOwned = Mesh->NEdgesOwned;
   NEdgesAll   = Mesh->NEdgesAll;
   NEdgesSize  = Mesh->NEdgesSize;

   NVertLevels = NVertLevels_;
   NTimeLevels = NTimeLevels_;
   CurLevel    = std::max(0, NTimeLevels - 2);
   NewLevel    = NTimeLevels - 1;

   MeshHalo = MeshHalo_;

   Name = Name_;

   // Allocate state host arrays
   for (int I = 0; I < NTimeLevels; I++) {
      LayerThicknessH[I] = HostArray2DReal("LayerThickness" + std::to_string(I),
                                           NCellsSize, NVertLevels);
      NormalVelocityH[I] = HostArray2DReal("NormalVelocity" + std::to_string(I),
                                           NEdgesSize, NVertLevels);
   }

   // Create device arrays and copy host data
   for (int I = 0; I < NTimeLevels; I++) {
      LayerThickness[I] = createDeviceMirrorCopy(LayerThicknessH[I]);
      NormalVelocity[I] = createDeviceMirrorCopy(NormalVelocityH[I]);
   }

   // Register fields and metadata for IO
   defineFields();

} // end state constructor

/// Create a new state by calling the constructor and put it in the
/// AllOceanStates map
OceanState *
OceanState::create(const std::string &Name, //< [in] Name for new state
                   HorzMesh *Mesh,          //< [in] HorzMesh for state
                   Halo *MeshHalo,          //< [in] Halo for Mesh
                   const int NVertLevels,   //< [in] number of vertical levels
                   const int NTimeLevels    //< [in] number of time levels
) {

   // Check to see if a state of the same name already exists and
   // if so, exit with an error
   if (AllOceanStates.find(Name) != AllOceanStates.end()) {
      LOG_ERROR(
          "Attempted to create an OceanState with name {} but an OceanState of "
          "that name already exists",
          Name);
      return nullptr;
   }

   // create a new state on the heap and put it in a map of
   // unique_ptrs, which will manage its lifetime
   auto *NewOceanState =
       new OceanState(Name, Mesh, MeshHalo, NVertLevels, NTimeLevels);
   AllOceanStates.emplace(Name, NewOceanState);

   return NewOceanState;
} // end state create

//------------------------------------------------------------------------------
// Destroys a local mesh and deallocates all arrays
OceanState::~OceanState() {

   // Kokkos arrays removed when no longer in scope

   int Err;
   Err = FieldGroup::destroy(StateGroupName);
   if (Err != 0)
      LOG_ERROR("Error removing FieldGroup {}", StateGroupName);
   Err = Field::destroy(LayerThicknessFldName);
   if (Err != 0)
      LOG_ERROR("Error removing Field {}", LayerThicknessFldName);
   Err = Field::destroy(NormalVelocityFldName);
   if (Err != 0)
      LOG_ERROR("Error removing Field {}", NormalVelocityFldName);

} // end destructor

//------------------------------------------------------------------------------
// Removes a state from list by name
void OceanState::erase(std::string InName // [in] name of state to remove
) {

   AllOceanStates.erase(InName); // remove the state from the list and in
                                 // the process, calls the destructor

} // end state erase
//------------------------------------------------------------------------------
// Removes all states to clean up before exit
void OceanState::clear() {

   AllOceanStates.clear(); // removes all states from the list and in
                           // the process, calls the destructors for each

} // end clear

//------------------------------------------------------------------------------
// Load state from file
void OceanState::loadStateFromFile(const std::string &StateFileName,
                                   Decomp *MeshDecomp) {

   int StateFileID;
   I4 CellDecompR8;
   I4 EdgeDecompR8;

   // Open the state file for reading (assume IO has already been initialized)
   I4 Err;
   Err = IO::openFile(StateFileID, StateFileName, IO::ModeRead);
   if (Err != 0)
      LOG_CRITICAL("OceanState: error opening state file");

   // Create the parallel IO decompositions required to read in state variables
   initParallelIO(CellDecompR8, EdgeDecompR8, MeshDecomp);

   // Read layerThickness and normalVelocity
   read(StateFileID, CellDecompR8, EdgeDecompR8);

   // Destroy the parallel IO decompositions
   finalizeParallelIO(CellDecompR8, EdgeDecompR8);

   // Sync with device
   copyToDevice(CurLevel);
} // end loadStateFromFile

//------------------------------------------------------------------------------
// Initialize the parallel IO decompositions for the mesh variables
void OceanState::initParallelIO(I4 &CellDecompR8, I4 &EdgeDecompR8,
                                Decomp *MeshDecomp) {

   I4 Err;
   I4 NDims             = 3;
   IO::Rearranger Rearr = IO::RearrBox;

   // Create the IO decomp for arrays with (NCells) dimensions
   std::vector<I4> CellDims{1, MeshDecomp->NCellsGlobal, NVertLevels};
   std::vector<I4> CellID(NCellsAll * NVertLevels, -1);
   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      for (int Level = 0; Level < NVertLevels; ++Level) {
         I4 GlobalID = (MeshDecomp->CellIDH(Cell) - 1) * NVertLevels + Level;
         CellID[Cell * NVertLevels + Level] = GlobalID;
      }
   }

   Err = IO::createDecomp(CellDecompR8, IO::IOTypeR8, NDims, CellDims,
                          NCellsAll * NVertLevels, CellID, Rearr);
   if (Err != 0)
      LOG_CRITICAL("OceanState: error creating cell IO decomposition");

   // Create the IO decomp for arrays with (NEdges) dimensions
   std::vector<I4> EdgeDims{1, MeshDecomp->NEdgesGlobal, NVertLevels};
   std::vector<I4> EdgeID(NEdgesAll * NVertLevels, -1);
   for (int Edge = 0; Edge < NEdgesAll; ++Edge) {
      for (int Level = 0; Level < NVertLevels; ++Level) {
         I4 GlobalID = (MeshDecomp->EdgeIDH(Edge) - 1) * NVertLevels + Level;
         EdgeID[Edge * NVertLevels + Level] = GlobalID;
      }
   }

   Err = IO::createDecomp(EdgeDecompR8, IO::IOTypeR8, NDims, EdgeDims,
                          NEdgesAll * NVertLevels, EdgeID, Rearr);
   if (Err != 0)
      LOG_CRITICAL("OceanState: error creating edge IO decomposition");

} // end initParallelIO

//------------------------------------------------------------------------------
// Define IO fields and metadata
void OceanState::defineFields() {

   int Err = 0;

   LayerThicknessFldName = "LayerThickness";
   NormalVelocityFldName = "NormalVelocity";
   if (Name != "Default") {
      LayerThicknessFldName.append(Name);
      NormalVelocityFldName.append(Name);
   }

   // Create fields for state variables
   int NDims = 2;
   std::vector<std::string> DimNames(NDims);
   DimNames[0] = "NEdges";
   DimNames[1] = "NVertLevels";
   auto NormalVelocityField =
       Field::create(NormalVelocityFldName,               // field name
                     "Velocity component normal to edge", // long Name
                     "m/s",                               // units
                     "sea_water_velocity",                // CF standard Name
                     -9.99E+10,                           // min valid value
                     9.99E+10,                            // max valid value
                     -9.99E+30, // scalar for undefined entries
                     NDims,     // number of dimensions
                     DimNames   // dimension names
       );

   DimNames[0] = "NCells";
   auto LayerThicknessField =
       Field::create(LayerThicknessFldName,               // Field name
                     "Thickness of layer on cell center", /// long Name
                     "m",                                 // units
                     "cell_thickness",                    // CF standard Name
                     0.0,                                 // min valid value
                     9.99E+30,                            // max valid value
                     -9.99E+30, // scalar used for undefined entries
                     NDims,     // number of dimensions
                     DimNames   // dimension names
       );

   // Create a field group for state fields
   StateGroupName = "State";
   if (Name != "Default") {
      StateGroupName.append(Name);
   }
   auto StateGroup = FieldGroup::create(StateGroupName);

   Err = StateGroup->addField(NormalVelocityFldName);
   if (Err != 0)
      LOG_ERROR("Error adding {} to field group {}", NormalVelocityFldName,
                StateGroupName);
   Err = StateGroup->addField(LayerThicknessFldName);
   if (Err != 0)
      LOG_ERROR("Error adding {} to field group {}", LayerThicknessFldName,
                StateGroupName);

   // Associate Field with data
   Err = NormalVelocityField->attachData<Array2DReal>(NormalVelocity[CurLevel]);
   if (Err != 0)
      LOG_ERROR("Error attaching data array to field {}",
                NormalVelocityFldName);
   Err = LayerThicknessField->attachData<Array2DReal>(LayerThickness[CurLevel]);
   if (Err != 0)
      LOG_ERROR("Error attaching data array to field {}",
                LayerThicknessFldName);

} // end defineIOFields

//------------------------------------------------------------------------------
// Destroy parallel decompositions
void OceanState::finalizeParallelIO(I4 CellDecompR8, I4 EdgeDecompR8) {

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
void OceanState::read(int StateFileID, I4 CellDecompR8, I4 EdgeDecompR8) {

   I4 Err;

   // Read LayerThickness into a temporary double-precision array
   int LayerThicknessID;
   HostArray2DR8 TmpLayerThicknessR8("TmpLayerThicknessR8", NCellsSize,
                                     NVertLevels);
   Err = IO::readArray(TmpLayerThicknessR8.data(), NCellsAll, "layerThickness",
                       StateFileID, CellDecompR8, LayerThicknessID);
   if (Err != 0)
      LOG_CRITICAL("OceanState: error reading layerThickness");

   // Copy the thickness data into the final state array of user-specified
   // precision
   deepCopy(LayerThicknessH[CurLevel], TmpLayerThicknessR8);

   // Read NormalVelocity  into a temporary double-precision array
   int NormalVelocityID;
   HostArray2DR8 TmpNormalVelocityR8("TmpNormalVelocityR8", NEdgesSize,
                                     NVertLevels);
   Err = IO::readArray(TmpNormalVelocityR8.data(), NEdgesAll, "normalVelocity",
                       StateFileID, EdgeDecompR8, NormalVelocityID);
   if (Err != 0)
      LOG_CRITICAL("OceanState: error reading normalVelocity");

   // Copy the velocity data into the final state array of user-specified
   // precision
   deepCopy(NormalVelocityH[CurLevel], TmpNormalVelocityR8);

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
// Perform state halo exchange
void OceanState::exchangeHalo(int TimeLevel) {
   copyToHost(TimeLevel);
   MeshHalo->exchangeFullArrayHalo(LayerThicknessH[TimeLevel], OnCell);
   MeshHalo->exchangeFullArrayHalo(NormalVelocityH[TimeLevel], OnEdge);
   copyToDevice(TimeLevel);
} // end exchangeHalo

//------------------------------------------------------------------------------
// Perform time level update
void OceanState::updateTimeLevels() {

   int NewLevel = NTimeLevels - 1;

   // Exchange halo
   exchangeHalo(NewLevel);

   // Update time levels for layer thickness and normal velocity
   for (int Level = 0; Level < NTimeLevels - 1; Level++) {
      std::swap(LayerThickness[Level + 1], LayerThickness[Level]);
      std::swap(LayerThicknessH[Level + 1], LayerThicknessH[Level]);

      std::swap(NormalVelocity[Level + 1], NormalVelocity[Level]);
      std::swap(NormalVelocityH[Level + 1], NormalVelocityH[Level]);
   }

   // Update IOField data associations
   int Err = 0;

   Err = Field::attachFieldData<Array2DReal>(NormalVelocityFldName,
                                             NormalVelocity[CurLevel]);
   Err = Field::attachFieldData<Array2DReal>(LayerThicknessFldName,
                                             LayerThickness[CurLevel]);

} // end updateTimeLevels

//------------------------------------------------------------------------------
// Get default state
OceanState *OceanState::getDefault() { return OceanState::DefaultOceanState; }

//------------------------------------------------------------------------------
// Get state by name
OceanState *OceanState::get(const std::string Name ///< [in] Name of state
) {

   // look for an instance of this name
   auto it = AllOceanStates.find(Name);

   // if found, return the state pointer
   if (it != AllOceanStates.end()) {
      return it->second.get();

      // otherwise print error and return null pointer
   } else {
      LOG_ERROR("OceanState::get: Attempt to retrieve non-existent state:");
      LOG_ERROR("{} has not been defined or has been removed", Name);
      return nullptr;
   }
} // end get state

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
