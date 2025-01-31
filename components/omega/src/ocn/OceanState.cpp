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

   if (NTimeLevels < 2) {
      LOG_ERROR("OceanState: the number of time level is lower than 2");
      return -2;
   }

   // Create the default state and set pointer to it
   OceanState::DefaultOceanState =
       create("Default", DefHorzMesh, DefHalo, NVertLevels, NTimeLevels);

   // State values are filled by a later read of the initial state or
   // restart file

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

   MeshHalo = MeshHalo_;

   Name = Name_;

   CurTimeIndex = 0;

   // Allocate state host arrays
   LayerThicknessH.resize(NTimeLevels);
   NormalVelocityH.resize(NTimeLevels);

   for (int I = 0; I < NTimeLevels; I++) {
      LayerThicknessH[I] = HostArray2DReal("LayerThickness" + std::to_string(I),
                                           NCellsSize, NVertLevels);
      NormalVelocityH[I] = HostArray2DReal("NormalVelocity" + std::to_string(I),
                                           NEdgesSize, NVertLevels);
   }

   // Allocate state device arrays
   LayerThickness.resize(NTimeLevels);
   NormalVelocity.resize(NTimeLevels);

   // Create device arrays and copy host data
   for (int I = 0; I < NTimeLevels; I++) {
      LayerThickness[I] = Array2DReal("LayerThickness" + std::to_string(I),
                                      NCellsSize, NVertLevels);
      NormalVelocity[I] = Array2DReal("NormalVelocity" + std::to_string(I),
                                      NEdgesSize, NVertLevels);
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

   // Add restart group if needed
   if (!FieldGroup::exists("Restart"))
      auto RestartGroup = FieldGroup::create("Restart");

   Err = StateGroup->addField(NormalVelocityFldName);
   if (Err != 0)
      LOG_ERROR("Error adding {} to field group {}", NormalVelocityFldName,
                StateGroupName);
   Err = StateGroup->addField(LayerThicknessFldName);
   if (Err != 0)
      LOG_ERROR("Error adding {} to field group {}", LayerThicknessFldName,
                StateGroupName);

   Err = FieldGroup::addFieldToGroup(NormalVelocityFldName, "Restart");
   if (Err != 0)
      LOG_ERROR("Error adding {} to Restart field group",
                NormalVelocityFldName);
   Err = FieldGroup::addFieldToGroup(LayerThicknessFldName, "Restart");
   if (Err != 0)
      LOG_ERROR("Error adding {} to Restart field group",
                LayerThicknessFldName);

   // Associate Field with data
   I4 TimeIndex;
   Err = getTimeIndex(TimeIndex, 0);

   Err =
       NormalVelocityField->attachData<Array2DReal>(NormalVelocity[TimeIndex]);
   if (Err != 0)
      LOG_ERROR("Error attaching data array to field {}",
                NormalVelocityFldName);
   Err =
       LayerThicknessField->attachData<Array2DReal>(LayerThickness[TimeIndex]);
   if (Err != 0)
      LOG_ERROR("Error attaching data array to field {}",
                LayerThicknessFldName);

} // end defineIOFields

//------------------------------------------------------------------------------
// Get layer thickness device array
I4 OceanState::getLayerThickness(Array2DReal &LayerThick,
                                 const I4 TimeLevel) const {
   I4 Err = 0;
   I4 TimeIndex;

   Err        = getTimeIndex(TimeIndex, TimeLevel);
   LayerThick = LayerThickness[TimeIndex];

   return Err;
}

//------------------------------------------------------------------------------
// Get layer thickness host array
I4 OceanState::getLayerThicknessH(HostArray2DReal &LayerThickH,
                                  const I4 TimeLevel) const {
   I4 Err = 0;
   I4 TimeIndex;

   Err         = getTimeIndex(TimeIndex, TimeLevel);
   LayerThickH = LayerThicknessH[TimeIndex];

   return Err;
}

//------------------------------------------------------------------------------
// Get normal velocity device array
I4 OceanState::getNormalVelocity(Array2DReal &NormVel,
                                 const I4 TimeLevel) const {

   I4 Err = 0;
   I4 TimeIndex;

   Err     = getTimeIndex(TimeIndex, TimeLevel);
   NormVel = NormalVelocity[TimeIndex];

   return Err;
}

//------------------------------------------------------------------------------
// Get normal velocity host array
I4 OceanState::getNormalVelocityH(HostArray2DReal &NormVelH,
                                  const I4 TimeLevel) const {
   I4 Err = 0;
   I4 TimeIndex;

   Err      = getTimeIndex(TimeIndex, TimeLevel);
   NormVelH = NormalVelocityH[TimeIndex];

   return Err;
}

//------------------------------------------------------------------------------
// Perform copy to device for state variables
// TimeLevel == [1:new, 0:current, -1:previous, -2:two times ago, ...]
I4 OceanState::copyToDevice(const I4 TimeLevel) {

   I4 Err = 0;
   I4 TimeIndex;

   Err = getTimeIndex(TimeIndex, TimeLevel);

   deepCopy(LayerThickness[TimeIndex], LayerThicknessH[TimeIndex]);
   deepCopy(NormalVelocity[TimeIndex], NormalVelocityH[TimeIndex]);

   return 0;
} // end copyToDevice

//------------------------------------------------------------------------------
// Perform copy to host for state variables
// TimeLevel == [1: new, 0:current, -1:previous, -2:two times ago, ...]
I4 OceanState::copyToHost(const I4 TimeLevel) {

   I4 Err = 0;
   I4 TimeIndex;

   Err = getTimeIndex(TimeIndex, TimeLevel);

   deepCopy(LayerThicknessH[TimeIndex], LayerThickness[TimeIndex]);
   deepCopy(NormalVelocityH[TimeIndex], NormalVelocity[TimeIndex]);

   return 0;
} // end copyToHost

//------------------------------------------------------------------------------
// Perform state halo exchange
// TimeLevel == [1:new, 0:current, -1:previous, -2:two times ago, ...]
I4 OceanState::exchangeHalo(const I4 TimeLevel) {

   I4 Err = 0;
   I4 TimeIndex;
   Err = getTimeIndex(TimeIndex, TimeLevel);

   MeshHalo->exchangeFullArrayHalo(LayerThickness[TimeIndex], OnCell);
   MeshHalo->exchangeFullArrayHalo(NormalVelocity[TimeIndex], OnEdge);

   return 0;

} // end exchangeHalo

//------------------------------------------------------------------------------
// Perform time level update
I4 OceanState::updateTimeLevels() {

   if (NTimeLevels == 1) {
      LOG_ERROR("OceanState: can't update time levels for NTimeLevels == 1");
      return -1;
   }

   // Exchange halo
   exchangeHalo(1);

   // Update current time index for layer thickness and normal velocity
   CurTimeIndex = (CurTimeIndex + 1) % NTimeLevels;

   I4 Err = 0;

   // Update IOField data associations
   Err = Field::attachFieldData<Array2DReal>(NormalVelocityFldName,
                                             NormalVelocity[CurTimeIndex]);
   Err = Field::attachFieldData<Array2DReal>(LayerThicknessFldName,
                                             LayerThickness[CurTimeIndex]);

   return 0;

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

//------------------------------------------------------------------------------
// Get time index from time level
// TimeLevel == [1:new, 0:current, -1:previous, -2:two times ago, ...]
I4 OceanState::getTimeIndex(I4 &TimeIndex, const I4 TimeLevel) const {

   // Check if time level is valid
   if (NTimeLevels > 1 && (TimeLevel > 1 || (TimeLevel + NTimeLevels) <= 1)) {
      LOG_ERROR("OceanState: Time level {} is out of range for NTimeLevels {}",
                TimeLevel, NTimeLevels);
      return -1;
   }

   TimeIndex = (TimeLevel + CurTimeIndex + NTimeLevels) % NTimeLevels;

   return 0;

} // end get time index

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
