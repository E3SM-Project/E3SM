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
#include "Halo.h"
#include "IO.h"
#include "IOField.h"
#include "Logging.h"
#include "MachEnv.h"
#include "MetaData.h"
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

   // These hard-wired variables need to be updated
   // with retrivals/config options
   int NTimeLevels = 2;
   int NVertLevels = 60;

   // Create the default state
   OceanState DefOceanState("Default", DefHorzMesh, DefDecomp, DefHalo,
                            NVertLevels, NTimeLevels);

   // Retrieve this mesh and set pointer to DefaultOceanState
   OceanState::DefaultOceanState = OceanState::get("Default");
   return Err;
}

//------------------------------------------------------------------------------
// Construct a new local state given a decomposition

OceanState::OceanState(
    const std::string &Name_, //< [in] Name for new state
    HorzMesh *Mesh,           //< [in] HorzMesh for state
    Decomp *MeshDecomp,       //< [in] Decomp for Mesh
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

   StateFileName = Mesh->MeshFileName;
   Name = Name_;

   // Allocate state host arrays
   for (int I = 0; I < NTimeLevels; I++) {
      LayerThicknessH[I] = HostArray2DR8("LayerThickness" + std::to_string(I),
                                         NCellsSize, NVertLevels);
      NormalVelocityH[I] = HostArray2DR8("NormalVelocity" + std::to_string(I),
                                         NEdgesSize, NVertLevels);
   }


   // Open the state file for reading (assume IO has already been initialized)
   I4 Err;
   Err = OMEGA::IO::openFile(StateFileID, StateFileName, IO::ModeRead);
   if (Err != 0)
      LOG_CRITICAL("OceanState: error opening state file");

   // Create the parallel IO decompositions required to read in state variables
   initParallelIO(MeshDecomp);

   // Read layerThickness and normalVelocity
   read();

   // Destroy the parallel IO decompositions
   finalizeParallelIO();

   // Register fields and metadata for IO
   defineIOFields();

   // Copy host data to device
   for (int I = 0; I < NTimeLevels; I++) {
      LayerThickness[I] = createDeviceMirrorCopy(LayerThicknessH[I]);
      NormalVelocity[I] = createDeviceMirrorCopy(NormalVelocityH[I]);  
   }

   // Associate this instance with a name
   AllOceanStates.emplace(Name, *this);

} // end state constructor

//------------------------------------------------------------------------------
// Destroys a local mesh and deallocates all arrays
OceanState::~OceanState() {

   // No operations needed, Kokkos arrays removed when no longer in scope

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
// Initialize the parallel IO decompositions for the mesh variables
void OceanState::initParallelIO(Decomp *MeshDecomp) {

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
void OceanState::defineIOFields() {

   int Err = 0;

   LayerThicknessIOName = "LayerThickness";
   NormalVelocityIOName = "NormalVelocity";
   std::string NCellsIOName = "NCells";
   std::string NEdgesIOName = "NEdges";
   std::string NVertLevelsIOName = "NVertLevels";
   if (Name != "Default") {
      LayerThicknessIOName.append(Name);
      NormalVelocityIOName.append(Name);
      NCellsIOName.append(Name);
      NEdgesIOName.append(Name);
      NVertLevelsIOName.append(Name);
   }


   // Create metadata dimensions
   auto CellDim = OMEGA::MetaDim::create(NCellsIOName, NCellsSize);
   auto EdgeDim = OMEGA::MetaDim::create(NEdgesIOName, NEdgesSize);
   auto VertDim = OMEGA::MetaDim::create(NVertLevelsIOName, NVertLevels);

   std::vector<std::shared_ptr<OMEGA::MetaDim>> LayerThicknessDim{CellDim,
                                                                  VertDim};
   std::vector<std::shared_ptr<OMEGA::MetaDim>> NormalVelocityDim{EdgeDim,
                                                                  VertDim};
   // Create metadate for variables
   auto NormalVelocityMeta = OMEGA::ArrayMetaData::create(
       NormalVelocityIOName,
       "Velocity component normal to edge", /// long Name
       "m/s",                               /// units
       "",                                  /// CF standard Name
       -9.99E+10,                           /// min valid value
       9.99E+10,                            /// max valid value
       -9.99E+30,        /// scalar used for undefined entries
       2,                /// number of dimensions
       NormalVelocityDim /// dim pointers
   );

   auto LayerThicknessMeta = OMEGA::ArrayMetaData::create(
       LayerThicknessIOName,
       "Thickness of layer on cell center", /// long Name
       "m",                                 /// units
       "cell_thickness",                    /// CF standard Name
       0.0,                                 /// min valid value
       9.99E+30,                            /// max valid value
       -9.99E+30,        /// scalar used for undefined entries
       2,                /// number of dimensions
       LayerThicknessDim /// dim pointers
   );

   // Group state metadata
   auto StateMetaGroup = OMEGA::MetaGroup::create("State");

   Err = StateMetaGroup->addField(NormalVelocityIOName);
   Err = StateMetaGroup->addField(LayerThicknessIOName);

   // Define IOFields for state variables
   Err = OMEGA::IOField::define(NormalVelocityIOName);
   Err = OMEGA::IOField::define(LayerThicknessIOName);

   // Associate IOField with data
   int CurLevel = NTimeLevels - 2;

   Err = OMEGA::IOField::attachData<OMEGA::Array2DR8>(NormalVelocityIOName,
                                                      NormalVelocity[CurLevel]);
   Err = OMEGA::IOField::attachData<OMEGA::Array2DR8>(LayerThicknessIOName,
                                                      LayerThickness[CurLevel]);
} // end defineIOFields

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
// Perform time level update
void OceanState::updateTimeLevels() {

   int NewLevel = NTimeLevels - 1;

   // Update time levels for layer thickness
   copyToHost(NewLevel);
   MeshHalo->exchangeFullArrayHalo(LayerThicknessH[NewLevel], OMEGA::OnCell);
   copyToDevice(NewLevel);

   Array2DR8 Temp;
   HostArray2DR8 TempH;

   for (int Level = 0; Level < NTimeLevels - 1; Level++) {
      std::swap(LayerThickness[Level + 1], LayerThickness[Level]);
      std::swap(LayerThicknessH[Level + 1], LayerThicknessH[Level]);
   }

   // Update time levels for normal velocity
   copyToHost(NewLevel);
   MeshHalo->exchangeFullArrayHalo(NormalVelocityH[NewLevel], OMEGA::OnEdge);
   copyToDevice(NewLevel);

   for (int Level = 0; Level < NTimeLevels - 1; Level++) {
      std::swap(NormalVelocity[Level + 1], NormalVelocity[Level]);
      std::swap(NormalVelocityH[Level + 1], NormalVelocityH[Level]);
   }

   // Update IOField data associations
   int Err      = 0;
   int CurLevel = NTimeLevels - 2;

   Err = OMEGA::IOField::attachData<OMEGA::Array2DR8>(NormalVelocityIOName,
                                                      NormalVelocity[CurLevel]);
   Err = OMEGA::IOField::attachData<OMEGA::Array2DR8>(LayerThicknessIOName,
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
