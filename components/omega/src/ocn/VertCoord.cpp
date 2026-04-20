//===-- base/VertCoord.cpp - vertical coordinate ----------------*- C++ -*-===//
//
// The VertCoord class contains member variables related to the vertical mesh
// and methods for computing these variables in Omega. The class will also serve
// as the container for information defining the extent of the vertical
// dimension and the number of active vertical layers in each column.
//
//===----------------------------------------------------------------------===//

#include "VertCoord.h"
#include "Dimension.h"
#include "Field.h"
#include "GlobalConstants.h"
#include "IO.h"
#include "IOStream.h"
#include "OmegaKokkos.h"

#include <limits>

namespace OMEGA {

// create static class members
VertCoord *VertCoord::DefaultVertCoord = nullptr;
std::map<std::string, std::unique_ptr<VertCoord>> VertCoord::AllVertCoords;

//------------------------------------------------------------------------------
// Initialize the default vertical coordinate, requires prior initialization
// of Decomp and Halo. The optional arguments simplify the use of the VertCoord
// in some unit tests. If ReadStream is false, the InitialVertCoord stream will
// not be read during construction. If a value for NVertLayers is passed as
// argument, the dimension will not be read from the mesh file.
void VertCoord::init(
    const bool
        ReadStream, //< [in] optional argument to read stream, true by default
    const int
        InNVertLayers //< [in] optional argument to set NVertLayers explicitly
                      //< instead of reading dimension from mesh file
) {

   Decomp *DefDecomp = Decomp::getDefault();

   Halo *DefHalo = Halo::getDefault();

   Config *OmegaConfig = Config::getOmegaConfig();

   VertCoord::DefaultVertCoord = create("Default", DefDecomp, DefHalo,
                                        OmegaConfig, ReadStream, InNVertLayers);

} // end init

//------------------------------------------------------------------------------
// Construct a new VertCoord instance given a Decomp.
VertCoord::VertCoord(const std::string &Name_, //< [in] Name for new VertCoord
                     const Decomp *Decomp,     //< [in] associated Decomp
                     Halo *MeshHalo,           //< [in] mesh halo exchanger
                     Config *Options,          //< [in] configuration options
                     const bool ReadStream,    //< [in] logical to read stream
                     const int InNVertLayers   //< [in] int to set vertical dim
) {
   Error Err; // Error code

   // If ReadStream is true, a prescribed value for NVertLayers is not valid
   if (ReadStream == true and InNVertLayers != 0) {
      ABORT_ERROR("VertCoord: ReadStream is true but a value for NVertLayers "
                  "is explicitly provided, which is not a valid combination");
   }

   // Read Config for movement weight type, store in enum
   Config VCoordConfig("VertCoord");
   Err += Options->get(VCoordConfig);
   CHECK_ERROR_ABORT(Err, "VertCoord: VertCoord group not found in Config");

   std::string MovementWeightStr;
   Err += VCoordConfig.get("MovementWeightType", MovementWeightStr);
   CHECK_ERROR_ABORT(Err,
                     "VertCoord: MovementWeightType not found in VertCoord");

   if (MovementWeightStr == "Fixed") {
      MvmtWgtChoice = MovementWeightType::Fixed;
   } else if (MovementWeightStr == "Uniform") {
      MvmtWgtChoice = MovementWeightType::Uniform;
   } else {
      ABORT_ERROR("VertCoord: Unknown MovementWeightType requested");
   }

   // Store name suffix
   Name = Name_;

   // Retrieve mesh filename from Decomp
   MeshFileName = Decomp->MeshFileName;

   // If InNVertlayers is 0 (default), attempt to read the dimension from the
   // mesh file. If the mesh file does not define a vertical dimension, use
   // NVertLayers = 1. If a value for InNVertLayers is explicitly provided,
   // use that value is instead.
   if (InNVertLayers == 0) {

      std::string OmegaDimName = "NVertLayers";
      std::string MPASDimName  = "nVertLevels";

      // Open the mesh file for reading (assume IO has already been initialized)
      IO::openFileRead(MeshFileID, MeshFileName);

      // Set NVertLayers
      I4 NVertLayersID;
      Err = IO::getDimFromFile(MeshFileID, OmegaDimName, NVertLayersID,
                               NVertLayers);
      if (Err.isFail()) { // dim not found, try again with older MPAS name
         Err = IO::getDimFromFile(MeshFileID, MPASDimName, NVertLayersID,
                                  NVertLayers);
         if (Err.isFail()) {
            LOG_INFO("VertCoord: vertical dimension not found in mesh file, "
                     "using NVertLayers = 1");
            NVertLayers = 1;
         }
      }
   } else {
      NVertLayers = InNVertLayers;
   }

   // Create the dimensions
   NVertLayersP1 = NVertLayers + 1;

   std::string DimName   = "NVertLayers";
   std::string DimP1Name = "NVertLayersP1";
   if (Name != "Default") {
      DimName.append(Name);
      DimP1Name.append(Name);
   }

   auto VertDim   = Dimension::create(DimName, NVertLayers);
   auto VertDimP1 = Dimension::create(DimP1Name, NVertLayersP1);

   // Retrieve mesh variables from Decomp
   NCellsOwned = Decomp->NCellsOwned;
   NCellsAll   = Decomp->NCellsAll;
   NCellsSize  = Decomp->NCellsSize;

   NEdgesOwned = Decomp->NEdgesOwned;
   NEdgesHalo0 = Decomp->NEdgesHaloH(0);
   NEdgesAll   = Decomp->NEdgesAll;
   NEdgesSize  = Decomp->NEdgesSize;

   NVerticesOwned = Decomp->NVerticesOwned;
   NVerticesHalo0 = Decomp->NVerticesHaloH(0);
   NVerticesAll   = Decomp->NVerticesAll;
   NVerticesSize  = Decomp->NVerticesSize;
   VertexDegree   = Decomp->VertexDegree;

   // Retrieve connectivity arrays from Decomp
   CellsOnEdge   = Decomp->CellsOnEdge;
   CellsOnVertex = Decomp->CellsOnVertex;

   // Allocate device arrays
   MaxLayerCell = Array1DI4("MaxLayerCell", NCellsSize);
   MinLayerCell = Array1DI4("MinLayerCell", NCellsSize);
   BottomDepth  = Array1DReal("BottomDepth", NCellsSize);
   PressureInterface =
       Array2DReal("PressureInterface", NCellsSize, NVertLayersP1);
   PressureMid     = Array2DReal("PressureMid", NCellsSize, NVertLayers);
   ZInterface      = Array2DReal("ZInterface", NCellsSize, NVertLayersP1);
   ZMid            = Array2DReal("ZMid", NCellsSize, NVertLayers);
   GeopotentialMid = Array2DReal("GeopotentialMid", NCellsSize, NVertLayers);
   LayerThicknessTarget =
       Array2DReal("LayerThicknessTarget", NCellsSize, NVertLayers);
   RefLayerThickness =
       Array2DReal("RefLayerThickness", NCellsSize, NVertLayers);

   // TODO: Temporary handling of SurfacePressure
   SurfacePressure = Array1DReal("SurfacePressure", NCellsSize);
   deepCopy(SurfacePressure, 0);

   // Make host copies for device arrays not being read from file
   PressureInterfaceH    = createHostMirrorCopy(PressureInterface);
   PressureMidH          = createHostMirrorCopy(PressureMid);
   ZInterfaceH           = createHostMirrorCopy(ZInterface);
   ZMidH                 = createHostMirrorCopy(ZMid);
   GeopotentialMidH      = createHostMirrorCopy(GeopotentialMid);
   LayerThicknessTargetH = createHostMirrorCopy(LayerThicknessTarget);
   RefLayerThicknessH    = createHostMirrorCopy(RefLayerThickness);

   // Define field metadata
   defineFields();

   // Set MinLayerCell, MaxLayerCell, and BottomDepth arrays
   setStreamArrays(ReadStream, MeshHalo);

   // Compute Edge and Vertex vertical ranges
   minMaxLayerEdge(MeshHalo);
   minMaxLayerVertex(MeshHalo);

   // Set computational masks
   setMasks();

   // Initialize movement weights
   initMovementWeights();

} // end constructor

//------------------------------------------------------------------------------
// Calls the VertCoord constructor and places it in the AllVertCoords map
VertCoord *VertCoord::create(
    const std::string &Name, // [in] name for new VertCoord
    const Decomp *Decomp,    // [in] associated Decomp
    Halo *MeshHalo,          // [in] mesh halo exchanger
    Config *Options,         // [in] configuration options
    const bool ReadStream,   // [in] optional logical to read stream
    const int InNVertLayers  // [in] optional int to set vertical dim
) {
   // Check to see if a VertCoord of the same name already exists and, if so,
   // exit with an error
   if (AllVertCoords.find(Name) != AllVertCoords.end()) {
      LOG_ERROR("Attempted to create a VertCoord with name {} but a VertCoord "
                "of that name already exists",
                Name);
      return nullptr;
   }

   // create a new VertCoord on the heap and put it in a map of unique_ptrs,
   // which will manage its lifetime
   auto *NewVertCoord = new VertCoord(Name, Decomp, MeshHalo, Options,
                                      ReadStream, InNVertLayers);
   AllVertCoords.emplace(Name, NewVertCoord);

   return NewVertCoord;
} // end create

//------------------------------------------------------------------------------
// Define IO fields and metadata
void VertCoord::defineFields() {

   // Set field names (append Name if not default)
   MinLayerCellFldName   = "MinLayerCell";
   MaxLayerCellFldName   = "MaxLayerCell";
   BottomDepthFldName    = "BottomDepth";
   PressInterfFldName    = "PressureInterface";
   PressMidFldName       = "PressureMid";
   ZInterfFldName        = "ZInterface";
   ZMidFldName           = "ZMid";
   GeopotFldName         = "GeopotentialMid";
   LyrThickTargetFldName = "LayerThicknessTarget";

   if (Name != "Default") {
      MinLayerCellFldName.append(Name);
      MaxLayerCellFldName.append(Name);
      BottomDepthFldName.append(Name);
      PressInterfFldName.append(Name);
      PressMidFldName.append(Name);
      ZInterfFldName.append(Name);
      ZMidFldName.append(Name);
      GeopotFldName.append(Name);
      LyrThickTargetFldName.append(Name);
   }

   // Create fields for VertCoord variables
   const I4 FillValueI4     = -999;
   const Real FillValueReal = -9.99e30;
   int NDims                = 1;
   std::vector<std::string> DimNames(NDims);
   DimNames[0] = "NCells";

   auto MinLayerCellField = Field::create(
       MinLayerCellFldName,                         // field name
       "Index to first active cell in each column", // long name or description
       "",                                          // units
       "",                                          // CF standard Name
       0,                                           // min valid value
       std::numeric_limits<I4>::max(),              // max valid value
       FillValueI4, // scalar for undefined entries
       NDims,       // number of dimensions
       DimNames     // dimension names
   );

   auto MaxLayerCellField = Field::create(
       MaxLayerCellFldName,                        // field name
       "Index to last active cell in each column", // long name or description
       "",                                         // units
       "",                                         // CF standard Name
       -1,                                         // min valid value
       std::numeric_limits<I4>::max(),             // max valid value
       FillValueI4, // scalar for undefined entries
       NDims,       // number of dimensions
       DimNames     // dimension names
   );

   auto BottomDepthField = Field::create(
       BottomDepthFldName, // field name
       "Depth of the bottom of the ocean. Given as a positive distance from"
       "sea level",                      // long name or description
       "m",                              // units
       "sea_floor_depth_below_geoid",    // CF standard Name
       0.0,                              // min valid value
       std::numeric_limits<Real>::max(), // max valid value
       FillValueReal,                    // scalar for undefined entries
       NDims,                            // number of dimensions
       DimNames                          // dimension names
   );

   NDims = 2;
   DimNames.resize(NDims);
   DimNames[1] = "NVertLayersP1";

   auto PressureInterfaceField = Field::create(
       PressInterfFldName,                      // field name
       "Pressure at vertical layer interfaces", // long name or description
       "Pa",                                    // units
       "sea_water_pressure",                    // CF standard Name
       0.0,                                     // min valid value
       std::numeric_limits<Real>::max(),        // max valid value
       FillValueReal,                           // scalar for undefined entries
       NDims,                                   // number of dimensions
       DimNames                                 // dimension names
   );

   auto ZInterfaceField = Field::create(
       ZInterfFldName,                               // field name
       "Cartesian Z coordinate at layer interfaces", // long name or description
       "m",                                          // units
       "height",                                     // CF standard Name
       std::numeric_limits<Real>::min(),             // min valid value
       std::numeric_limits<Real>::max(),             // max valid value
       FillValueReal, // scalar for undefined entries
       NDims,         // number of dimensions
       DimNames       // dimension names
   );

   DimNames[1] = "NVertLayers";

   auto PressureMidField = Field::create(
       PressMidFldName,                        // field name
       "Pressure at vertical layer midpoints", // long name or description
       "Pa",                                   // units
       "sea_water_pressure",                   // CF standard Name
       0.0,                                    // min valid value
       std::numeric_limits<Real>::max(),       // max valid value
       FillValueReal,                          // scalar for undefined entries
       NDims,                                  // number of dimensions
       DimNames                                // dimension names
   );

   auto ZMidField = Field::create(
       ZMidFldName,                                 // field name
       "Cartesian Z coordinate at layer midpoints", // long name or description
       "m",                                         // units
       "height",                                    // CF standard Name
       std::numeric_limits<Real>::min(),            // min valid value
       std::numeric_limits<Real>::max(),            // max valid value
       FillValueReal, // scalar for undefined entries
       NDims,         // number of dimensions
       DimNames       // dimension names
   );

   auto GeopotentialMidField = Field::create(
       GeopotFldName,                     // field name
       "Geopotential at layer midpoints", // long name or description
       "m^2 s^-2",                        // units
       "geopotential",                    // CF standard Name
       std::numeric_limits<Real>::min(),  // min valid value
       std::numeric_limits<Real>::max(),  // max valid value
       FillValueReal,                     // scalar for undefined entries
       NDims,                             // number of dimensions
       DimNames                           // dimension names
   );

   auto LayerThicknessTargetField =
       Field::create(LyrThickTargetFldName, // field name
                     "desired layer thickness based on total perturbation from "
                     "the reference thickness", // long name or description
                     "m",                       // units
                     "",                        // CF standard Name
                     0.0,                       // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     FillValueReal, // scalar for undefined entries
                     NDims,         // number of dimensions
                     DimNames       // dimension names
       );

   // Create a field group for initial VertCoord fields
   InitGroupName = "InitVertCoord";
   if (Name != "Default") {
      InitGroupName.append(Name);
   }
   auto InitVCoordGroup = FieldGroup::create(InitGroupName);

   InitVCoordGroup->addField(MinLayerCellFldName);
   InitVCoordGroup->addField(MaxLayerCellFldName);
   InitVCoordGroup->addField(BottomDepthFldName);

   MinLayerCellField->attachData<Array1DI4>(MinLayerCell);
   MaxLayerCellField->attachData<Array1DI4>(MaxLayerCell);
   BottomDepthField->attachData<Array1DReal>(BottomDepth);

   // Create a field group for VertCoord fields
   GroupName = "VertCoord";
   if (Name != "Default") {
      GroupName.append(Name);
   }
   auto VCoordGroup = FieldGroup::create(GroupName);

   VCoordGroup->addField(PressInterfFldName);
   VCoordGroup->addField(PressMidFldName);
   VCoordGroup->addField(ZInterfFldName);
   VCoordGroup->addField(ZMidFldName);
   VCoordGroup->addField(GeopotFldName);
   VCoordGroup->addField(LyrThickTargetFldName);

   // Associate Field with data
   PressureInterfaceField->attachData<Array2DReal>(PressureInterface);
   PressureMidField->attachData<Array2DReal>(PressureMid);
   ZInterfaceField->attachData<Array2DReal>(ZInterface);
   ZMidField->attachData<Array2DReal>(ZMid);
   GeopotentialMidField->attachData<Array2DReal>(GeopotentialMid);
   LayerThicknessTargetField->attachData<Array2DReal>(LayerThicknessTarget);

} // end defineFields

//------------------------------------------------------------------------------
// Destroys a local VertCoord and deallocates all arrays
VertCoord::~VertCoord() {

   if (FieldGroup::exists(InitGroupName)) {
      Field::destroy(MinLayerCellFldName);
      Field::destroy(MaxLayerCellFldName);
      Field::destroy(BottomDepthFldName);
      FieldGroup::destroy(InitGroupName);
   }

   if (FieldGroup::exists(GroupName)) {
      Field::destroy(PressInterfFldName);
      Field::destroy(PressMidFldName);
      Field::destroy(ZInterfFldName);
      Field::destroy(ZMidFldName);
      Field::destroy(GeopotFldName);
      Field::destroy(LyrThickTargetFldName);
      FieldGroup::destroy(GroupName);
   }

} // end destructor

//------------------------------------------------------------------------------
// Removes a VertCoord from map by name
void VertCoord::erase(std::string Name) {
   AllVertCoords.erase(Name); // removes the VertCoord from the list and in the
                              // process, calls the destructor
} // end erase

//------------------------------------------------------------------------------
// Removes all VertCoords to clean up before exit
void VertCoord::clear() {

   AllVertCoords.clear(); // removes all VertCoords from the list and in the
                          // process, calls the destructors for each

   DefaultVertCoord = nullptr; // prevent dangling pointer
} // end clear

//------------------------------------------------------------------------------
// If ReadStream = true, read MinLayerCell, MaxLayerCell, and BottomDepth from
// the initial stream. If ReadStream = false, default values are used.
void VertCoord::setStreamArrays(const bool ReadStream, Halo *MeshHalo) {

   Error Err; // Error code

   OMEGA_SCOPE(LocMinLayerCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, MaxLayerCell);
   OMEGA_SCOPE(LocBottomDepth, BottomDepth);

   // If ReadStream is true (default) attempt to read values for MinLayerCell,
   // MaxLayerCell, and BottomDepth from the InitialVertCoord stream. Otherwise,
   // MinLayerCell and MaxLayerCell will be set to the first and last indices of
   // the vertical range, BottomDepth will remain uninitialized and will need
   // to be initialized explicitly if needed.
   if (ReadStream) {

      I4 FillValueI4     = -1;
      Real FillValueReal = -999._Real;

      deepCopy(MinLayerCell, FillValueI4);
      deepCopy(MaxLayerCell, FillValueI4);
      deepCopy(BottomDepth, FillValueReal);

      // Fetch input stream and validate
      std::string StreamName = "InitialVertCoord";
      if (Name != "Default") {
         StreamName.append(Name);
      }

      // Validate InitalVertCoord stream
      auto VCoordStream = IOStream::get(StreamName);
      bool IsValidated  = VCoordStream->validate();

      // Read InitialVertCoord stream
      if (IsValidated) {
         // Attempt to read stream, an error will be raised if any field fails
         // to be read. Determine which fields may not have been read properly.
         // If MinLayerCell or MaxLayerCell were not read properly default
         // values will be used. If BottomDepth was not read properly, abort
         // with error.
         Err = IOStream::read(StreamName);
         if (Err.isFail()) {
            LOG_INFO("VertCoord: Error while reading {} stream", StreamName);
            I4 Sum1 = 0;
            parallelReduce(
                {MinLayerCell.extent_int(0)},
                KOKKOS_LAMBDA(int I, int &Accum) {
                   Accum += LocMinLayerCell(I);
                },
                Sum1);
            if (Sum1 < 0) {
               LOG_INFO("VertCoord: Error reading MinLayerCell from {}, "
                        "using MinLayerCell = 0",
                        StreamName);
               deepCopy(MinLayerCell, 1);
            }
            I4 Sum2 = 0;
            parallelReduce(
                {MaxLayerCell.extent_int(0)},
                KOKKOS_LAMBDA(int I, int &Accum) {
                   Accum += LocMaxLayerCell(I);
                },
                Sum2);
            if (Sum2 < 0) {
               LOG_INFO("VertCoord: Error reading MaxLayerCell from {}, "
                        "using MaxLayerCell = NVertLayers - 1",
                        StreamName);
               deepCopy(MaxLayerCell, NVertLayers);
            }
            Real Sum3 = 0.;
            parallelReduce(
                {BottomDepth.extent_int(0)},
                KOKKOS_LAMBDA(int I, Real &Accum) {
                   Accum += LocBottomDepth(I);
                },
                Sum3);
            if (Sum3 < 0.) {
               ABORT_ERROR("VertCoord: Error reading bottomDepth from {}",
                           StreamName);
            }
         }
      } else {
         ABORT_ERROR("Error validating IO stream {}", StreamName);
      }
   } else {
      deepCopy(MinLayerCell, 1);
      deepCopy(MaxLayerCell, NVertLayers);
   }

   // Subtract 1 to convert to zero-based indexing
   parallelFor(
       {NCellsAll}, KOKKOS_LAMBDA(int ICell) {
          LocMinLayerCell(ICell) -= 1;
          LocMaxLayerCell(ICell) -= 1;
       });

   // Exchange halos since stream only reads owned cells
   MeshHalo->exchangeFullArrayHalo(MinLayerCell, OnCell);
   MeshHalo->exchangeFullArrayHalo(MaxLayerCell, OnCell);
   MeshHalo->exchangeFullArrayHalo(BottomDepth, OnCell);

   // The index ICell = NCellsAll represents an inactive cell
   OMEGA_SCOPE(LocNCellsAll, NCellsAll);
   OMEGA_SCOPE(LocNVertLayersP1, NVertLayersP1);
   parallelFor(
       {1}, KOKKOS_LAMBDA(const int &) {
          LocMinLayerCell(LocNCellsAll) = LocNVertLayersP1;
          LocMaxLayerCell(LocNCellsAll) = -1;
       });

   // Make host copies for device arrays read from mesh file
   MaxLayerCellH = createHostMirrorCopy(MaxLayerCell);
   MinLayerCellH = createHostMirrorCopy(MinLayerCell);
   BottomDepthH  = createHostMirrorCopy(BottomDepth);
}

//------------------------------------------------------------------------------
// Compute min and max layer indices for edges based on MinLayerCell and
// MaxLayerCell
void VertCoord::minMaxLayerEdge(Halo *MeshHalo) {

   MinLayerEdgeTop = Array1DI4("MinLayerEdgeTop", NEdgesSize);
   MinLayerEdgeBot = Array1DI4("MinLayerEdgeBot", NEdgesSize);
   MaxLayerEdgeTop = Array1DI4("MaxLayerEdgeTop", NEdgesSize);
   MaxLayerEdgeBot = Array1DI4("MaxLayerEdgeBot", NEdgesSize);

   OMEGA_SCOPE(LocNVertLayersP1, NVertLayersP1);
   OMEGA_SCOPE(LocCellsOnEdge, CellsOnEdge);
   OMEGA_SCOPE(LocMinLayerCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, MaxLayerCell);
   OMEGA_SCOPE(LocMinLayerEdgeTop, MinLayerEdgeTop);
   OMEGA_SCOPE(LocMinLayerEdgeBot, MinLayerEdgeBot);
   OMEGA_SCOPE(LocMaxLayerEdgeTop, MaxLayerEdgeTop);
   OMEGA_SCOPE(LocMaxLayerEdgeBot, MaxLayerEdgeBot);
   parallelFor(
       {NEdgesOwned}, KOKKOS_LAMBDA(int IEdge) {
          I4 Lyr1;
          I4 Lyr2;
          const I4 ICell1 = LocCellsOnEdge(IEdge, 0);
          const I4 ICell2 = LocCellsOnEdge(IEdge, 1);

          Lyr1 = LocMaxLayerCell(ICell1) == -1 ? LocNVertLayersP1
                                               : LocMinLayerCell(ICell1);
          Lyr2 = LocMaxLayerCell(ICell2) == -1 ? LocNVertLayersP1
                                               : LocMinLayerCell(ICell2);
          LocMinLayerEdgeTop(IEdge) = Kokkos::min(Lyr1, Lyr2);

          Lyr1 = LocMaxLayerCell(ICell1) == -1 ? 0 : LocMinLayerCell(ICell1);
          Lyr2 = LocMaxLayerCell(ICell2) == -1 ? 0 : LocMinLayerCell(ICell2);
          LocMinLayerEdgeBot(IEdge) = Kokkos::max(Lyr1, Lyr2);

          LocMaxLayerEdgeTop(IEdge) =
              Kokkos::min(LocMaxLayerCell(ICell1), LocMaxLayerCell(ICell2));
          LocMaxLayerEdgeBot(IEdge) =
              Kokkos::max(LocMaxLayerCell(ICell1), LocMaxLayerCell(ICell2));
       });

   MeshHalo->exchangeFullArrayHalo(MinLayerEdgeTop, OnEdge);
   MeshHalo->exchangeFullArrayHalo(MinLayerEdgeBot, OnEdge);
   MeshHalo->exchangeFullArrayHalo(MaxLayerEdgeTop, OnEdge);
   MeshHalo->exchangeFullArrayHalo(MaxLayerEdgeBot, OnEdge);

   OMEGA_SCOPE(LocNEdgesAll, NEdgesAll);
   parallelFor(
       {1}, KOKKOS_LAMBDA(const int &) {
          LocMinLayerEdgeTop(LocNEdgesAll) = LocNVertLayersP1;
          LocMinLayerEdgeBot(LocNEdgesAll) = LocNVertLayersP1;
          LocMaxLayerEdgeTop(LocNEdgesAll) = -1;
          LocMaxLayerEdgeBot(LocNEdgesAll) = -1;
       });

   // Set the [MinLayerEdgeBot, MaxLayerEdgeTop] range to be invalid
   // on the outermost edges of the halo layer.
   OMEGA_SCOPE(LocNEdgesHalo0, NEdgesHalo0);
   OMEGA_SCOPE(LocNCellsAll, NCellsAll);
   // start from NEdgesHalo(0) to exclude any edges of the owned cells
   parallelFor(
       {NEdgesAll - NEdgesHalo0}, KOKKOS_LAMBDA(int EdgeOffset) {
          const int IEdge = LocNEdgesHalo0 + EdgeOffset;

          // Is any cell on this edge an invalid (remote) cell ?
          if (LocCellsOnEdge(IEdge, 0) == LocNCellsAll ||
              LocCellsOnEdge(IEdge, 1) == LocNCellsAll) {
             LocMinLayerEdgeBot(IEdge) = LocNVertLayersP1;
             LocMaxLayerEdgeTop(IEdge) = -1;
          }
       });

   MinLayerEdgeTopH = createHostMirrorCopy(MinLayerEdgeTop);
   MinLayerEdgeBotH = createHostMirrorCopy(MinLayerEdgeBot);
   MaxLayerEdgeTopH = createHostMirrorCopy(MaxLayerEdgeTop);
   MaxLayerEdgeBotH = createHostMirrorCopy(MaxLayerEdgeBot);
} // end MinMaxLayerEdge

//------------------------------------------------------------------------------
// Compute min and max layer indices for vertices based on MinLayerCell and
// MaxLayerCell
void VertCoord::minMaxLayerVertex(Halo *MeshHalo) {

   MinLayerVertexTop = Array1DI4("MinLayerVertexTop", NVerticesSize);
   MinLayerVertexBot = Array1DI4("MinLayerVertexBot", NVerticesSize);
   MaxLayerVertexTop = Array1DI4("MaxLayerVertexTop", NVerticesSize);
   MaxLayerVertexBot = Array1DI4("MaxLayerVertexBot", NVerticesSize);

   OMEGA_SCOPE(LocNVertLayersP1, NVertLayersP1);
   OMEGA_SCOPE(LocVertexDegree, VertexDegree);
   OMEGA_SCOPE(LocCellsOnVertex, CellsOnVertex);
   OMEGA_SCOPE(LocMinLayerCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, MaxLayerCell);
   OMEGA_SCOPE(LocMinLayerVertexTop, MinLayerVertexTop);
   OMEGA_SCOPE(LocMinLayerVertexBot, MinLayerVertexBot);
   OMEGA_SCOPE(LocMaxLayerVertexTop, MaxLayerVertexTop);
   OMEGA_SCOPE(LocMaxLayerVertexBot, MaxLayerVertexBot);

   parallelFor(
       {NVerticesOwned}, KOKKOS_LAMBDA(int IVertex) {
          I4 Lyr;
          I4 ICell = LocCellsOnVertex(IVertex, 0);
          Lyr      = LocMaxLayerCell(ICell) == -1 ? 0 : LocMinLayerCell(ICell);
          LocMinLayerVertexBot(IVertex) = Lyr;
          for (int I = 1; I < LocVertexDegree; ++I) {
             ICell = LocCellsOnVertex(IVertex, I);
             Lyr   = LocMaxLayerCell(ICell) == -1 ? 0 : LocMinLayerCell(ICell);
             LocMinLayerVertexBot(IVertex) =
                 Kokkos::max(LocMinLayerVertexBot(IVertex), Lyr);
          }

          ICell = LocCellsOnVertex(IVertex, 0);
          Lyr   = LocMaxLayerCell(ICell) == -1 ? LocNVertLayersP1
                                               : LocMinLayerCell(ICell);
          LocMinLayerVertexTop(IVertex) = Lyr;
          for (int I = 1; I < LocVertexDegree; ++I) {
             ICell = LocCellsOnVertex(IVertex, I);
             Lyr   = LocMaxLayerCell(ICell) == -1 ? LocNVertLayersP1
                                                  : LocMinLayerCell(ICell);
             LocMinLayerVertexTop(IVertex) =
                 Kokkos::min(LocMinLayerVertexTop(IVertex), Lyr);
          }

          ICell                         = LocCellsOnVertex(IVertex, 0);
          LocMaxLayerVertexBot(IVertex) = LocMaxLayerCell(ICell);
          for (int I = 1; I < LocVertexDegree; ++I) {
             ICell                         = LocCellsOnVertex(IVertex, I);
             LocMaxLayerVertexBot(IVertex) = Kokkos::max(
                 LocMaxLayerVertexBot(IVertex), LocMaxLayerCell(ICell));
          }

          ICell                         = LocCellsOnVertex(IVertex, 0);
          LocMaxLayerVertexTop(IVertex) = LocMaxLayerCell(ICell);
          for (int I = 1; I < LocVertexDegree; ++I) {
             ICell                         = LocCellsOnVertex(IVertex, I);
             LocMaxLayerVertexTop(IVertex) = Kokkos::min(
                 LocMaxLayerVertexTop(IVertex), LocMaxLayerCell(ICell));
          }
       });

   MeshHalo->exchangeFullArrayHalo(MinLayerVertexTop, OnVertex);
   MeshHalo->exchangeFullArrayHalo(MinLayerVertexBot, OnVertex);
   MeshHalo->exchangeFullArrayHalo(MaxLayerVertexTop, OnVertex);
   MeshHalo->exchangeFullArrayHalo(MaxLayerVertexBot, OnVertex);

   OMEGA_SCOPE(LocNVerticesAll, NVerticesAll);
   parallelFor(
       {1}, KOKKOS_LAMBDA(const int &) {
          LocMinLayerVertexTop(LocNVerticesAll) = LocNVertLayersP1;
          LocMinLayerVertexBot(LocNVerticesAll) = LocNVertLayersP1;
          LocMaxLayerVertexTop(LocNVerticesAll) = -1;
          LocMaxLayerVertexBot(LocNVerticesAll) = -1;
       });

   // Set the [MinLayerVertexBot, MaxLayerVertexTop] range to be invalid
   // on the outermost vertices of the halo layer.
   OMEGA_SCOPE(LocNVerticesHalo0, NVerticesHalo0);
   OMEGA_SCOPE(LocNCellsAll, NCellsAll);
   // start from NVerticesHalo(0) to exclude any vertices of the owned cells
   parallelFor(
       {NVerticesAll - NVerticesHalo0}, KOKKOS_LAMBDA(int VertexOffset) {
          const int IVertex = LocNVerticesHalo0 + VertexOffset;

          // Is any cell on this vertex an invalid (remote) cell ?
          bool IsOuterMost = LocCellsOnVertex(IVertex, 0) == LocNCellsAll;
          for (int I = 1; I < LocVertexDegree; ++I) {
             IsOuterMost =
                 IsOuterMost || LocCellsOnVertex(IVertex, I) == LocNCellsAll;
          }

          if (IsOuterMost) {
             LocMinLayerVertexBot(IVertex) = LocNVertLayersP1;
             LocMaxLayerVertexTop(IVertex) = -1;
          }
       });

   MinLayerVertexTopH = createHostMirrorCopy(MinLayerVertexTop);
   MinLayerVertexBotH = createHostMirrorCopy(MinLayerVertexBot);
   MaxLayerVertexTopH = createHostMirrorCopy(MaxLayerVertexTop);
   MaxLayerVertexBotH = createHostMirrorCopy(MaxLayerVertexBot);
} // end MinMaxLayerVertex

//------------------------------------------------------------------------------
// set computational masks for mesh elements
void VertCoord::setMasks() {

   EdgeMask = Array2DReal("EdgeMask", NEdgesSize, NVertLayers);

   OMEGA_SCOPE(LocEdgeMask, EdgeMask);
   OMEGA_SCOPE(LocMinLyrEdgeBot, MinLayerEdgeBot);
   OMEGA_SCOPE(LocMaxLyrEdgeTop, MaxLayerEdgeTop);

   // EdgeMask = 1 if active layers on both sides, 0 otherwise.
   deepCopy(EdgeMask, 0.);
   parallelForOuter(
       {NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          const I4 KMin = LocMinLyrEdgeBot(IEdge);
          const I4 KMax = LocMaxLyrEdgeTop(IEdge);

          parallelForInner(
              Team, KMax - KMin + 1, INNER_LAMBDA(int K) {
                 I4 KLyr = KMin + K;

                 LocEdgeMask(IEdge, KLyr) = 1._Real;
              });
       });

   EdgeMaskH = createHostMirrorCopy(EdgeMask);

   CellMask = Array2DReal("CellMask", NCellsSize, NVertLayers);

   OMEGA_SCOPE(LocCellMask, CellMask);
   OMEGA_SCOPE(LocMinLyrCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLyrCell, MaxLayerCell);

   // CellMask = 1 in active layers, 0 otherwise.
   deepCopy(CellMask, 0.);
   parallelForOuter(
       {NCellsAll}, KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          const I4 KMin = LocMinLyrCell(ICell);
          const I4 KMax = LocMaxLyrCell(ICell);

          parallelForInner(
              Team, KMax - KMin + 1, INNER_LAMBDA(int K) {
                 I4 KLyr = KMin + K;

                 LocCellMask(ICell, KLyr) = 1._Real;
              });
       });

   CellMaskH = createHostMirrorCopy(CellMask);

   VertexMask = Array2DReal("VertexMask", NVerticesSize, NVertLayers);

   OMEGA_SCOPE(LocVrtxMask, VertexMask);
   OMEGA_SCOPE(LocMinLyrVrtxTop, MinLayerVertexTop);
   OMEGA_SCOPE(LocMaxLyrVrtxBot, MaxLayerVertexBot);

   // VertexMask = 1 if at least 1 surrounding cell layer is active,
   // 0 otherwise.
   deepCopy(VertexMask, 0.);
   parallelForOuter(
       {NVerticesAll}, KOKKOS_LAMBDA(int IVertex, const TeamMember &Team) {
          const I4 KMin = LocMinLyrVrtxTop(IVertex);
          const I4 KMax = LocMaxLyrVrtxBot(IVertex);

          parallelForInner(
              Team, KMax - KMin + 1, INNER_LAMBDA(int K) {
                 I4 KLyr = KMin + K;

                 LocVrtxMask(IVertex, KLyr) = 1._Real;
              });
       });

   VertexMaskH = createHostMirrorCopy(VertexMask);

} // end setMasks()

//------------------------------------------------------------------------------
// Store VertCoordMovementWeights based on config choice
void VertCoord::initMovementWeights() {

   Error Err; // default successful error code

   VertCoordMovementWeights =
       Array1DReal("VertCoordMovementWeights", NVertLayers);

   OMEGA_SCOPE(LocVertCoordMovementWeights, VertCoordMovementWeights);
   if (MvmtWgtChoice == MovementWeightType::Fixed) {
      deepCopy(VertCoordMovementWeights, 0._Real);
      parallelFor(
          {1}, KOKKOS_LAMBDA(const int &) {
             LocVertCoordMovementWeights(0) = 1._Real;
          });
   } else if (MvmtWgtChoice == MovementWeightType::Uniform) {
      deepCopy(VertCoordMovementWeights, 1._Real);
   }

   VertCoordMovementWeightsH = createHostMirrorCopy(VertCoordMovementWeights);
}

//------------------------------------------------------------------------------
// Compute the pressure at each layer interface and midpoint given the
// LayerThickness and SurfacePressure. Hierarchical parallelism is used with a
// parallel_for loop over all cells and a parallel_scan performing a prefix sum
// in each column to compute pressure from the top-most active layer to the
// bottom-most active layer.
void VertCoord::computePressure(
    const Array2DReal &LayerThickness, // [in] pseudo thickness
    const Array1DReal &SurfacePressure // [in] surface pressure
) {

   OMEGA_SCOPE(LocMinLayerCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, MaxLayerCell);
   OMEGA_SCOPE(LocPressInterf, PressureInterface);
   OMEGA_SCOPE(LocPressMid, PressureMid);

   parallelForOuter(
       "computePressure", {NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          const I4 KMin   = LocMinLayerCell(ICell);
          const I4 KMax   = LocMaxLayerCell(ICell);
          const I4 KRange = vertRange(KMin, KMax);

          LocPressInterf(ICell, KMin) = SurfacePressure(ICell);
          parallelScanInner(
              Team, KRange, INNER_LAMBDA(int K, Real &Accum, bool IsFinal) {
                 const I4 KLyr  = K + KMin;
                 Real Increment = Gravity * RhoSw * LayerThickness(ICell, KLyr);
                 Accum += Increment;

                 if (IsFinal) {
                    LocPressInterf(ICell, KLyr + 1) =
                        SurfacePressure(ICell) + Accum;
                    LocPressMid(ICell, KLyr) =
                        SurfacePressure(ICell) + Accum - 0.5 * Increment;
                 }
              });
       });
} // end computePressure

//------------------------------------------------------------------------------
// Compute geometric height z at layer interfaces and midpoints given the
// LayerThickness, SpecVol, and BottomDepth. Hierarchical parallelism is used
// with a parallel_for loop over cells and a parallel_scan performing a prefix
// sum in each column to compute z from the bottom-most active layer to the
// top-most active layer
void VertCoord::computeZHeight(
    const Array2DReal &LayerThickness, // [in] pseudo thickness
    const Array2DReal &SpecVol         // [in] specific volume
) {

   OMEGA_SCOPE(LocMinLayerCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, MaxLayerCell);
   OMEGA_SCOPE(LocZInterf, ZInterface);
   OMEGA_SCOPE(LocZMid, ZMid);
   OMEGA_SCOPE(LocBotDepth, BottomDepth);

   parallelForOuter(
       "computeZHeight", {NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          const I4 KMin   = LocMinLayerCell(ICell);
          const I4 KMax   = LocMaxLayerCell(ICell);
          const I4 KRange = vertRange(KMin, KMax);

          LocZInterf(ICell, KMax + 1) = -LocBotDepth(ICell);
          parallelScanInner(
              Team, KRange, INNER_LAMBDA(int K, Real &Accum, bool IsFinal) {
                 const I4 KLyr = KMax - K;
                 Real DZ =
                     RhoSw * SpecVol(ICell, KLyr) * LayerThickness(ICell, KLyr);
                 Accum += DZ;
                 if (IsFinal) {
                    LocZInterf(ICell, KLyr) = -LocBotDepth(ICell) + Accum;
                    LocZMid(ICell, KLyr) =
                        -LocBotDepth(ICell) + Accum - 0.5 * DZ;
                 }
              });
       });
} // end computeZHeight

//------------------------------------------------------------------------------
// Compute geopotential given Zmid, TidalPotential, and SelfAttractionLoading.
// Nested parallel_fors loop over all cells and all active layers in a column to
// compute the geopotential at the midpoint of each layer. The tidal potential
// and SAL are configurable, default-off features. When off these arrays will
// just be zeroes.
void VertCoord::computeGeopotential(
    const Array1DReal &TidalPotential,       // [in] tidal potential
    const Array1DReal &SelfAttractionLoading // [in] self attraction and loading
) {

   OMEGA_SCOPE(LocMinLayerCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, MaxLayerCell);
   OMEGA_SCOPE(LocGeopotMid, GeopotentialMid);
   OMEGA_SCOPE(LocZMid, ZMid);

   parallelForOuter(
       "computeGeopotential", {NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          const I4 KMin    = LocMinLayerCell(ICell);
          const I4 KMax    = LocMaxLayerCell(ICell);
          const I4 NChunks = vertRangeChunked(KMin, KMax);
          parallelForInner(
              Team, NChunks, INNER_LAMBDA(const int KChunk) {
                 const I4 KStart = chunkStart(KChunk, KMin);
                 const I4 KLen   = chunkLength(KChunk, KStart, KMax);
                 for (int KVec = 0; KVec < KLen; ++KVec) {
                    const I4 K             = KStart + KVec;
                    LocGeopotMid(ICell, K) = Gravity * LocZMid(ICell, K) +
                                             TidalPotential(ICell) +
                                             SelfAttractionLoading(ICell);
                 }
              });
       });
} // end compute Geopotential

//------------------------------------------------------------------------------
// Compute the desired target thickness, given PressureInterface,
// RefLayerThickness, and VertCoordMovementWeights. Hierarchical parallelsim is
// used with an outer parallel_for loop over cells, and 2 paralel_reduce
// reductions and a parallel_for over the active layers within a column.
void VertCoord::computeTargetThickness() {

   OMEGA_SCOPE(LocMinLayerCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, MaxLayerCell);
   OMEGA_SCOPE(LocLayerThickTarget, LayerThicknessTarget);
   OMEGA_SCOPE(LocPressInterf, PressureInterface);
   OMEGA_SCOPE(LocRefLayerThick, RefLayerThickness);
   OMEGA_SCOPE(LocVertCoordMvmtWgts, VertCoordMovementWeights);

   parallelForOuter(
       "computeTargetThickness", {NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          const I4 KMin   = LocMinLayerCell(ICell);
          const I4 KMax   = LocMaxLayerCell(ICell);
          const I4 KRange = vertRange(KMin, KMax);

          Real Coeff =
              (LocPressInterf(ICell, KMax + 1) - LocPressInterf(ICell, KMin)) /
              (Gravity * RhoSw);

          Real SumWh   = 0;
          Real SumRefH = 0;
          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(const int K, Real &LocalWh, Real &LocalSum) {
                 const I4 KLyr            = K + KMin;
                 const Real RefLayerThick = LocRefLayerThick(ICell, KLyr);
                 LocalWh += LocVertCoordMvmtWgts(KLyr) * RefLayerThick;
                 LocalSum += RefLayerThick;
              },
              SumWh, SumRefH);
          Coeff -= SumRefH;

          const I4 NChunks = vertRangeChunked(KMin, KMax);

          parallelForInner(
              Team, NChunks, INNER_LAMBDA(const int KChunk) {
                 const I4 KStart = chunkStart(KChunk, KMin);
                 const I4 KLen   = chunkLength(KChunk, KStart, KMax);
                 for (int KVec = 0; KVec < KLen; ++KVec) {
                    const I4 K = KStart + KVec;
                    LocLayerThickTarget(ICell, K) =
                        LocRefLayerThick(ICell, K) *
                        (1._Real + Coeff * LocVertCoordMvmtWgts(K) / SumWh);
                 }
              });
       });
}

//------------------------------------------------------------------------------
// Perform deepCopy for each variable array from device to host
void VertCoord::copyToHost() {

   deepCopy(PressureInterfaceH, PressureInterface);
   deepCopy(PressureMidH, PressureMid);
   deepCopy(ZInterfaceH, ZInterface);
   deepCopy(ZMidH, ZMid);
   deepCopy(GeopotentialMidH, GeopotentialMid);
   deepCopy(LayerThicknessTargetH, LayerThicknessTarget);
   deepCopy(RefLayerThicknessH, RefLayerThickness);
}

//------------------------------------------------------------------------------
// Perform deepCopy for each variable array from host to device
void VertCoord::copyToDevice() {

   deepCopy(PressureInterface, PressureInterfaceH);
   deepCopy(PressureMid, PressureMidH);
   deepCopy(ZInterface, ZInterfaceH);
   deepCopy(ZMid, ZMidH);
   deepCopy(GeopotentialMid, GeopotentialMidH);
   deepCopy(LayerThicknessTarget, LayerThicknessTargetH);
   deepCopy(RefLayerThickness, RefLayerThicknessH);
}

//------------------------------------------------------------------------------
// Get default VertCoord
VertCoord *VertCoord::getDefault() { return VertCoord::DefaultVertCoord; }

//------------------------------------------------------------------------------
// Get VertCoord by name
VertCoord *VertCoord::get(const std::string Name ///< [in] Name of VertCoord
) {

   // look for an instance of this name
   auto it = AllVertCoords.find(Name);

   // if found, return the VertCoord pointer
   if (it != AllVertCoords.end()) {
      return it->second.get();

      // otherwise print error and return null pointer
   } else {
      LOG_ERROR("VertCoord::get: Attempt to retrieve non-existant VertCoord:");
      LOG_ERROR("{} has not been defined or has been removed", Name);
      return nullptr;
   }

} // end get VertCoord

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
