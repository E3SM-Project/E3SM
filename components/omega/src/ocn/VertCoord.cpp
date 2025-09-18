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
// Begin initialization of default vertical coordinate, requires prior
// initialization of Decomp.
void VertCoord::init1() {

   Decomp *DefDecomp = Decomp::getDefault();

   VertCoord::DefaultVertCoord = create("Default", DefDecomp);

} // end init1

//------------------------------------------------------------------------------
// Complete initialization of default vertical coordinate, requires prior
// initialization of HorzMesh.
void VertCoord::init2() {

   Config *OmegaConfig = Config::getOmegaConfig();

   DefaultVertCoord->completeSetup(OmegaConfig);

} // end init2

//------------------------------------------------------------------------------
// Construct a new VertCoord instance given a Decomp. New object is incomplete
// and completeSetup must be called afterwards
VertCoord::VertCoord(const std::string &Name_, //< [in] Name for new VertCoord
                     const Decomp *Decomp      //< [in] associated Decomp
) {

   // Store name suffix
   Name = Name_;

   // Retrieve mesh filename from Decomp
   MeshFileName = Decomp->MeshFileName;

   // Open the mesh file for reading (assume IO has already been initialized)
   IO::openFile(MeshFileID, MeshFileName, IO::ModeRead);

   // Set NVertLayers and NVertLayersP1 and create the vertical dimensions
   Error Err; // Error code
   I4 NVertLayersID;
   Err = IO::getDimFromFile(MeshFileID, "nVertLevels", NVertLayersID,
                            NVertLayers);
   if (!Err.isSuccess()) {
      LOG_WARN("VertCoord: error reading nVertLevels from mesh file, "
               "using NVertLayers = 1");
      NVertLayers = 1;
   }
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
   NEdgesAll   = Decomp->NEdgesAll;
   NEdgesSize  = Decomp->NEdgesSize;

   NVerticesOwned = Decomp->NVerticesOwned;
   NVerticesAll   = Decomp->NVerticesAll;
   NVerticesSize  = Decomp->NVerticesSize;
   VertexDegree   = Decomp->VertexDegree;

   // Retrieve connectivity arrays from HorzMesh
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

   // Make host copies for device arrays not being read from file
   PressureInterfaceH    = createHostMirrorCopy(PressureInterface);
   PressureMidH          = createHostMirrorCopy(PressureMid);
   ZInterfaceH           = createHostMirrorCopy(ZInterface);
   ZMidH                 = createHostMirrorCopy(ZMid);
   GeopotentialMidH      = createHostMirrorCopy(GeopotentialMid);
   LayerThicknessTargetH = createHostMirrorCopy(LayerThicknessTarget);
   RefLayerThicknessH    = createHostMirrorCopy(RefLayerThickness);

} // end constructor

//------------------------------------------------------------------------------
// Complete construction of new VertCoord instance
void VertCoord::completeSetup(Config *Options //< [in] configuration options
) {

   // Define field metadata
   defineFields();

   I4 FillValueI4     = -1;
   Real FillValueReal = -999._Real;

   deepCopy(MinLayerCell, FillValueI4);
   deepCopy(MaxLayerCell, FillValueI4);
   deepCopy(BottomDepth, FillValueReal);

   OMEGA_SCOPE(LocMinLayerCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, MaxLayerCell);
   OMEGA_SCOPE(LocBottomDepth, BottomDepth);

   // Fetch input stream and validate
   std::string StreamName = "InitialVertCoord";
   if (Name != "Default") {
      StreamName.append(Name);
   }

   auto VCoordStream = IOStream::get(StreamName);

   bool IsValidated = VCoordStream->validate();

   // Read InitialVertCoord stream
   Error Err; // Error code
   if (IsValidated) {
      Err = IOStream::read(StreamName);
      if (!Err.isSuccess()) {
         LOG_WARN("VertCoord: Error reading {} stream", StreamName);
         I4 Sum1 = 0;
         parallelReduce(
             {MinLayerCell.extent_int(0)},
             KOKKOS_LAMBDA(int I, int &Accum) { Accum += LocMinLayerCell(I); },
             Sum1);
         if (Sum1 < 0) {
            LOG_WARN("VertCoord: Error reading minLevelCell from {}, "
                     "using MinLayerCell = 0",
                     StreamName);
            deepCopy(MinLayerCell, 1);
         }
         I4 Sum2 = 0;
         parallelReduce(
             {MaxLayerCell.extent_int(0)},
             KOKKOS_LAMBDA(int I, int &Accum) { Accum += LocMaxLayerCell(I); },
             Sum2);
         if (Sum2 < 0) {
            LOG_WARN("VertCoord: Error reading maxLevelCell from {}, "
                     "using MaxLayerCell = NVertLayers - 1",
                     StreamName);
            deepCopy(MaxLayerCell, NVertLayers);
         }
         Real Sum3 = 0.;
         parallelReduce(
             {BottomDepth.extent_int(0)},
             KOKKOS_LAMBDA(int I, Real &Accum) { Accum += LocBottomDepth(I); },
             Sum3);
         if (Sum3 < 0.) {
            ABORT_ERROR("VertCoord: Error reading bottomDepth from {}",
                        StreamName);
         }
      }
   } else {
      ABORT_ERROR("Error validating IO stream {}", StreamName);
   }

   // Subtract 1 to convert to zero-based indexing
   parallelFor(
       {NCellsAll}, KOKKOS_LAMBDA(int ICell) {
          LocMinLayerCell(ICell) -= 1;
          LocMaxLayerCell(ICell) -= 1;
       });

   // Compute Edge and Vertex vertical ranges
   minMaxLayerEdge();
   minMaxLayerVertex();

   // Initialize movement weights
   initMovementWeights(Options);

   // Make host copies for device arrays read from mesh file
   MaxLayerCellH = createHostMirrorCopy(MaxLayerCell);
   MinLayerCellH = createHostMirrorCopy(MinLayerCell);
   BottomDepthH  = createHostMirrorCopy(BottomDepth);

   // Fetch reference desnity from Config
   Config TendConfig("Tendencies");
   Err.reset();
   Err += Options->get(TendConfig);
   CHECK_ERROR_ABORT(Err, "VertCoord: Tendencies group not found in Config");

   Err += TendConfig.get("Density0", Rho0);
   CHECK_ERROR_ABORT(Err, "VertCoord: Density0 not found in TendConfig");

} // end completeSetup

//------------------------------------------------------------------------------
// Calls the VertCoord constructor and places it in the AllVertCoords map
VertCoord *
VertCoord::create(const std::string &Name, // [in] name for new VertCoord
                  const Decomp *Decomp     // [in] associated Decomp
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
   auto *NewVertCoord = new VertCoord(Name, Decomp);
   AllVertCoords.emplace(Name, NewVertCoord);

   return NewVertCoord;
} // end create

//------------------------------------------------------------------------------
// Define IO fields and metadata
void VertCoord::defineFields() {

   // Set field names (append Name if not default)
   MinLayerCellFldName   = "MinLevelCell";
   MaxLayerCellFldName   = "MaxLevelCell";
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

} // end clear

//------------------------------------------------------------------------------
// Compute min and max layer indices for edges based on MinLayerCell and
// MaxLayerCell
void VertCoord::minMaxLayerEdge() {

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
       {NEdgesAll}, KOKKOS_LAMBDA(int IEdge) {
          I4 Lyr1;
          I4 Lyr2;
          const I4 ICell1 = LocCellsOnEdge(IEdge, 0);
          const I4 ICell2 = LocCellsOnEdge(IEdge, 1);
          Lyr1            = LocMaxLayerCell(ICell1) == -1 ? LocNVertLayersP1
                                                          : LocMinLayerCell(ICell1);
          Lyr2            = LocMaxLayerCell(ICell2) == -1 ? LocNVertLayersP1
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

   OMEGA_SCOPE(LocNEdgesAll, NEdgesAll);
   parallelFor(
       {1}, KOKKOS_LAMBDA(const int &) {
          LocMinLayerEdgeTop(LocNEdgesAll) = LocNVertLayersP1;
          LocMinLayerEdgeBot(LocNEdgesAll) = LocNVertLayersP1;
          LocMaxLayerEdgeTop(LocNEdgesAll) = -1;
          LocMaxLayerEdgeBot(LocNEdgesAll) = -1;
       });

   MinLayerEdgeTopH = createHostMirrorCopy(MinLayerEdgeTop);
   MinLayerEdgeBotH = createHostMirrorCopy(MinLayerEdgeBot);
   MaxLayerEdgeTopH = createHostMirrorCopy(MaxLayerEdgeTop);
   MaxLayerEdgeBotH = createHostMirrorCopy(MaxLayerEdgeBot);
} // end MinMaxLayerEdge

//------------------------------------------------------------------------------
// Compute min and max layer indices for vertices based on MinLayerCell and
// MaxLayerCell
void VertCoord::minMaxLayerVertex() {

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
       {NVerticesAll}, KOKKOS_LAMBDA(int IVertex) {
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
   OMEGA_SCOPE(LocNVerticesAll, NVerticesAll);
   parallelFor(
       {1}, KOKKOS_LAMBDA(const int &) {
          LocMinLayerVertexTop(LocNVerticesAll) = LocNVertLayersP1;
          LocMinLayerVertexBot(LocNVerticesAll) = LocNVertLayersP1;
          LocMaxLayerVertexTop(LocNVerticesAll) = -1;
          LocMaxLayerVertexBot(LocNVerticesAll) = -1;
       });

   MinLayerVertexTopH = createHostMirrorCopy(MinLayerVertexTop);
   MinLayerVertexBotH = createHostMirrorCopy(MinLayerVertexBot);
   MaxLayerVertexTopH = createHostMirrorCopy(MaxLayerVertexTop);
   MaxLayerVertexBotH = createHostMirrorCopy(MaxLayerVertexBot);
} // end MinMaxLayerVertex

//------------------------------------------------------------------------------
// Store VertCoordMovementWeights based on config choice
void VertCoord::initMovementWeights(
    Config *Options // [in] configuration options
) {

   Error Err; // default successful error code

   Config VCoordConfig("VertCoord");
   Err += Options->get(VCoordConfig);
   CHECK_ERROR_ABORT(Err, "VertCoord: VertCoord group not found in Config");

   std::string MovementWeightType;
   Err += VCoordConfig.get("MovementWeightType", MovementWeightType);
   CHECK_ERROR_ABORT(Err,
                     "VertCoord: MovementWeightType not found in VertCoord");

   VertCoordMovementWeights =
       Array1DReal("VertCoordMovementWeights", NVertLayers);

   OMEGA_SCOPE(LocVertCoordMovementWeights, VertCoordMovementWeights);
   if (MovementWeightType == "Fixed") {
      deepCopy(VertCoordMovementWeights, 0._Real);
      parallelFor(
          {1}, KOKKOS_LAMBDA(const int &) {
             LocVertCoordMovementWeights(0) = 1._Real;
          });
   } else if (MovementWeightType == "Uniform") {
      deepCopy(VertCoordMovementWeights, 1._Real);
   } else {
      ABORT_ERROR("VertCoord: Unknown MovementWeightType requested");
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

   OMEGA_SCOPE(LocRho0, Rho0);
   OMEGA_SCOPE(LocMinLayerCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, MaxLayerCell);
   OMEGA_SCOPE(LocPressInterf, PressureInterface);
   OMEGA_SCOPE(LocPressMid, PressureMid);

   const auto Policy = TeamPolicy(NCellsAll, OMEGA_TEAMSIZE, 1);
   Kokkos::parallel_for(
       "computePressure", Policy, KOKKOS_LAMBDA(const TeamMember &Member) {
          const I4 ICell = Member.league_rank();
          const I4 KMin  = LocMinLayerCell(ICell);
          const I4 KMax  = LocMaxLayerCell(ICell);
          const I4 Range = KMax - KMin + 1;

          LocPressInterf(ICell, KMin) = SurfacePressure(ICell);
          Kokkos::parallel_scan(TeamThreadRange(Member, Range),
                                [=](int K, Real &Accum, bool IsFinal) {
                                   const I4 KLyr  = K + KMin;
                                   Real Increment = Gravity * LocRho0 *
                                                    LayerThickness(ICell, KLyr);
                                   Accum += Increment;

                                   if (IsFinal) {
                                      LocPressInterf(ICell, KLyr + 1) =
                                          SurfacePressure(ICell) + Accum;
                                      LocPressMid(ICell, KLyr) =
                                          SurfacePressure(ICell) + Accum -
                                          0.5 * Increment;
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

   OMEGA_SCOPE(LocRho0, Rho0);
   OMEGA_SCOPE(LocMinLayerCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, MaxLayerCell);
   OMEGA_SCOPE(LocZInterf, ZInterface);
   OMEGA_SCOPE(LocZMid, ZMid);
   OMEGA_SCOPE(LocBotDepth, BottomDepth);

   const auto Policy = TeamPolicy(NCellsAll, OMEGA_TEAMSIZE, 1);
   Kokkos::parallel_for(
       "computeZHeight", Policy, KOKKOS_LAMBDA(const TeamMember &Member) {
          const I4 ICell = Member.league_rank();
          const I4 KMin  = LocMinLayerCell(ICell);
          const I4 KMax  = LocMaxLayerCell(ICell);
          const I4 Range = KMax - KMin + 1;

          LocZInterf(ICell, KMax + 1) = -LocBotDepth(ICell);
          Kokkos::parallel_scan(
              TeamThreadRange(Member, Range),
              [=](int K, Real &Accum, bool IsFinal) {
                 const I4 KLyr = KMax - K;
                 Real DZ       = LocRho0 * SpecVol(ICell, KLyr) *
                           LayerThickness(ICell, KLyr);
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

   Kokkos::parallel_for(
       "computeGeopotential", TeamPolicy(NCellsAll, OMEGA_TEAMSIZE),
       KOKKOS_LAMBDA(const TeamMember &Member) {
          const I4 ICell   = Member.league_rank();
          const I4 KMin    = LocMinLayerCell(ICell);
          const I4 KMax    = LocMaxLayerCell(ICell);
          const I4 KRange  = KMax - KMin + 1;
          const I4 NChunks = (KRange + VecLength - 1) / VecLength;
          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(Member, NChunks), [=](const int KChunk) {
                 const I4 KStart = KMin + KChunk * VecLength;
                 const I4 KEnd   = KStart + VecLength;

                 const I4 KLen =
                     KEnd > KMax + 1 ? KMax + 1 - KStart : VecLength;
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

   OMEGA_SCOPE(LocRho0, Rho0);
   OMEGA_SCOPE(LocMinLayerCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, MaxLayerCell);
   OMEGA_SCOPE(LocLayerThickTarget, LayerThicknessTarget);
   OMEGA_SCOPE(LocPressInterf, PressureInterface);
   OMEGA_SCOPE(LocRefLayerThick, RefLayerThickness);
   OMEGA_SCOPE(LocVertCoordMvmtWgts, VertCoordMovementWeights);

   Kokkos::parallel_for(
       "computeTargetThickness", TeamPolicy(NCellsAll, OMEGA_TEAMSIZE),
       KOKKOS_LAMBDA(const TeamMember &Member) {
          const I4 ICell = Member.league_rank();
          const I4 KMin  = LocMinLayerCell(ICell);
          const I4 KMax  = LocMaxLayerCell(ICell);

          Real Coeff =
              (LocPressInterf(ICell, KMax + 1) - LocPressInterf(ICell, KMin)) /
              (Gravity * LocRho0);

          Real SumWh   = 0;
          Real SumRefH = 0;
          Kokkos::parallel_reduce(
              Kokkos::TeamThreadRange(Member, KMin, KMax + 1),
              [=](const int K, Real &LocalWh, Real &LocalSum) {
                 const Real RefLayerThick = LocRefLayerThick(ICell, K);
                 LocalWh += LocVertCoordMvmtWgts(K) * RefLayerThick;
                 LocalSum += RefLayerThick;
              },
              SumWh, SumRefH);
          Coeff -= SumRefH;

          const I4 KRange  = KMax - KMin + 1;
          const I4 NChunks = (KRange + VecLength - 1) / VecLength;

          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(Member, NChunks), [=](const int KChunk) {
                 const I4 KStart = KMin + KChunk * VecLength;
                 const I4 KEnd   = KStart + VecLength;

                 const I4 KLen =
                     KEnd > KMax + 1 ? KMax + 1 - KStart : VecLength;
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
