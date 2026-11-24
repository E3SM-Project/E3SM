//===-- base/HorzMesh.cpp - horizontal mesh methods -------------*- C++ -*-===//
//
// The mesh (Mesh) class initializes a local mesh domain based on a given
// domain decomposition. It retrieves the mesh count and connectivity
// information from the Decomp object and reads in all other mesh variables
// from the mesh file. It also manages the device copies of the mesh data.
// It is meant to provide a container for passing mesh variables throughout
// the OMEGA tendency computation routines.
//
//===----------------------------------------------------------------------===//

#include "HorzMesh.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Error.h"
#include "Field.h"
#include "GlobalConstants.h"
#include "Halo.h"
#include "IOStream.h"
#include "Logging.h"
#include "OmegaKokkos.h"
#include "Reductions.h"

#include <algorithm>
#include <limits>
#include <vector>

namespace OMEGA {

// create the static class members
HorzMesh *HorzMesh::DefaultHorzMesh = nullptr;
std::map<std::string, std::unique_ptr<HorzMesh>> HorzMesh::AllHorzMeshes;

//------------------------------------------------------------------------------
// Initialize the mesh. Assumes that Decomp and VertCoord have already been
// initialized.

void HorzMesh::init(const Clock *ModelClock //< [in] Model clock for IO alarms
) {

   OMEGA_REQUIRE(ModelClock, "Null ModelClock pointer in HorzMesh::init");

   // Retrieve the default decomposition
   Decomp *DefDecomp = Decomp::getDefault();
   OMEGA_REQUIRE(DefDecomp, "Null default Decomp pointer in HorzMesh::init");

   // Create the default mesh and set pointer to it
   HorzMesh::DefaultHorzMesh = create("Default", DefDecomp, ModelClock);
}

//------------------------------------------------------------------------------
// Construct a new local mesh given a decomposition

HorzMesh::HorzMesh(const std::string &Name, //< [in] Name for new mesh
                   Decomp *MeshDecomp,      //< [in] Decomp for the new mesh
                   const Clock *ModelClock  //< [in] Model clock for IO alarms
) {

   // Set mesh name based on input
   MeshName = Name;

   // Mesh filename should be the same as that used for the decomposition
   MeshFileName = MeshDecomp->MeshFileName;

   // Retrieve mesh cell/edge/vertex totals from Decomp
   NCellsHalo   = MeshDecomp->NCellsHalo;
   NCellsHaloH  = MeshDecomp->NCellsHaloH;
   NCellsOwned  = MeshDecomp->NCellsOwned;
   NCellsAll    = MeshDecomp->NCellsAll;
   NCellsSize   = MeshDecomp->NCellsSize;
   NCellsGlobal = MeshDecomp->NCellsGlobal;

   NEdgesHalo     = MeshDecomp->NEdgesHalo;
   NEdgesHaloH    = MeshDecomp->NEdgesHaloH;
   NEdgesOwned    = MeshDecomp->NEdgesOwned;
   NEdgesAll      = MeshDecomp->NEdgesAll;
   NEdgesSize     = MeshDecomp->NEdgesSize;
   MaxCellsOnEdge = MeshDecomp->MaxCellsOnEdge;
   MaxEdges       = MeshDecomp->MaxEdges;
   MaxEdges2      = MeshDecomp->MaxEdges2;
   NEdgesGlobal   = MeshDecomp->NEdgesGlobal;

   OMEGA_REQUIRE(MaxEdges <= MaxEdgesBound,
                 "HorzMesh: MaxEdges={} exceeds compile-time bound {}",
                 MaxEdges, MaxEdgesBound);

   NVerticesHalo  = MeshDecomp->NVerticesHalo;
   NVerticesHaloH = MeshDecomp->NVerticesHaloH;
   NVerticesOwned = MeshDecomp->NVerticesOwned;
   NVerticesAll   = MeshDecomp->NVerticesAll;
   NVerticesSize  = MeshDecomp->NVerticesSize;
   VertexDegree   = MeshDecomp->VertexDegree;

   // Retrieve connectivity arrays from Decomp

   CellID = MeshDecomp->CellID;

   CellsOnCellH       = MeshDecomp->CellsOnCellH;
   EdgesOnCellH       = MeshDecomp->EdgesOnCellH;
   NEdgesOnCellH      = MeshDecomp->NEdgesOnCellH;
   VerticesOnCellH    = MeshDecomp->VerticesOnCellH;
   CellsOnEdgeH       = MeshDecomp->CellsOnEdgeH;
   EdgesOnEdgeH       = MeshDecomp->EdgesOnEdgeH;
   NEdgesOnEdgeH      = MeshDecomp->NEdgesOnEdgeH;
   VerticesOnEdgeH    = MeshDecomp->VerticesOnEdgeH;
   CellsOnVertexH     = MeshDecomp->CellsOnVertexH;
   EdgesOnVertexH     = MeshDecomp->EdgesOnVertexH;
   NEdgesReconOnCellH = MeshDecomp->NEdgesReconOnCellH;
   ReconStencilCellH  = MeshDecomp->ReconStencilCellH;

   CellsOnCell       = MeshDecomp->CellsOnCell;
   EdgesOnCell       = MeshDecomp->EdgesOnCell;
   NEdgesOnCell      = MeshDecomp->NEdgesOnCell;
   VerticesOnCell    = MeshDecomp->VerticesOnCell;
   CellsOnEdge       = MeshDecomp->CellsOnEdge;
   EdgesOnEdge       = MeshDecomp->EdgesOnEdge;
   NEdgesOnEdge      = MeshDecomp->NEdgesOnEdge;
   VerticesOnEdge    = MeshDecomp->VerticesOnEdge;
   CellsOnVertex     = MeshDecomp->CellsOnVertex;
   EdgesOnVertex     = MeshDecomp->EdgesOnVertex;
   NEdgesReconOnCell = MeshDecomp->NEdgesReconOnCell;
   ReconStencilCell  = MeshDecomp->ReconStencilCell;

   // OnSphere is needed early (before the mesh field definitions below)
   // to decide whether to define the ReconWeightsCell field, so we
   // use the value Decomp already read from the mesh file. It is
   // redetermined below (with the rest of the sphere/plane attributes)
   // once the full mesh stream is read, which is redundant but harmless
   // since both come from the same file.
   OnSphere = MeshDecomp->OnSphere;

   // Create Omega Dimensions associated with this mesh
   createDimensions(MeshDecomp);

   // Allocate remaining device arrays and define all mesh fields for
   // associated I/O and initialization. Equivalent host arrays will be
   // allocated and initialized after the fields are filled on the device.
   // This must happen before the mesh stream is read below, since the
   // read populates the arrays attached to these fields.
   defineMeshFields();

   // Read the input mesh stream to fill most of the fields
   // Some global mesh attributes are needed - we define both current
   // and older MPAS names for back compatibility
   std::string UnknownStr = "unknown";
   Metadata ReqMetaData{{"OnSphere", UnknownStr},   {"on_a_sphere", UnknownStr},
                        {"IsPeriodic", UnknownStr}, {"is_periodic", UnknownStr},
                        {"SphereRadius", 0.0},      {"sphere_radius", 0.0},
                        {"XPeriod", 0.0},           {"YPeriod", 0.0},
                        {"x_period", 0.0},          {"y_period", 0.0}};

   Error Err = IOStream::read("HorzMeshIn", ModelClock, ReqMetaData);
   CHECK_ERROR_ABORT(Err, "HorzMesh: error reading input mesh stream");

   // Extract global attributes into mesh variables
   // Change strings to lower case for easier comparison

   std::string OnSphereStr =
       std::any_cast<std::string>(ReqMetaData["OnSphere"]);
   std::string OnSphereOld =
       std::any_cast<std::string>(ReqMetaData["on_a_sphere"]);
   std::transform(OnSphereStr.begin(), OnSphereStr.end(), OnSphereStr.begin(),
                  [](unsigned char c) { return std::tolower(c); });
   std::transform(OnSphereOld.begin(), OnSphereOld.end(), OnSphereOld.begin(),
                  [](unsigned char c) { return std::tolower(c); });

   std::string IsPeriodStr =
       std::any_cast<std::string>(ReqMetaData["IsPeriodic"]);
   std::string IsPeriodOld =
       std::any_cast<std::string>(ReqMetaData["is_periodic"]);
   std::transform(IsPeriodStr.begin(), IsPeriodStr.end(), IsPeriodStr.begin(),
                  [](unsigned char c) { return std::tolower(c); });
   std::transform(IsPeriodOld.begin(), IsPeriodOld.end(), IsPeriodOld.begin(),
                  [](unsigned char c) { return std::tolower(c); });

   R8 NewRadius = std::any_cast<R8>(ReqMetaData["SphereRadius"]);
   R8 OldRadius = std::any_cast<R8>(ReqMetaData["sphere_radius"]);
   R8 TmpRadius = std::max(NewRadius, OldRadius);

   R8 NewXPeriod = std::any_cast<R8>(ReqMetaData["XPeriod"]);
   R8 OldXPeriod = std::any_cast<R8>(ReqMetaData["x_period"]);
   R8 TmpXPeriod = std::max(NewXPeriod, OldXPeriod);

   R8 NewYPeriod = std::any_cast<R8>(ReqMetaData["YPeriod"]);
   R8 OldYPeriod = std::any_cast<R8>(ReqMetaData["y_period"]);
   R8 TmpYPeriod = std::max(NewYPeriod, OldYPeriod);

   if (OnSphereStr == "yes" or OnSphereOld == "yes") { // mesh on sphere
      OnSphere = true;
      OnPlane  = false;
      if (TmpRadius > 0.0) {
         SphereRadius = static_cast<Real>(TmpRadius);
      } else {
         ABORT_ERROR("Mesh is on sphere but sphere radius either missing or 0");
      }
      // The mesh must have been generated on the same sphere as the Earth
      // radius from the Physical Constants Dictionary. Compare in double
      // precision (before the cast to Real) so that the tolerance does not
      // depend on the build precision; the tolerance leaves room for a mesh
      // file that stores the radius with fewer significant digits.
      constexpr R8 RadiusRelTol = 1.e-8;
      OMEGA_REQUIRE(std::abs((pcd::mean_radius - TmpRadius) /
                             pcd::mean_radius) <= RadiusRelTol,
                    "Input mesh sphere radius {} does not match the Earth "
                    "radius {} from the Physical Constants Dictionary "
                    "(relative tolerance {})",
                    TmpRadius, pcd::mean_radius, RadiusRelTol);
      IsPeriodic = false;
      XPeriod    = 0.0;
      YPeriod    = 0.0;
   } else if (OnSphereStr == "no" or OnSphereOld == "no") { // planar mesh
      OnSphere     = false;
      OnPlane      = true;
      SphereRadius = 0.0;
      if (IsPeriodStr == "yes" or IsPeriodOld == "yes") { // periodic mesh
         IsPeriodic = true;
         if (TmpXPeriod > 0.0 or TmpYPeriod > 0.0) {
            XPeriod = static_cast<Real>(TmpXPeriod);
            YPeriod = static_cast<Real>(TmpYPeriod);
         } else {
            ABORT_ERROR(
                "Mesh is periodic but periodicity lengths missing or invalid");
         }

      } else { // non-periodic (also assumed if no is periodic in mesh file)
         IsPeriodic = false;
         XPeriod    = 0.0;
         YPeriod    = 0.0;
      }
   } else {
      ABORT_ERROR("OnSphere or on_a_sphere missing from mesh file");
   }

   // Complete read arrays (fill halos and duplicate on host)
   completeReadArrays();

   // Compute additional mesh quantities
   // Compute EdgeSignOnCells and EdgeSignOnVertex
   computeEdgeSign();

   // Compute mesh scaling coefficients
   computeMeshScaling();

} // end horizontal mesh constructor

//------------------------------------------------------------------------------
// Creates a new mesh by calling the constructor and adding to the map of
// all horizontal meshes
HorzMesh *HorzMesh::create(const std::string &Name, //< [in] Name for new mesh
                           Decomp *MeshDecomp,      //< [in] Decomp for new mesh
                           const Clock *ModelClock  //< [in] Model clock for IO
) {
   OMEGA_REQUIRE(MeshDecomp,
                 "Null Decomp pointer in HorzMesh::create with Name = {}",
                 Name);
   OMEGA_REQUIRE(ModelClock,
                 "Null ModelClock pointer in HorzMesh::create with Name = {}",
                 Name);

   // Check to see if a mesh of the same name already exists and
   // if so, exit with an error
   if (AllHorzMeshes.find(Name) != AllHorzMeshes.end())
      ABORT_ERROR("Attempted to create a HorzMesh with name {} but a HorzMesh "
                  "of that name already exists",
                  Name);

   // create a new mesh on the heap and put it in a map of
   // unique_ptrs, which will manage its lifetime
   auto *NewHorzMesh = new HorzMesh(Name, MeshDecomp, ModelClock);
   AllHorzMeshes.emplace(Name, NewHorzMesh);

   return NewHorzMesh;
}

//------------------------------------------------------------------------------
// Destroys a local mesh and deallocates all arrays
HorzMesh::~HorzMesh() {

   // No operations needed, Kokkos arrays removed when no longer in scope

} // end deconstructor

//------------------------------------------------------------------------------
// Removes a mesh from list by name
void HorzMesh::erase(std::string InName // [in] name of mesh to remove
) {

   AllHorzMeshes.erase(InName); // remove the mesh from the list and in
                                // the process, calls the destructor

} // end mesh erase
//------------------------------------------------------------------------------
// Removes all meshes to clean up before exit
void HorzMesh::clear() {

   AllHorzMeshes.clear(); // removes all meshes from the list and in
                          // the process, calls the destructors for each

   DefaultHorzMesh = nullptr; // prevent dangling pointer
} // end clear

//------------------------------------------------------------------------------
// Create Dimension instances for all dimensions associated with the mesh
void HorzMesh::createDimensions(Decomp *MeshDecomp) {

   // Create non-distributed dimensions
   // If these have already been created, skip since the dims should not change
   if (!Dimension::exists("MaxCellsOnEdge"))
      auto MaxCellsOnEdgeDim =
          Dimension::create("MaxCellsOnEdge", MaxCellsOnEdge);
   if (!Dimension::exists("MaxEdges"))
      auto MaxEdgesDim = Dimension::create("MaxEdges", MaxEdges);
   if (!Dimension::exists("MaxEdges2"))
      auto MaxEdgesDim = Dimension::create("MaxEdges2", MaxEdges2);
   if (!Dimension::exists("VertexDegree"))
      auto VertexDegreeDim = Dimension::create("VertexDegree", VertexDegree);
   if (!Dimension::exists("R3"))
      auto R3Dim = Dimension::create("R3", 3);

   // For distributed dimensions we need to compute offsets along each
   // dimension (the global offset at each local point with -1 for non-owned
   // points)

   // Distinguish these dimensions by name if this is not the default mesh
   std::string NameSuffix;
   if (MeshName == "Default") {
      NameSuffix = "";
   } else {
      NameSuffix = MeshName;
   }

   // Create the offset and Dimension for NCells
   std::string DimName = "NCells" + NameSuffix;
   HostArray1DI4 CellOffset("CellOffSet", NCellsSize);
   I4 NCellsGlobal = MeshDecomp->NCellsGlobal;
   // NCellsOwned, NCellsSize already set
   for (int Cell = 0; Cell < NCellsSize; ++Cell) {
      if (Cell < NCellsOwned) {
         CellOffset(Cell) = MeshDecomp->CellIDH(Cell) - 1;
      } else {
         CellOffset(Cell) = -1;
      }
   }
   // Create the distributed dimension
   auto NCellsDim =
       Dimension::create(DimName, NCellsGlobal, NCellsSize, CellOffset);

   // Create the offset and Dimension for NEdges
   DimName = "NEdges" + NameSuffix;
   HostArray1DI4 EdgeOffset("EdgeOffSet", NEdgesSize);
   I4 NEdgesGlobal = MeshDecomp->NEdgesGlobal;
   // NEdgesOwned, NEdgesSize already set
   for (int Edge = 0; Edge < NEdgesSize; ++Edge) {
      if (Edge < NEdgesOwned) {
         EdgeOffset(Edge) = MeshDecomp->EdgeIDH(Edge) - 1;
      } else {
         EdgeOffset(Edge) = -1;
      }
   }
   // Create the distributed dimension
   auto NEdgesDim =
       Dimension::create(DimName, NEdgesGlobal, NEdgesSize, EdgeOffset);

   // Create the offset and Dimension for NVertices
   DimName = "NVertices" + NameSuffix;
   HostArray1DI4 VrtxOffset("VertexOffSet", NVerticesSize);
   I4 NVerticesGlobal = MeshDecomp->NVerticesGlobal;
   // NVerticesOwned, NVerticesSize already set
   for (int Vrtx = 0; Vrtx < NVerticesSize; ++Vrtx) {
      if (Vrtx < NVerticesOwned) {
         VrtxOffset(Vrtx) = MeshDecomp->VertexIDH(Vrtx) - 1;
      } else {
         VrtxOffset(Vrtx) = -1;
      }
   }
   // Create the distributed dimension
   auto NVerticesDim =
       Dimension::create(DimName, NVerticesGlobal, NVerticesSize, VrtxOffset);

} // end createDimensions

//------------------------------------------------------------------------------
// Compute the sign of edge contributions to a cell/vertex for each edge
void HorzMesh::computeEdgeSign() {

   OMEGA_SCOPE(o_NEdgesOnCell, NEdgesOnCell);
   OMEGA_SCOPE(o_EdgesOnCell, EdgesOnCell);
   OMEGA_SCOPE(o_CellsOnEdge, CellsOnEdge);
   OMEGA_SCOPE(o_EdgeSignOnCell, EdgeSignOnCell);

   parallelFor(
       {NCellsAll}, KOKKOS_LAMBDA(int Cell) {
          for (int i = 0; i < o_NEdgesOnCell(Cell); i++) {
             int Edge = o_EdgesOnCell(Cell, i);

             // Vector points from cell 0 to cell 1
             if (Cell == o_CellsOnEdge(Edge, 0)) {
                o_EdgeSignOnCell(Cell, i) = -1.0;
             } else {
                o_EdgeSignOnCell(Cell, i) = 1.0;
             }
          }
       });

   EdgeSignOnCellH = createHostMirrorCopy(EdgeSignOnCell);

   EdgeSignOnVertex =
       Array2DReal("EdgeSignOnVertex", NVerticesSize, VertexDegree);

   OMEGA_SCOPE(o_VertexDegree, VertexDegree);
   OMEGA_SCOPE(o_EdgesOnVertex, EdgesOnVertex);
   OMEGA_SCOPE(o_VerticesOnEdge, VerticesOnEdge);
   OMEGA_SCOPE(o_EdgeSignOnVertex, EdgeSignOnVertex);

   parallelFor(
       {NVerticesAll}, KOKKOS_LAMBDA(int Vertex) {
          for (int i = 0; i < o_VertexDegree; i++) {
             int Edge = o_EdgesOnVertex(Vertex, i);

             // Vector points from vertex 0 to vertex 1
             if (Vertex == o_VerticesOnEdge(Edge, 0)) {
                o_EdgeSignOnVertex(Vertex, i) = -1.0;
             } else {
                o_EdgeSignOnVertex(Vertex, i) = 1.0;
             }
          }
       });

   EdgeSignOnVertexH = createHostMirrorCopy(EdgeSignOnVertex);
} // end computeEdgeSign

//------------------------------------------------------------------------------
// Set mesh scaling coefficients for mixing terms in momentum and tracer
// equations so viscosity and diffusion scale with mesh
void HorzMesh::computeMeshScaling() {

   Config *OmegaConfig = Config::getOmegaConfig();
   Config HorzMeshConfig("HorzMesh");
   Error Err = OmegaConfig->get(HorzMeshConfig);
   CHECK_ERROR_ABORT(Err,
                     "HorzMesh: HorzMesh group not found in configuration");

   bool MeshScalingEnable;
   Err = HorzMeshConfig.get("MeshScalingEnable", MeshScalingEnable);
   CHECK_ERROR_ABORT(Err, "HorzMesh: MeshScalingEnable not found in HorzMesh "
                          "configuration");

   // Reference width (median cell width) for the scaling
   computeMedianCellWidth();

   OMEGA_SCOPE(o_MeshScalingDel2, MeshScalingDel2);
   OMEGA_SCOPE(o_MeshScalingDel4, MeshScalingDel4);

   if (MeshScalingEnable) {
      OMEGA_SCOPE(o_MedianCellWidth, MedianCellWidth);
      OMEGA_SCOPE(o_AreaCell, AreaCell);
      OMEGA_SCOPE(o_CellsOnEdge, CellsOnEdge);

      // Mesh scaling is computed using CellWidth derived from AreaCell
      // relative to the global median cell width.
      parallelFor(
          {NEdgesAll}, KOKKOS_LAMBDA(int Edge) {
             const int Cell0 = o_CellsOnEdge(Edge, 0);
             const int Cell1 = o_CellsOnEdge(Edge, 1);
             const Real CellWidth =
                 2.0_Real *
                 Kokkos::sqrt((o_AreaCell(Cell0) + o_AreaCell(Cell1)) /
                              (2.0_Real * Pi));
             const Real Del2Scale = CellWidth / o_MedianCellWidth;

             o_MeshScalingDel2(Edge) = Del2Scale;
             o_MeshScalingDel4(Edge) = Del2Scale * Del2Scale * Del2Scale;
          });

   } else {
      // Without mesh scaling, the mesh scaling is set as 1 and
      // is the same everywhere
      parallelFor(
          {NEdgesAll}, KOKKOS_LAMBDA(int Edge) {
             o_MeshScalingDel2(Edge) = 1.0_Real;
             o_MeshScalingDel4(Edge) = 1.0_Real;
          });
   } // if MeshScalingEnable

   MeshScalingDel2H = createHostMirrorCopy(MeshScalingDel2);
   MeshScalingDel4H = createHostMirrorCopy(MeshScalingDel4);

   // Write the median cell width and the range of the scaling coefficients
   // to the log. Boundary edges have only one valid cell (the other is
   // NCellsAll) and are masked out of the mixing terms, so they are skipped.
   Real MinDel2 = std::numeric_limits<Real>::max();
   Real MaxDel2 = std::numeric_limits<Real>::lowest();
   Real MinDel4 = std::numeric_limits<Real>::max();
   Real MaxDel4 = std::numeric_limits<Real>::lowest();
   for (int Edge = 0; Edge < NEdgesOwned; ++Edge) {
      if (CellsOnEdgeH(Edge, 0) == NCellsAll or
          CellsOnEdgeH(Edge, 1) == NCellsAll)
         continue;
      MinDel2 = std::min(MinDel2, MeshScalingDel2H(Edge));
      MaxDel2 = std::max(MaxDel2, MeshScalingDel2H(Edge));
      MinDel4 = std::min(MinDel4, MeshScalingDel4H(Edge));
      MaxDel4 = std::max(MaxDel4, MeshScalingDel4H(Edge));
   }

   MPI_Comm Comm = Halo::get(MeshName)->getComm();
   std::vector<Real> MinVals =
       globalMinVal(std::vector<Real>{MinDel2, MinDel4}, Comm);
   std::vector<Real> MaxVals =
       globalMaxVal(std::vector<Real>{MaxDel2, MaxDel4}, Comm);

   LOG_INFO("HorzMesh: MeshScalingEnable = {}; median cell width = {} m",
            MeshScalingEnable, MedianCellWidth);
   LOG_INFO("HorzMesh: MeshScalingDel2 min/max = {}/{}", MinVals[0],
            MaxVals[0]);
   LOG_INFO("HorzMesh: MeshScalingDel4 min/max = {}/{}", MinVals[1],
            MaxVals[1]);

} // end computeMeshScaling

//------------------------------------------------------------------------------
// Compute the global median cell width, where the width of each cell is the
// diameter of a circle with the same area. The widths of the owned cells are
// gathered on one task, so the result does not depend on the decomposition.
void HorzMesh::computeMedianCellWidth() {

   MPI_Comm Comm  = Halo::get(MeshName)->getComm();
   const int Root = 0;
   int MyTask;
   int NTasks;
   MPI_Comm_rank(Comm, &MyTask);
   MPI_Comm_size(Comm, &NTasks);

   std::vector<Real> LocWidths(NCellsOwned);
   for (int Cell = 0; Cell < NCellsOwned; ++Cell)
      LocWidths[Cell] = 2.0_Real * std::sqrt(AreaCellH(Cell) / Pi);

   // Gather the number of owned cells on each task to set the offsets
   I4 LocCount = NCellsOwned;
   std::vector<I4> Counts(NTasks, 0);
   std::vector<I4> Offsets(NTasks, 0);
   int Err = MPI_Gather(&LocCount, 1, MPI_INT32_T, Counts.data(), 1,
                        MPI_INT32_T, Root, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("HorzMesh: MPI error gathering cell counts for the median "
                  "cell width");

   std::vector<Real> AllWidths;
   if (MyTask == Root) {
      I4 NCellsTot = 0;
      for (int Task = 0; Task < NTasks; ++Task) {
         Offsets[Task] = NCellsTot;
         NCellsTot += Counts[Task];
      }
      OMEGA_REQUIRE(NCellsTot == NCellsGlobal,
                    "HorzMesh: gathered {} cells for the median cell width, "
                    "expected {}",
                    NCellsTot, NCellsGlobal);
      AllWidths.resize(NCellsTot);
   }

   Err = MPI_Gatherv(LocWidths.data(), LocCount, MPI_RealKind, AllWidths.data(),
                     Counts.data(), Offsets.data(), MPI_RealKind, Root, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("HorzMesh: MPI error gathering cell widths for the median "
                  "cell width");

   // Median of all cell widths (mean of the two middle values when the number
   // of cells is even)
   Real Median = 0.0_Real;
   if (MyTask == Root) {
      auto Mid = AllWidths.begin() + AllWidths.size() / 2;
      std::nth_element(AllWidths.begin(), Mid, AllWidths.end());
      Median = *Mid;
      if (AllWidths.size() % 2 == 0)
         Median =
             0.5_Real * (Median + *std::max_element(AllWidths.begin(), Mid));
   }

   Err = MPI_Bcast(&Median, 1, MPI_RealKind, Root, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("HorzMesh: MPI error broadcasting the median cell width");

   OMEGA_REQUIRE(Median > 0.0_Real,
                 "HorzMesh: median cell width must be positive, got {}",
                 Median);
   MedianCellWidth = Median;

} // end computeMedianCellWidth

//------------------------------------------------------------------------------
// Finish filling arrays read from a file by filling halos and copying to host
void HorzMesh::completeReadArrays() {

   // Get precomputed halo information by name - mesh name is assumed to be
   // same as halo name (including Default)
   Halo *HorzMeshHalo = Halo::get(MeshName);

   // Fill halos for all arrays read from file
   HorzMeshHalo->exchangeFullArrayHalo(XCell, OnCell);
   HorzMeshHalo->exchangeFullArrayHalo(YCell, OnCell);
   HorzMeshHalo->exchangeFullArrayHalo(ZCell, OnCell);
   HorzMeshHalo->exchangeFullArrayHalo(XEdge, OnEdge);
   HorzMeshHalo->exchangeFullArrayHalo(YEdge, OnEdge);
   HorzMeshHalo->exchangeFullArrayHalo(ZEdge, OnEdge);
   HorzMeshHalo->exchangeFullArrayHalo(XVertex, OnVertex);
   HorzMeshHalo->exchangeFullArrayHalo(YVertex, OnVertex);
   HorzMeshHalo->exchangeFullArrayHalo(ZVertex, OnVertex);
   HorzMeshHalo->exchangeFullArrayHalo(LatCell, OnCell);
   HorzMeshHalo->exchangeFullArrayHalo(LonCell, OnCell);
   HorzMeshHalo->exchangeFullArrayHalo(LatEdge, OnEdge);
   HorzMeshHalo->exchangeFullArrayHalo(LonEdge, OnEdge);
   HorzMeshHalo->exchangeFullArrayHalo(LatVertex, OnVertex);
   HorzMeshHalo->exchangeFullArrayHalo(LonVertex, OnVertex);
   HorzMeshHalo->exchangeFullArrayHalo(AreaCell, OnCell);
   HorzMeshHalo->exchangeFullArrayHalo(MeshDensity, OnCell);
   HorzMeshHalo->exchangeFullArrayHalo(AreaTriangle, OnVertex);
   HorzMeshHalo->exchangeFullArrayHalo(KiteAreasOnVertex, OnVertex);
   HorzMeshHalo->exchangeFullArrayHalo(DcEdge, OnEdge);
   HorzMeshHalo->exchangeFullArrayHalo(DvEdge, OnEdge);
   HorzMeshHalo->exchangeFullArrayHalo(AngleEdge, OnEdge);
   HorzMeshHalo->exchangeFullArrayHalo(WeightsOnEdge, OnEdge);
   HorzMeshHalo->exchangeFullArrayHalo(FCell, OnCell);
   HorzMeshHalo->exchangeFullArrayHalo(FEdge, OnEdge);
   HorzMeshHalo->exchangeFullArrayHalo(FVertex, OnVertex);
   // ReconWeightsCell is (NCells, R3, MaxEdges2) to match the mesh
   // file's on-disk dimension order. Halo's generic 3D exchange expects
   // the mesh dimension in the middle (as in Tracers' [Tracer, Cell,
   // Vert]), so exchange each R3 component as a 2D (Cell, MaxEdges2)
   // slice instead of the full 3D array.
   if (OnSphere) {
      for (int IComp = 0; IComp < 3; ++IComp) {
         auto ReconWeightsCellSlice =
             Kokkos::subview(ReconWeightsCell, Kokkos::ALL, IComp, Kokkos::ALL);
         HorzMeshHalo->exchangeFullArrayHalo(ReconWeightsCellSlice, OnCell);
      }
   }

   // Create host copies
   XCellH             = createHostMirrorCopy(XCell);
   YCellH             = createHostMirrorCopy(YCell);
   ZCellH             = createHostMirrorCopy(ZCell);
   XEdgeH             = createHostMirrorCopy(XEdge);
   YEdgeH             = createHostMirrorCopy(YEdge);
   ZEdgeH             = createHostMirrorCopy(ZEdge);
   XVertexH           = createHostMirrorCopy(XVertex);
   YVertexH           = createHostMirrorCopy(YVertex);
   ZVertexH           = createHostMirrorCopy(ZVertex);
   LatCellH           = createHostMirrorCopy(LatCell);
   LonCellH           = createHostMirrorCopy(LonCell);
   LatEdgeH           = createHostMirrorCopy(LatEdge);
   LonEdgeH           = createHostMirrorCopy(LonEdge);
   LatVertexH         = createHostMirrorCopy(LatVertex);
   LonVertexH         = createHostMirrorCopy(LonVertex);
   AreaCellH          = createHostMirrorCopy(AreaCell);
   MeshDensityH       = createHostMirrorCopy(MeshDensity);
   AreaTriangleH      = createHostMirrorCopy(AreaTriangle);
   KiteAreasOnVertexH = createHostMirrorCopy(KiteAreasOnVertex);
   DcEdgeH            = createHostMirrorCopy(DcEdge);
   DvEdgeH            = createHostMirrorCopy(DvEdge);
   AngleEdgeH         = createHostMirrorCopy(AngleEdge);
   WeightsOnEdgeH     = createHostMirrorCopy(WeightsOnEdge);
   FCellH             = createHostMirrorCopy(FCell);
   FEdgeH             = createHostMirrorCopy(FEdge);
   FVertexH           = createHostMirrorCopy(FVertex);
   if (OnSphere) {
      ReconWeightsCellH = createHostMirrorCopy(ReconWeightsCell);
   }

} // end completeReadArrays

//------------------------------------------------------------------------------
// Get default mesh
HorzMesh *HorzMesh::getDefault() { return HorzMesh::DefaultHorzMesh; }

//------------------------------------------------------------------------------
// Get mesh by name
HorzMesh *HorzMesh::get(const std::string Name ///< [in] Name of mesh
) {

   // look for an instance of this name
   auto it = AllHorzMeshes.find(Name);

   // if not found, abort
   if (it == AllHorzMeshes.end())
      ABORT_ERROR("HorzMesh::get: Mesh {} not found. Mesh has not been defined "
                  "or has been removed",
                  Name);

   // found the mesh, return the pointer
   return it->second.get();

} // end get mesh

//------------------------------------------------------------------------------
// Define all Horz mesh fields and associated mesh stream for any mesh I/O
void HorzMesh::defineMeshFields() {

   // First create a field group for mesh fields - simply HorzMesh for default,
   // but add name of mesh for any other instances. We do not add the suffix to
   // the field names because we expect any file containing the mesh fields
   // to use the standard field names.
   // We currently define a mesh group only for fields that are read in after
   // the decomposition (HorzMeshIn). Other mesh groups may be defined later
   // if fields are required in output, but connectivity arrays are not
   // available after the mesh Decomposition has been defined.

   std::string MeshSuffix = "";
   if (MeshName != "Default")
      MeshSuffix = MeshName;

   std::string MeshGroupName = "HorzMeshIn" + MeshSuffix;
   auto MeshGroupIn          = FieldGroup::create(MeshGroupName);

   // Check whether the file name in the input stream matches that of
   // the Decomp. They are typically different only for unit testing where
   // an internal variable overrides the default configuration. The stream
   // name is assumed to be the same as the constructed group name above.

   IOStream::changeFilename(MeshGroupName, MeshFileName);

   // The connectivity arrays are computed in decomp and here contain only local
   // connectivity and local addresses that are dependent on that decomposition.
   // They should not be read or written and are therefore not included in
   // mesh field groups.

   // CellID;                          < global cell ID for each local cell
   // CellsOnCell, CellsOnCellH;       < Indx of cells that neighbor each cell
   // EdgesOnCell, EdgesOnCellH;       < Indx of edges that border each cell
   // NEdgesOnCell, NEdgesOnCellH;     < Num of active edges around each cell
   // VerticesOnCell, VerticesOnCellH; < Indx of vertices bordering each cell
   // CellsOnEdge, CellsOnEdgeH;       < Indx of cells straddling each edge
   // EdgesOnEdge, EdgesOnEdgeH;       < Indx of edges around cells across edge
   // NEdgesOnEdge, NEdgesOnEdgeH;     < Num of edges around cells across edge
   // VerticesOnEdge, VerticesOnEdgeH; < Indx of vertices straddling each edge
   // CellsOnVertex, CellsOnVertexH;   < Indx of cells that share a vertex
   // EdgesOnVertex, EdgesOnVertexH;   < Indx of edges sharing vertex as endpt

   // These fields are all read from the mesh file so are added to both the
   // full mesh group and the input mesh group

   // Coordinate arrays
   // Cell coords
   std::string FieldName = "XCell";
   XCell = Array1DReal("XCell", NCellsSize); // allocate space and init to zero
   int NDims = 1;
   std::vector<std::string> DimNames(NDims);
   DimNames[0] = "NCells";
   auto XCellField =
       Field::create(FieldName,                           // field name
                     "X Coordinates of cell centers (m)", // long Name
                     "m",                                 // units
                     "",                                  // CF standard Name
                     -7.0E+6,                             // min valid value
                     7.0E+6,                              // max valid value
                     NDims,                               // num of dimensions
                     DimNames,                            // dimension names
                     false                                // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, XCell);

   FieldName = "YCell";
   YCell = Array1DReal("YCell", NCellsSize); // allocate space and init to zero
   auto YCellField =
       Field::create(FieldName,                           // field name
                     "Y Coordinates of cell centers (m)", // long Name
                     "m",                                 // units
                     "",                                  // CF standard Name
                     -7.0E+6,                             // min valid value
                     7.0E+6,                              // max valid value
                     NDims,                               // num of dimensions
                     DimNames,                            // dimension names
                     false                                // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, YCell);

   FieldName = "ZCell";
   ZCell = Array1DReal("ZCell", NCellsSize); // allocate space and init to zero
   auto ZCellField =
       Field::create(FieldName,                           // field name
                     "Z Coordinates of cell centers (m)", // long Name
                     "m",                                 // units
                     "",                                  // CF standard Name
                     -7.0E+6,                             // min valid value
                     7.0E+6,                              // max valid value
                     NDims,                               // num of dimensions
                     DimNames,                            // dimension names
                     false                                // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, ZCell);

   FieldName = "LatCell";
   LatCell = Array1DReal("LatCell", NCellsSize); // allocate space init to zero
   auto LatCellField =
       Field::create(FieldName,                              // field name
                     "Latitude coordinates of cell centers", // long Name
                     "radians",                              // units
                     "latitude",                             // CF standard Name
                     -3.1415927,                             // min valid value
                     3.1415927,                              // max valid value
                     NDims,    // num of dimensions
                     DimNames, // dimension names
                     false     // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, LatCell);

   FieldName = "LonCell";
   LonCell = Array1DReal("LonCell", NCellsSize); // allocate space init to zero
   auto LonCellField =
       Field::create(FieldName,                               // field name
                     "Longitude coordinates of cell centers", // long Name
                     "radians",                               // units
                     "longitude", // CF standard Name
                     0.0,         // min valid value
                     6.28319,     // max valid value
                     NDims,       // num of dimensions
                     DimNames,    // dimension names
                     false        // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, LonCell);

   // Edge coords
   DimNames[0] = "NEdges";
   FieldName   = "XEdge";
   XEdge = Array1DReal("XEdge", NEdgesSize); // allocate space init to zero
   auto XEdgeField =
       Field::create(FieldName,                         // field name
                     "X Coordinates of cell edges (m)", // long Name
                     "m",                               // units
                     "",                                // CF standard Name
                     -7.0E+6,                           // min valid value
                     7.0E+6,                            // max valid value
                     NDims,                             // num of dimensions
                     DimNames,                          // dimension names
                     false                              // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, XEdge);

   FieldName = "YEdge";
   YEdge     = Array1DReal("YEdge", NEdgesSize); // allocate space init to zero
   auto YEdgeField =
       Field::create(FieldName,                         // field name
                     "Y Coordinates of cell edges (m)", // long Name
                     "m",                               // units
                     "",                                // CF standard Name
                     -7.0E+6,                           // min valid value
                     7.0E+6,                            // max valid value
                     NDims,                             // num of dimensions
                     DimNames,                          // dimension names
                     false                              // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, YEdge);

   FieldName = "ZEdge";
   ZEdge     = Array1DReal("ZEdge", NEdgesSize); // allocate space init to zero
   auto ZEdgeField =
       Field::create(FieldName,                         // field name
                     "Z Coordinates of cell edges (m)", // long Name
                     "m",                               // units
                     "",                                // CF standard Name
                     -7.0E+6,                           // min valid value
                     7.0E+6,                            // max valid value
                     NDims,                             // num of dimensions
                     DimNames,                          // dimension names
                     false                              // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, ZEdge);

   FieldName = "LatEdge";
   LatEdge = Array1DReal("LatEdge", NEdgesSize); // allocate space init to zero
   auto LatEdgeField =
       Field::create(FieldName,                            // field name
                     "Latitude coordinates of cell edges", // long Name
                     "radians",                            // units
                     "latitude",                           // CF standard Name
                     -3.1415927,                           // min valid value
                     3.1415927,                            // max valid value
                     NDims,                                // num of dimensions
                     DimNames,                             // dimension names
                     false                                 // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, LatEdge);

   FieldName = "LonEdge";
   LonEdge = Array1DReal("LonEdge", NEdgesSize); // allocate space init to zero
   auto LonEdgeField =
       Field::create(FieldName,                             // field name
                     "Longitude coordinates of cell edges", // long Name
                     "radians",                             // units
                     "longitude",                           // CF standard Name
                     0.0,                                   // min valid value
                     6.28319,                               // max valid value
                     NDims,                                 // num of dimensions
                     DimNames,                              // dimension names
                     false // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, LonEdge);

   // Vertex coordinates
   DimNames[0] = "NVertices";
   FieldName   = "XVertex";
   XVertex     = Array1DReal("XVertex", NVerticesSize); // allocate space
   auto XVertexField =
       Field::create(FieldName,                            // field name
                     "X Coordinates of cell vertices (m)", // long Name
                     "m",                                  // units
                     "",                                   // CF standard Name
                     -7.0E+6,                              // min valid value
                     7.0E+6,                               // max valid value
                     NDims,                                // num of dimensions
                     DimNames,                             // dimension names
                     false                                 // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, XVertex);

   FieldName = "YVertex";
   YVertex   = Array1DReal("YVertex", NVerticesSize); // allocate space
   auto YVertexField =
       Field::create(FieldName,                            // field name
                     "Y Coordinates of cell vertices (m)", // long Name
                     "m",                                  // units
                     "",                                   // CF standard Name
                     -7.0E+6,                              // min valid value
                     7.0E+6,                               // max valid value
                     NDims,                                // num of dimensions
                     DimNames,                             // dimension names
                     false                                 // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, YVertex);

   FieldName = "ZVertex";
   ZVertex   = Array1DReal("ZVertex", NVerticesSize); // allocate space
   auto ZVertexField =
       Field::create(FieldName,                            // field name
                     "Z Coordinates of cell vertices (m)", // long Name
                     "m",                                  // units
                     "",                                   // CF standard Name
                     -7.0E+6,                              // min valid value
                     7.0E+6,                               // max valid value
                     NDims,                                // num of dimensions
                     DimNames,                             // dimension names
                     false                                 // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, ZVertex);

   FieldName = "LatVertex";
   LatVertex = Array1DReal("LatVertex", NVerticesSize); // allocate space
   auto LatVertexField =
       Field::create(FieldName,                               // field name
                     "Latitude coordinates of cell vertices", // long Name
                     "radians",                               // units
                     "latitude", // CF standard Name
                     -3.1415927, // min valid value
                     3.1415927,  // max valid value
                     NDims,      // num of dimensions
                     DimNames,   // dimension names
                     false       // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, LatVertex);

   FieldName = "LonVertex";
   LonVertex = Array1DReal("LonVertex", NVerticesSize); // allocate space
   auto LonVertexField =
       Field::create(FieldName,                                // field name
                     "Longitude coordinates of cell vertices", // long Name
                     "radians",                                // units
                     "longitude", // CF standard Name
                     0.0,         // min valid value
                     6.28319,     // max valid value
                     NDims,       // num of dimensions
                     DimNames,    // dimension names
                     false        // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, LonVertex);

   // Other mesh properties
   // Mesh areas, lengths, and angles
   DimNames[0]        = "NCells";
   FieldName          = "AreaCell";
   AreaCell           = Array1DReal("AreaCell", NCellsSize); // allocate space
   auto AreaCellField = Field::create(FieldName,             // field name
                                      "Area of each cell (m2)", // long name
                                      "m2",                     // units
                                      "cell_area", // CF standard name
                                      0.0,         // min valid value
                                      9.99E+30,    // max valid value
                                      NDims,       // num of dimensions
                                      DimNames,    // dimension names
                                      false        // not time dependent
   );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, AreaCell);

   FieldName   = "MeshDensity";
   MeshDensity = Array1DReal("MeshDensity", NCellsSize); // allocate space
   auto MeshDensityField =
       Field::create(FieldName, // field name
                     "Value of density function used to generate mesh cell at"
                     " cell centers", // long name
                     "",              // units
                     "",              // CF standard name
                     0.0,             // min valid value
                     9.99E+30,        // max valid value
                     NDims,           // num of dimensions
                     DimNames,        // dimension names
                     false            // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, MeshDensity);

   DimNames[0]  = "NVertices";
   FieldName    = "AreaTriangle";
   AreaTriangle = Array1DReal("AreaTriangle", NVerticesSize); // allocate space
   auto AreaTriangleField =
       Field::create(FieldName, // field name
                     "Area of each triangle in the dual grid (m2)", // lng name
                     "m2",                                          // units
                     "cell_area", // CF standard name
                     0.0,         // min valid value
                     9.99E+30,    // max valid value
                     NDims,       // num of dimensions
                     DimNames,    // dimension names
                     false        // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, AreaTriangle);

   DimNames[0] = "NEdges";
   FieldName   = "DvEdge";
   DvEdge      = Array1DReal("DvEdge", NEdgesSize); // allocate space
   auto DvEdgeField =
       Field::create(FieldName, // field name
                     "Length of each edge, computed as the distance between"
                     " verticesOnEdge (m)", // long name
                     "m",                   // units
                     "",                    // CF standard name
                     0.0,                   // min valid value
                     9.99E+30,              // max valid value
                     NDims,                 // num of dimensions
                     DimNames,              // dimension names
                     false                  // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, DvEdge);

   FieldName = "DcEdge";
   DcEdge    = Array1DReal("DcEdge", NEdgesSize); // allocate space
   auto DcEdgeField =
       Field::create(FieldName, // field name
                     "Length of each edge, computed as the distance between"
                     " CellsOnEdge (m)", // long name
                     "m",                // units
                     "",                 // CF standard name
                     0.0,                // min valid value
                     9.99E+30,           // max valid value
                     NDims,              // num of dimensions
                     DimNames,           // dimension names
                     false               // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, DcEdge);

   FieldName = "AngleEdge";
   AngleEdge = Array1DReal("AngleEdge", NEdgesSize); // allocate space
   auto AngleEdgeField =
       Field::create(FieldName, // field name
                     "Angle the edge normal makes with local eastward direction"
                     " (radians)", // long name
                     "radians",    // units
                     "",           // CF standard name
                     -3.1415927,   // min valid value
                     3.1415927,    // max valid value
                     NDims,        // num of dimensions
                     DimNames,     // dimension names
                     false         // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, AngleEdge);

   NDims = 2;
   DimNames.resize(2);
   DimNames[0] = "NVertices";
   DimNames[1] = "VertexDegree";
   FieldName   = "KiteAreasOnVertex";
   KiteAreasOnVertex =
       Array2DReal("KiteAreasOnVertex", NVerticesSize, VertexDegree);
   auto KiteAreasOnVertexField =
       Field::create(FieldName, // field name
                     "Area of the portions of each dual cell that are part of "
                     "each cellsOnVertex (m2)", // long name
                     "m2",                      // units
                     "",                        // CF standard name
                     0.0,                       // min valid value
                     9.99E+30,                  // max valid value
                     NDims,                     // num of dimensions
                     DimNames,                  // dimension names
                     false                      // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array2DReal>(FieldName, KiteAreasOnVertex);

   // Mesh weights
   DimNames[0]   = "NEdges";
   DimNames[1]   = "MaxEdges2";
   FieldName     = "WeightsOnEdge";
   WeightsOnEdge = Array2DReal("KiteAreasOnVertex", NEdgesSize, MaxEdges2);
   auto WeightsOnEdgeField =
       Field::create(FieldName, // field name
                     "Reconstruction weights associated with each of the"
                     " edgesOnEdge", // long name
                     "",             // units
                     "",             // CF standard name
                     -1.0,           // min valid value
                     1.0,            // max valid value
                     NDims,          // num of dimensions
                     DimNames,       // dimension names
                     false           // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array2DReal>(FieldName, WeightsOnEdge);

   // Coriolis parameter at the cells, edges, and vertices
   NDims = 1;
   DimNames.resize(1);
   DimNames[0] = "NEdges";
   FieldName   = "FEdge";
   FEdge       = Array1DReal("FEdge", NEdgesSize);
   auto FEdgeField =
       Field::create(FieldName,                                   // field name
                     "Coriolis parameter at edges (radians s-1)", // long name
                     "radians s-1",                               // units
                     "coriolis_parameter", // CF standard name
                     -1.0E-3,              // min valid value
                     1.0E-3,               // max valid value
                     NDims,                // num of dimensions
                     DimNames,             // dimension names
                     false                 // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, FEdge);

   DimNames[0] = "NVertices";
   FieldName   = "FVertex";
   FVertex     = Array1DReal("FVertex", NVerticesSize);
   auto FVertexField =
       Field::create(FieldName, // field name
                     "Coriolis parameter at vertices (radians s-1)",
                     "radians s-1",        // units
                     "coriolis_parameter", // CF standard name
                     -1.0E-3,              // min valid value
                     1.0E-3,               // max valid value
                     NDims,                // num of dimensions
                     DimNames,             // dimension names
                     false                 // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, FVertex);

   DimNames[0] = "NCells";
   FieldName   = "FCell";
   FCell       = Array1DReal("FCell", NCellsSize);
   auto FCellField =
       Field::create(FieldName, // field name
                     "Coriolis parameter at cell centers (radians s-1)",
                     "radians s-1",        // units
                     "coriolis_parameter", // CF standard name
                     -1.0E-3,              // min valid value
                     1.0E-3,               // max valid value
                     NDims,                // num of dimensions
                     DimNames,             // dimension names
                     false                 // not time dependent
       );
   MeshGroupIn->addField(FieldName);
   Field::attachFieldData<Array1DReal>(FieldName, FCell);

   if (OnSphere) {
      // Vector Reconstruction
      // NEdgesReconOnCell/ReconStencilCell come from Decomp (see
      // constructor), not read here. ReconWeightsCell is pure data and
      // is read as a normal mesh field below.
      NDims = 3;
      DimNames.resize(NDims);
      DimNames[0] = "NCells";
      DimNames[1] = "R3";
      DimNames[2] = "MaxEdges2";
      FieldName   = "ReconWeightsCell";
      ReconWeightsCell =
          Array3DReal("ReconWeightsCell", NCellsSize, 3, MaxEdges2);
      auto ReconWeightsCellField = Field::create(
          FieldName, // field name
          "weights to reconstruct Cartesian vector field from edge-"
          "normal values at cell centers",     // long name
          "",                                  // units
          "",                                  // CF standard name
          std::numeric_limits<Real>::lowest(), // min valid value
          std::numeric_limits<Real>::max(),    // max valid value
          NDims,                               // num of dimensions
          DimNames,                            // dimension names
          false                                // not time dependent
      );
      MeshGroupIn->addField(FieldName);
      Field::attachFieldData<Array3DReal>(FieldName, ReconWeightsCell);
   }

   // The sign and scaling fields are computed internally so are allocated
   // here. Because they are not read from a file or typically written to a
   // file, we do not create field metadata, though the definitions are
   // included here in comments in case they are needed later.

   EdgeSignOnCell = Array2DReal("EdgeSignOnCell", NCellsSize, MaxEdges);
   EdgeSignOnVertex =
       Array2DReal("EdgeSignOnVertex", NVerticesSize, VertexDegree);
   MeshScalingDel2 = Array1DReal("MeshScalingDel2", NEdgesSize);
   MeshScalingDel4 = Array1DReal("MeshScalingDel4", NEdgesSize);

   // Sign vectors
   //
   // NDims = 2;
   // DimNames.resize(2);
   // DimNames[0] = "NCells";
   // DimNames[1] = "MaxEdges";
   // FieldName = "EdgeSignOnCell";
   // auto EdgeSignOnCellField =
   //    Field::create(FieldName,            // field name
   //                  "Sign of vector connecting cells", // long name
   //                  "",                   // units
   //                  "",                   // CF standard name
   //                  -1.0,                 // min valid value
   //                  1.0,                  // max valid value
   //                  -9.99E+30,            // scalar for undefined entries
   //                  NDims,                // num of dimensions
   //                  DimNames,             // dimension names
   //                  false                 // not time dependent
   //    );
   // MeshGroupIn->addField(FieldName);
   // Field::attachFieldData<Array2DReal>(FieldName, EdgeSignOnCell);
   //
   // DimNames[0] = "NVertices";
   // DimNames[1] = "VertexDegree";
   // FieldName = "EdgeSignOnVertex";
   // auto EdgeSignOnVertexField =
   //    Field::create(FieldName,            // field name
   //                  "Sign of vector connecting vertices", // long name
   //                  "",                   // units
   //                  "",                   // CF standard name
   //                  -1.0,                 // min valid value
   //                  1.0,                  // max valid value
   //                  -9.99E+30,            // scalar for undefined entries
   //                  NDims,                // num of dimensions
   //                  DimNames,             // dimension names
   //                  false                 // not time dependent
   //    );
   // MeshGroupIn->addField(FieldName);
   // Field::attachFieldData<Array2DReal>(FieldName, EdgeSignOnVertex);
   //
   // Mesh scaling
   //
   // NDims = 1;
   // DimNames.resize(1);
   // DimNames[0] = "NEdges";
   // FieldName = "MeshScalingDel2";
   // auto MeshScalingDel2Field =
   //    Field::create(FieldName,            // field name
   //                  "Coefficient to Laplacian mixing terms", // long name
   //                  "",                   // units
   //                  "",                   // CF standard name
   //                  -9.99E+30,            // min valid value
   //                  9.99E+30,             // max valid value
   //                  -9.99E+30,            // scalar for undefined entries
   //                  NDims,                // num of dimensions
   //                  DimNames,             // dimension names
   //                  false                 // not time dependent
   //    );
   // MeshGroupIn->addField(FieldName);
   // Field::attachFieldData<Array1DReal>(FieldName, MeshScalingDel2);
   //
   // FieldName = "MeshScalingDel4";
   // auto MeshScalingDel4Field =
   //    Field::create(FieldName,            // field name
   //                  "Coefficient to biharmonic mixing terms", // long name
   //                  "",                   // units
   //                  "",                   // CF standard name
   //                  -9.99E+30,            // min valid value
   //                  9.99E+30,             // max valid value
   //                  -9.99E+30,            // scalar for undefined entries
   //                  NDims,                // num of dimensions
   //                  DimNames,             // dimension names
   //                  false                 // not time dependent
   //    );
   // MeshGroupIn->addField(FieldName);
   // Field::attachFieldData<Array1DReal>(FieldName, MeshScalingDel4);

} // end defineMeshFields

//------------------------------------------------------------------------------

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
