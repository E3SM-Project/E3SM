//===-- base/Decomp.cpp - domain decomposition methods ----------*- C++ -*-===//
//
// The decomposition (Decomp) class partitions an OMEGA horizontal domain into
// a number of sub-domains that can be distributed across a parallel
// environment. Currently, this relies on the (Par)Metis partitioning package
// to partition a set of nodes (cells) based on the adjacency graph created
// during the OMEGA mesh generators. A default decomposition is created
// from the default Machine Env that specifies the MPI details and number
// of MPI tasks. Other decompositions can be created with the same mesh
// on any of the subset environments that are possible in MachEnv.
// The Decomp class stores a number of index-space arrays that describe
// the partition, neighbor information, global IDs for cell, edge and
// vertex points in an Omega mesh.
//
//===----------------------------------------------------------------------===//

#include "Decomp.h"
#include "Broadcast.h"
#include "Config.h"
#include "DataTypes.h"
#include "Error.h"
#include "IO.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "mpi.h"
#include "parmetis.h"

#include <algorithm>
#include <set>
#include <string>
#include <vector>

namespace OMEGA {

// create the static class members
Decomp *Decomp::DefaultDecomp = nullptr;
std::map<std::string, std::unique_ptr<Decomp>> Decomp::AllDecomps;

// Some useful local utility routines
//------------------------------------------------------------------------------
// checks a global cell/edge/vertex ID value to see if it is within range.
// The readMesh function must have been called to define global sizes before
// calling any of these.

bool Decomp::validCellID(I4 InCellID) {
   return (InCellID > 0 && InCellID <= NCellsGlobal);
}
bool Decomp::validEdgeID(I4 InEdgeID) {
   return (InEdgeID > 0 && InEdgeID <= NEdgesGlobal);
}
bool Decomp::validVertexID(I4 InVertexID) {
   return (InVertexID > 0 && InVertexID <= NVerticesGlobal);
}

//------------------------------------------------------------------------------
// Local routine that searches a std::vector<I4> for a particular entry and
// returns the index of that entry. If not found, the size is returned
// (corresponding to the last index + 1)

I4 srchVector(const std::vector<I4> &InVector, // vector to search
              I4 Value                         // value to search for
) {

   // first use the std::find routine to determine the iterator location
   auto it = std::find(InVector.begin(), InVector.end(), Value);
   // now translate the iterator into an actual vector index
   I4 LocIndx = std::distance(InVector.begin(), it);

   return LocIndx;

} // end function srchVector (std::vector)

//------------------------------------------------------------------------------
// A search routine for vectors in which the vector is a Kokkos array rather
// than a std::vector. It searches for a value and returns the first index of
// that // entry. If not found, the size is returned (corresponding to the
// last index + 1)

I4 srchVector(HostArray1DI4 InVector, // vector to search
              I4 Value                // value to search for
) {

   // extract the vector length
   // I4 VecSize  = InVector.totElems();
   I4 VecSize  = InVector.size();
   I4 LocIndex = VecSize; // set default to size (return value if not found)

   // Loop over elements, searching for Value
   for (int n = 0; n < VecSize; ++n) {
      if (InVector(n) == Value) { // found value
         LocIndex = n;
         break;
      }
   }

   return LocIndex;

} // end function srchVector (Kokkos)

// Routines needed for creating the decomposition
//------------------------------------------------------------------------------
// Reads mesh adjacency, index and size information from a file. This includes
// all neighbor information and connections between cells, edges, vertices.
// These are initially read into a uniform linear distribution across MPI tasks
// and redistributed later after the decomposition is complete.

void readMesh(const int MeshFileID, // file ID for open mesh file
              const MachEnv *InEnv, // input machine environment for MPI layout
              I4 &NCellsGlobal,     // total number of cells
              I4 &NEdgesGlobal,     // total number of edges
              I4 &NVerticesGlobal,  // total number of vertices
              I4 &MaxEdges,         // max number of edges on a cell
              I4 &MaxCellsOnEdge,   // max number of cells sharing edge
              I4 &VertexDegree,     // number of cells/edges sharing vrtx
              std::vector<I4> &CellsOnCellInit, // cell neighbors for each cell
              std::vector<I4> &EdgesOnCellInit, // edge IDs for each cell edge
              std::vector<I4> &VerticesOnCellInit, // vertices around each cell
              std::vector<I4> &CellsOnEdgeInit,    // cell IDs sharing edge
              std::vector<I4> &EdgesOnEdgeInit,    // all edges neighboring edge
              std::vector<I4> &VerticesOnEdgeInit, // vertices on ends of edge
              std::vector<I4> &CellsOnVertexInit,  // cells meeting at each vrtx
              std::vector<I4> &EdgesOnVertexInit   // edges meeting at each vrtx
) {

   Error Err; // error code for IO calls

   // Retrieve some info on the MPI layout
   MPI_Comm Comm = InEnv->getComm();
   I4 NumTasks   = InEnv->getNumTasks();
   I4 MyTask     = InEnv->getMyTask();

   // Read in mesh size information - these are dimension lengths in
   // the input mesh file. Check both the name under Omega name conventions
   // and the older MPAS name.
   std::string DimName    = "NCells";
   std::string DimNameOld = "nCells";
   I4 NCellsID;
   Err = IO::getDimFromFile(MeshFileID, DimName, NCellsID, NCellsGlobal);
   if (Err.isFail()) { // dim not found, try again with older MPAS name
      Err = IO::getDimFromFile(MeshFileID, DimNameOld, NCellsID, NCellsGlobal);
      CHECK_ERROR_ABORT(Err, "Decomp: error reading nCells");
      if (NCellsGlobal <= 0)
         ABORT_ERROR("Decomp: Bad NCells ({}) read from file", NCellsGlobal);
   }

   DimName    = "NEdges";
   DimNameOld = "nEdges";
   I4 NEdgesID;
   Err = IO::getDimFromFile(MeshFileID, DimName, NEdgesID, NEdgesGlobal);
   if (Err.isFail()) { // dim not found, try again with older MPAS name
      Err = IO::getDimFromFile(MeshFileID, DimNameOld, NEdgesID, NEdgesGlobal);
      CHECK_ERROR_ABORT(Err, "Decomp: error reading NEdges");
      if (NEdgesGlobal <= 0)
         ABORT_ERROR("Decomp: Bad NEdges ({}) read from file", NEdgesGlobal);
   }

   DimName    = "NVertices";
   DimNameOld = "nVertices";
   I4 NVerticesID;
   Err = IO::getDimFromFile(MeshFileID, DimName, NVerticesID, NVerticesGlobal);
   if (Err.isFail()) { // dim not found, try again with older MPAS name
      Err = IO::getDimFromFile(MeshFileID, DimNameOld, NVerticesID,
                               NVerticesGlobal);
      CHECK_ERROR_ABORT(Err, "Decomp: error reading NVertices");
      if (NVerticesGlobal <= 0)
         ABORT_ERROR("Decomp: Bad NVertices ({}) read from file",
                     NVerticesGlobal);
   }

   DimName    = "MaxEdges";
   DimNameOld = "maxEdges";
   I4 MaxEdgesID;
   Err = IO::getDimFromFile(MeshFileID, DimName, MaxEdgesID, MaxEdges);
   if (Err.isFail()) { // dim not found, try again with older MPAS name
      Err = IO::getDimFromFile(MeshFileID, DimNameOld, MaxEdgesID, MaxEdges);
      CHECK_ERROR_ABORT(Err, "Decomp: error reading MaxEdges");
      if (MaxEdges <= 0)
         ABORT_ERROR("Decomp: Bad MaxEdges ({}) read from file", MaxEdges);
   }

   DimName    = "VertexDegree";
   DimNameOld = "vertexDegree";
   I4 VertexDegreeID;
   Err = IO::getDimFromFile(MeshFileID, DimName, VertexDegreeID, VertexDegree);
   if (Err.isFail()) { // dim not found, try again with older MPAS name
      Err = IO::getDimFromFile(MeshFileID, DimNameOld, VertexDegreeID,
                               VertexDegree);
      CHECK_ERROR_ABORT(Err, "Decomp: error reading VertexDegree");
      if (VertexDegree <= 0)
         ABORT_ERROR("Decomp: Bad VertexDegree ({}) read from file",
                     VertexDegree);
   }

   DimName    = "MaxCellsOnEdge";
   DimNameOld = "TWO";
   I4 MaxCellsOnEdgeID;
   Err = IO::getDimFromFile(MeshFileID, DimName, MaxCellsOnEdgeID,
                            MaxCellsOnEdge);
   if (Err.isFail()) { // dim not found, try again with older MPAS name
      Err = IO::getDimFromFile(MeshFileID, DimNameOld, MaxCellsOnEdgeID,
                               MaxCellsOnEdge);
      CHECK_ERROR_ABORT(Err, "Decomp: error reading MaxCellsOnEdge");
      if (MaxCellsOnEdge <= 0)
         ABORT_ERROR("Decomp: Bad MaxCellsOnEdge/TWO ({}) read from file",
                     MaxCellsOnEdge);
   }
   I4 MaxEdgesOnEdge = 2 * MaxEdges; // 2*MaxCellsOnEdge

   // Create the linear decompositions for parallel IO
   // Determine the size of each block, divided as evenly as possible
   I4 NCellsChunk    = (NCellsGlobal - 1) / NumTasks + 1;
   I4 NEdgesChunk    = (NEdgesGlobal - 1) / NumTasks + 1;
   I4 NVerticesChunk = (NVerticesGlobal - 1) / NumTasks + 1;

   // If global size did not divide evenly over processors, the last processor
   // only has the remaining indices and not a full block (chunk)
   I4 NCellsLocal;
   I4 NEdgesLocal;
   I4 NVerticesLocal;
   if (MyTask < NumTasks - 1) {
      NCellsLocal    = NCellsChunk;
      NEdgesLocal    = NEdgesChunk;
      NVerticesLocal = NVerticesChunk;
   } else { // last MPI task
      I4 StartAdd    = NCellsChunk * (NumTasks - 1);
      NCellsLocal    = NCellsGlobal - StartAdd;
      StartAdd       = NEdgesChunk * (NumTasks - 1);
      NEdgesLocal    = NEdgesGlobal - StartAdd;
      StartAdd       = NVerticesChunk * (NumTasks - 1);
      NVerticesLocal = NVerticesGlobal - StartAdd;
   }

   // Define offset index for each array (essentially the global index)
   // for the parallel IO decomposition. The XxOnCell arrays are all the
   // same size and decomposition

   I4 OnCellSize = NCellsChunk * MaxEdges;
   I4 NDims      = 2;
   std::vector<I4> OnCellDims{NCellsGlobal, MaxEdges};
   std::vector<I4> OnCellOffset(OnCellSize, -1);
   for (int Cell = 0; Cell < NCellsLocal; ++Cell) {
      I4 CellGlob = MyTask * NCellsChunk + Cell;
      for (int Edge = 0; Edge < MaxEdges; ++Edge) {
         OnCellOffset[Cell * MaxEdges + Edge] = CellGlob * MaxEdges + Edge;
      }
   }

   // Now we do the same for XXOnEdge and XXOnVertex arrays
   // The EdgesOnEdge is a different size so needs its own
   I4 OnEdgeSize   = NEdgesChunk * MaxCellsOnEdge;
   I4 OnEdgeSize2  = NEdgesChunk * MaxEdgesOnEdge;
   I4 OnVertexSize = NVerticesChunk * VertexDegree;
   std::vector<I4> OnEdgeDims{NEdgesGlobal, MaxCellsOnEdge};
   std::vector<I4> OnEdgeDims2{NEdgesGlobal, MaxEdgesOnEdge};
   std::vector<I4> OnVertexDims{NVerticesGlobal, VertexDegree};
   std::vector<I4> OnEdgeOffset(OnEdgeSize, -1);
   std::vector<I4> OnEdgeOffset2(OnEdgeSize2, -1);
   std::vector<I4> OnVertexOffset(OnVertexSize, -1);
   for (int Edge = 0; Edge < NEdgesLocal; ++Edge) {
      I4 EdgeGlob = MyTask * NEdgesChunk + Edge;
      for (int Cell = 0; Cell < MaxCellsOnEdge; ++Cell) {
         OnEdgeOffset[Edge * MaxCellsOnEdge + Cell] =
             EdgeGlob * MaxCellsOnEdge + Cell;
      }
      for (int Cell = 0; Cell < MaxEdgesOnEdge; ++Cell) {
         OnEdgeOffset2[Edge * MaxEdgesOnEdge + Cell] =
             EdgeGlob * MaxEdgesOnEdge + Cell;
      }
   }
   for (int Vrtx = 0; Vrtx < NVerticesLocal; ++Vrtx) {
      I4 VertexGlob = MyTask * NVerticesChunk + Vrtx;
      for (int Cell = 0; Cell < VertexDegree; ++Cell) {
         OnVertexOffset[Vrtx * VertexDegree + Cell] =
             VertexGlob * VertexDegree + Cell;
      } // end loop VertexDegree
   } // end loop NVerticesLocal

   // Create the parallel IO decompositions
   IO::Rearranger Rearr = IO::RearrBox;
   I4 OnCellDecomp      = IO::createDecomp(IO::IOTypeI4, NDims, OnCellDims,
                                           OnCellSize, OnCellOffset, Rearr);
   I4 OnEdgeDecomp      = IO::createDecomp(IO::IOTypeI4, NDims, OnEdgeDims,
                                           OnEdgeSize, OnEdgeOffset, Rearr);
   I4 OnEdgeDecomp2     = IO::createDecomp(IO::IOTypeI4, NDims, OnEdgeDims2,
                                           OnEdgeSize2, OnEdgeOffset2, Rearr);
   I4 OnVertexDecomp    = IO::createDecomp(IO::IOTypeI4, NDims, OnVertexDims,
                                           OnVertexSize, OnVertexOffset, Rearr);

   // Now read the connectivity arrays. Try reading under the new Omega
   // name convention and the older MPAS mesh names.
   CellsOnCellInit.resize(OnCellSize);
   EdgesOnCellInit.resize(OnCellSize);
   VerticesOnCellInit.resize(OnCellSize);
   CellsOnEdgeInit.resize(OnEdgeSize);
   EdgesOnEdgeInit.resize(OnEdgeSize2);
   VerticesOnEdgeInit.resize(OnEdgeSize);
   CellsOnVertexInit.resize(OnVertexSize);
   EdgesOnVertexInit.resize(OnVertexSize);

   std::string VarName    = "CellsOnCell";
   std::string VarNameOld = "cellsOnCell";
   int CellsOnCellID;
   Err = IO::readArray(&CellsOnCellInit[0], OnCellSize, VarName, MeshFileID,
                       OnCellDecomp, CellsOnCellID);
   if (Err.isFail()) { // not found, try again under older name
      Err = IO::readArray(&CellsOnCellInit[0], OnCellSize, VarNameOld,
                          MeshFileID, OnCellDecomp, CellsOnCellID);
      CHECK_ERROR_ABORT(Err, "Decomp: error reading CellsOnCell");
   }

   VarName    = "EdgesOnCell";
   VarNameOld = "edgesOnCell";
   int EdgesOnCellID;
   Err = IO::readArray(&EdgesOnCellInit[0], OnCellSize, VarName, MeshFileID,
                       OnCellDecomp, EdgesOnCellID);
   if (Err.isFail()) { // not found, try again under older name
      Err = IO::readArray(&EdgesOnCellInit[0], OnCellSize, VarNameOld,
                          MeshFileID, OnCellDecomp, EdgesOnCellID);
      CHECK_ERROR_ABORT(Err, "Decomp: error reading EdgesOnCell");
   }

   VarName    = "VerticesOnCell";
   VarNameOld = "verticesOnCell";
   int VerticesOnCellID;
   Err = IO::readArray(&VerticesOnCellInit[0], OnCellSize, VarName, MeshFileID,
                       OnCellDecomp, VerticesOnCellID);
   if (Err.isFail()) { // not found, try again under older name
      Err = IO::readArray(&VerticesOnCellInit[0], OnCellSize, VarNameOld,
                          MeshFileID, OnCellDecomp, VerticesOnCellID);
      CHECK_ERROR_ABORT(Err, "Decomp: error reading EdgesOnCell");
   }

   VarName    = "CellsOnEdge";
   VarNameOld = "cellsOnEdge";
   int CellsOnEdgeID;
   Err = IO::readArray(&CellsOnEdgeInit[0], OnEdgeSize, VarName, MeshFileID,
                       OnEdgeDecomp, CellsOnEdgeID);
   if (Err.isFail()) { // not found, try again under older name
      Err = IO::readArray(&CellsOnEdgeInit[0], OnEdgeSize, VarNameOld,
                          MeshFileID, OnEdgeDecomp, CellsOnEdgeID);
      CHECK_ERROR_ABORT(Err, "Decomp: error reading CellsOnEdge");
   }

   VarName    = "EdgesOnEdge";
   VarNameOld = "edgesOnEdge";
   int EdgesOnEdgeID;
   Err = IO::readArray(&EdgesOnEdgeInit[0], OnEdgeSize2, VarName, MeshFileID,
                       OnEdgeDecomp2, EdgesOnEdgeID);
   if (Err.isFail()) { // not found, try again under older name
      Err = IO::readArray(&EdgesOnEdgeInit[0], OnEdgeSize2, VarNameOld,
                          MeshFileID, OnEdgeDecomp2, EdgesOnEdgeID);
      CHECK_ERROR_ABORT(Err, "Decomp: error reading EdgesOnEdge");
   }

   VarName    = "VerticesOnEdge";
   VarNameOld = "verticesOnEdge";
   int VerticesOnEdgeID;
   Err = IO::readArray(&VerticesOnEdgeInit[0], OnEdgeSize, VarName, MeshFileID,
                       OnEdgeDecomp, VerticesOnEdgeID);
   if (Err.isFail()) { // not found, try again under older name
      Err = IO::readArray(&VerticesOnEdgeInit[0], OnEdgeSize, VarNameOld,
                          MeshFileID, OnEdgeDecomp, VerticesOnEdgeID);
      CHECK_ERROR_ABORT(Err, "Decomp: error reading VerticesOnEdge");
   }

   VarName    = "CellsOnVertex";
   VarNameOld = "cellsOnVertex";
   int CellsOnVertexID;
   Err = IO::readArray(&CellsOnVertexInit[0], OnVertexSize, VarName, MeshFileID,
                       OnVertexDecomp, CellsOnVertexID);
   if (Err.isFail()) { // not found, try again under older name
      Err = IO::readArray(&CellsOnVertexInit[0], OnVertexSize, VarNameOld,
                          MeshFileID, OnVertexDecomp, CellsOnVertexID);
      CHECK_ERROR_ABORT(Err, "Decomp: error reading CellsOnVertex");
   }

   VarName    = "EdgesOnVertex";
   VarNameOld = "edgesOnVertex";
   int EdgesOnVertexID;
   Err = IO::readArray(&EdgesOnVertexInit[0], OnVertexSize, VarName, MeshFileID,
                       OnVertexDecomp, EdgesOnVertexID);
   if (Err.isFail()) { // not found, try again under older name
      Err = IO::readArray(&EdgesOnVertexInit[0], OnVertexSize, VarNameOld,
                          MeshFileID, OnVertexDecomp, EdgesOnVertexID);
      CHECK_ERROR_ABORT(Err, "Decomp: error reading EdgesOnVertex");
   }

   // Initial decompositions are no longer needed so remove them now
   IO::destroyDecomp(OnCellDecomp);
   IO::destroyDecomp(OnEdgeDecomp);
   IO::destroyDecomp(OnEdgeDecomp2);
   IO::destroyDecomp(OnVertexDecomp);

} // end readMesh

//------------------------------------------------------------------------------
// Initialize the decomposition and create the default decomposition with
// (currently) one partition per MPI task using selected partition method

void Decomp::init(const std::string &MeshFileName) {

   bool TimerFlag = Pacer::start("Decomp init", 0);
   Error Err; // default successful return code

   I4 InHaloWidth;
   std::string DecompMethodStr;

   // Retrieve options from Config
   Config *OmegaConfig = Config::getOmegaConfig();

   Config DecompConfig("Decomp");
   Err = OmegaConfig->get(DecompConfig);
   CHECK_ERROR_ABORT(Err, "Decomp: Decomp group not found in Config");

   Err = DecompConfig.get("HaloWidth", InHaloWidth);
   CHECK_ERROR_ABORT(Err, "Decomp: HaloWidth not found in Decomp Config");

   Err = DecompConfig.get("DecompMethod", DecompMethodStr);
   CHECK_ERROR_ABORT(Err, "Decomp: DecompMethod not found in Decomp Config");

   PartMethod Method = getPartMethodFromStr(DecompMethodStr);

   // Retrieve the default machine environment
   MachEnv *DefEnv = MachEnv::getDefault();

   // Use one partition per MPI task as the default
   I4 NParts = DefEnv->getNumTasks();

   // Create the default decomposition and set pointer to it
   Decomp::DefaultDecomp = Decomp::create("Default", DefEnv, NParts, Method,
                                          InHaloWidth, MeshFileName);

   TimerFlag = Pacer::stop("Decomp init", 0) && TimerFlag;
   if (!TimerFlag)
      LOG_WARN("Decomp::init: Error in timers");

} // End init decomposition

//------------------------------------------------------------------------------
// Construct a new decomposition across an input MachEnv using
// NPart partitions of the mesh.

Decomp::Decomp(
    const std::string &Name,         //< [in] Name for new decomposition
    const MachEnv *InEnv,            //< [in] MachEnv for the new partition
    I4 NParts,                       //< [in] num of partitions for new decomp
    PartMethod Method,               //< [in] method for partitioning
    I4 InHaloWidth,                  //< [in] width of halo in new decomp
    const std::string &MeshFileName_ //< [in] name of file with mesh info
) {

   bool TimerFlag = Pacer::start("Decomp construct", 1);

   // Retrieve some info on the MPI layout
   MPI_Comm Comm = InEnv->getComm();
   I4 NumTasks   = InEnv->getNumTasks();
   I4 MyTask     = InEnv->getMyTask();

   // Open the mesh file for reading (assume IO has already been initialized)
   TimerFlag = Pacer::start("Decomp read mesh", 2) && TimerFlag;
   int FileID;
   MeshFileName = MeshFileName_;
   IO::openFileRead(FileID, MeshFileName);

   // Read mesh size and connectivity information
   std::vector<I4> CellsOnCellInit;
   std::vector<I4> EdgesOnCellInit;
   std::vector<I4> VerticesOnCellInit;
   std::vector<I4> CellsOnEdgeInit;
   std::vector<I4> EdgesOnEdgeInit;
   std::vector<I4> VerticesOnEdgeInit;
   std::vector<I4> CellsOnVertexInit;
   std::vector<I4> EdgesOnVertexInit;
   HaloWidth = InHaloWidth;

   readMesh(FileID, InEnv, NCellsGlobal, NEdgesGlobal, NVerticesGlobal,
            MaxEdges, MaxCellsOnEdge, VertexDegree, CellsOnCellInit,
            EdgesOnCellInit, VerticesOnCellInit, CellsOnEdgeInit,
            EdgesOnEdgeInit, VerticesOnEdgeInit, CellsOnVertexInit,
            EdgesOnVertexInit);

   // Close file
   IO::closeFile(FileID);
   TimerFlag = Pacer::stop("Decomp read mesh", 2) && TimerFlag;

   // In the case of single task avoid calling a full partitioning routine and
   // just set the needed variables directly. This is done because some METIS
   // functions can raise SIGFPE when numparts == 1 due to division by zero
   // See: https://github.com/KarypisLab/METIS/issues/67
   TimerFlag = Pacer::start("Decomp part cells", 2) && TimerFlag;
   if (NumTasks == 1) {
      partCellsSingleTask();
   } else {
      // Use the mesh adjacency information to create a partition of cells
      switch (Method) { // branch depending on method chosen

      //------------------------------------------------------------------------
      // Serial Metis KWay method
      case PartMethodMetisKWay:

         partCellsMetisKWay(InEnv, CellsOnCellInit);
         break;

      //------------------------------------------------------------------------
      // Parallel ParMetis KWay method
      case PartMethodParMetisKWay:

         partCellsParMetisKWay(InEnv, CellsOnCellInit);
         break;

      //------------------------------------------------------------------------
      // Unknown partitioning method
      default:
         ABORT_ERROR("Decomp: Unknown or unsupported decomposition method");

      } // End switch on Method
   }
   TimerFlag = Pacer::stop("Decomp part cells", 2) && TimerFlag;

   //---------------------------------------------------------------------------
   // Cell partitioning complete. Redistribute the initial XXOnCell arrays
   // to their final locations.
   TimerFlag = Pacer::start("Decomp rearrange cells", 2) && TimerFlag;
   rearrangeCellArrays(InEnv, CellsOnCellInit, EdgesOnCellInit,
                       VerticesOnCellInit);
   TimerFlag = Pacer::stop("Decomp rearrange cells", 2) && TimerFlag;

   // Partition the edges
   TimerFlag = Pacer::start("Decomp part edges", 2) && TimerFlag;
   partEdges(InEnv, CellsOnEdgeInit);
   TimerFlag = Pacer::stop("Decomp part edges", 2) && TimerFlag;

   // Edge partitioning complete. Redistribute the initial XXOnEdge arrays
   // to their final locations.
   TimerFlag = Pacer::start("Decomp rearrange edges", 2) && TimerFlag;
   rearrangeEdgeArrays(InEnv, CellsOnEdgeInit, EdgesOnEdgeInit,
                       VerticesOnEdgeInit);
   TimerFlag = Pacer::stop("Decomp rearrange edges", 2) && TimerFlag;

   // Partition the vertices
   TimerFlag = Pacer::start("Decomp part vertices", 2) && TimerFlag;
   partVertices(InEnv, CellsOnVertexInit);
   TimerFlag = Pacer::stop("Decomp part vertices", 2) && TimerFlag;

   // Vertex partitioning complete. Redistribute the initial XXOnVertex arrays
   // to their final locations.
   TimerFlag = Pacer::start("Decomp rearrange vertices", 2) && TimerFlag;
   rearrangeVertexArrays(InEnv, CellsOnVertexInit, EdgesOnVertexInit);
   TimerFlag = Pacer::stop("Decomp rearrange vertices", 2) && TimerFlag;

   // Convert global addresses to local addresses. Create the global to
   // local address ordered maps to simplify and optimize searches.
   // Invalid/non-existent edges have been assigned NXxGlobal+1 and we want
   // to map that to the local NXxAll+1 (NXxSize). We insert that value
   // first in the map so that later attempts to change will be ignored.

   TimerFlag = Pacer::start("Decomp construct global to local", 2) && TimerFlag;
   std::map<I4, I4> GlobToLocCell;
   GlobToLocCell[NCellsGlobal + 1] = NCellsAll;
   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      GlobToLocCell[CellIDH(Cell)] = Cell;
   }

   std::map<I4, I4> GlobToLocEdge;
   GlobToLocEdge[NEdgesGlobal + 1] = NEdgesAll;
   for (int Edge = 0; Edge < NEdgesAll; ++Edge) {
      GlobToLocEdge[EdgeIDH(Edge)] = Edge;
   }

   std::map<I4, I4> GlobToLocVrtx;
   GlobToLocVrtx[NVerticesGlobal + 1] = NVerticesAll;
   for (int Vrtx = 0; Vrtx < NVerticesAll; ++Vrtx) {
      GlobToLocVrtx[VertexIDH(Vrtx)] = Vrtx;
   }

   // CellsOnCell
   for (int Cell = 0; Cell < NCellsSize; ++Cell) {
      for (int NbrCell = 0; NbrCell < MaxEdges; ++NbrCell) {
         I4 GlobID = CellsOnCellH(Cell, NbrCell);
         // In some cases (eg in halo regions) there are still
         // some remaining CellIDs that are not in the local owned
         // or halo domain and are on another task.  Catch that case here
         auto it = GlobToLocCell.find(GlobID);
         I4 LocalAdd;
         if (it != GlobToLocCell.end()) {
            LocalAdd = it->second; // retrieve map entry if exists
         } else {
            LocalAdd = NCellsAll; // otherwise set to bndy/unknown address
            GlobToLocCell[GlobID] = NCellsAll; // add an entry to map
         }
         CellsOnCellH(Cell, NbrCell) = LocalAdd;
      }
   }

   // EdgesOnCell
   for (int Cell = 0; Cell < NCellsSize; ++Cell) {
      for (int Edge = 0; Edge < MaxEdges; ++Edge) {
         I4 GlobID = EdgesOnCellH(Cell, Edge);
         // In some cases (eg in halo regions) there are still
         // some remaining EdgeIDs that are not in the local owned
         // or halo domain and are on another task.  Catch that case here
         auto it = GlobToLocEdge.find(GlobID);
         I4 LocalAdd;
         if (it != GlobToLocEdge.end()) {
            LocalAdd = it->second; // retrieve map entry if exists
         } else {
            LocalAdd = NEdgesAll; // otherwise set to bndy/unknown address
            GlobToLocEdge[GlobID] = NEdgesAll; // add an entry to map
         }
         EdgesOnCellH(Cell, Edge) = LocalAdd;
      }
   }

   // VerticesOnCell
   for (int Cell = 0; Cell < NCellsSize; ++Cell) {
      for (int Vrtx = 0; Vrtx < MaxEdges; ++Vrtx) {
         I4 GlobID = VerticesOnCellH(Cell, Vrtx);
         // In some cases (eg in halo regions) there are still
         // some remaining EdgeIDs that are not in the local owned
         // or halo domain and are on another task.  Catch that case here
         auto it = GlobToLocVrtx.find(GlobID);
         I4 LocalAdd;
         if (it != GlobToLocVrtx.end()) {
            LocalAdd = it->second; // retrieve map entry if exists
         } else {
            LocalAdd = NVerticesAll; // otherwise set to bndy/unknown address
            GlobToLocVrtx[GlobID] = NVerticesAll; // add an entry to map
         }
         VerticesOnCellH(Cell, Vrtx) = LocalAdd;
      }
   }

   // CellsOnEdge
   for (int Edge = 0; Edge < NEdgesSize; ++Edge) {
      for (int Cell = 0; Cell < MaxCellsOnEdge; ++Cell) {
         I4 GlobID = CellsOnEdgeH(Edge, Cell);
         auto it   = GlobToLocCell.find(GlobID);
         I4 LocalAdd;
         if (it != GlobToLocCell.end()) {
            LocalAdd = it->second; // retrieve map entry if exists
         } else {
            LocalAdd = NCellsAll; // otherwise set to bndy/unknown address
            GlobToLocCell[GlobID] = NCellsAll; // add an entry to map
         }
         CellsOnEdgeH(Edge, Cell) = LocalAdd;
      }
   }

   // EdgesOnEdge
   for (int Edge = 0; Edge < NEdgesSize; ++Edge) {
      for (int NbrEdge = 0; NbrEdge < 2 * MaxEdges; ++NbrEdge) {
         I4 GlobID = EdgesOnEdgeH(Edge, NbrEdge);
         auto it   = GlobToLocEdge.find(GlobID);
         I4 LocalAdd;
         if (it != GlobToLocEdge.end()) {
            LocalAdd = it->second; // retrieve map entry if exists
         } else {
            LocalAdd = NEdgesAll; // otherwise set to bndy/unknown address
            GlobToLocEdge[GlobID] = NEdgesAll; // add an entry to map
         }
         EdgesOnEdgeH(Edge, NbrEdge) = LocalAdd;
      }
   }

   // VerticesOnEdge
   for (int Edge = 0; Edge < NEdgesSize; ++Edge) {
      for (int Vrtx = 0; Vrtx < 2; ++Vrtx) {
         I4 GlobID = VerticesOnEdgeH(Edge, Vrtx);
         auto it   = GlobToLocVrtx.find(GlobID);
         I4 LocalAdd;
         if (it != GlobToLocVrtx.end()) {
            LocalAdd = it->second; // retrieve map entry if exists
         } else {
            LocalAdd = NVerticesAll; // otherwise set to bndy/unknown address
            GlobToLocVrtx[GlobID] = NVerticesAll; // add an entry to map
         }
         VerticesOnEdgeH(Edge, Vrtx) = LocalAdd;
      }
   }

   // CellsOnVertex
   for (int Vrtx = 0; Vrtx < NVerticesSize; ++Vrtx) {
      for (int Cell = 0; Cell < VertexDegree; ++Cell) {
         I4 GlobID = CellsOnVertexH(Vrtx, Cell);
         auto it   = GlobToLocCell.find(GlobID);
         I4 LocalAdd;
         if (it != GlobToLocCell.end()) {
            LocalAdd = it->second; // retrieve map entry if exists
         } else {
            LocalAdd = NCellsAll; // otherwise set to bndy/unknown address
            GlobToLocCell[GlobID] = NCellsAll; // add an entry to map
         }
         CellsOnVertexH(Vrtx, Cell) = LocalAdd;
      }
   }

   // EdgesOnVertex
   for (int Vrtx = 0; Vrtx < NVerticesSize; ++Vrtx) {
      for (int Edge = 0; Edge < VertexDegree; ++Edge) {
         I4 GlobID = EdgesOnVertexH(Vrtx, Edge);
         auto it   = GlobToLocEdge.find(GlobID);
         I4 LocalAdd;
         if (it != GlobToLocEdge.end()) {
            LocalAdd = it->second; // retrieve map entry if exists
         } else {
            LocalAdd = NEdgesAll; // otherwise set to bndy/unknown address
            GlobToLocEdge[GlobID] = NEdgesAll; // add an entry to map
         }
         EdgesOnVertexH(Vrtx, Edge) = LocalAdd;
      }
   }
   TimerFlag = Pacer::stop("Decomp construct global to local", 2) && TimerFlag;

   // Create device copies of all arrays

   TimerFlag  = Pacer::start("Decomp construct device copy", 2) && TimerFlag;
   NCellsHalo = createDeviceMirrorCopy(NCellsHaloH);
   CellID     = createDeviceMirrorCopy(CellIDH);
   CellLoc    = createDeviceMirrorCopy(CellLocH);

   NEdgesHalo = createDeviceMirrorCopy(NEdgesHaloH);
   EdgeID     = createDeviceMirrorCopy(EdgeIDH);
   EdgeLoc    = createDeviceMirrorCopy(EdgeLocH);

   NVerticesHalo = createDeviceMirrorCopy(NVerticesHaloH);
   VertexID      = createDeviceMirrorCopy(VertexIDH);
   VertexLoc     = createDeviceMirrorCopy(VertexLocH);

   CellsOnCell    = createDeviceMirrorCopy(CellsOnCellH);
   EdgesOnCell    = createDeviceMirrorCopy(EdgesOnCellH);
   VerticesOnCell = createDeviceMirrorCopy(VerticesOnCellH);
   NEdgesOnCell   = createDeviceMirrorCopy(NEdgesOnCellH);

   CellsOnEdge    = createDeviceMirrorCopy(CellsOnEdgeH);
   EdgesOnEdge    = createDeviceMirrorCopy(EdgesOnEdgeH);
   VerticesOnEdge = createDeviceMirrorCopy(VerticesOnEdgeH);
   NEdgesOnEdge   = createDeviceMirrorCopy(NEdgesOnEdgeH);

   CellsOnVertex = createDeviceMirrorCopy(CellsOnVertexH);
   EdgesOnVertex = createDeviceMirrorCopy(EdgesOnVertexH);
   TimerFlag     = Pacer::stop("Decomp construct device copy", 2) && TimerFlag;
   TimerFlag     = Pacer::stop("Decomp construct", 1) && TimerFlag;
   if (!TimerFlag)
      LOG_WARN("Decomp constructor: Error encounterd in timers");
} // end decomposition constructor

//------------------------------------------------------------------------------
// Creates a new decomposition using the constructor and puts it in the
// AllDecomps map
Decomp *Decomp::create(
    const std::string &Name,        //< [in] Name for new decomposition
    const MachEnv *Env,             //< [in] MachEnv for the new partition
    I4 NParts,                      //< [in] num of partitions for new decomp
    PartMethod Method,              //< [in] method for partitioning
    I4 HaloWidth,                   //< [in] width of halo in new decomp
    const std::string &MeshFileName //< [in] name of file with mesh info
) {

   bool TimerFlag = Pacer::start("Decomp create", 1);
   // Check to see if a decomposition of the same name already exists and
   // if so, exit with an error
   if (AllDecomps.find(Name) != AllDecomps.end()) {
      ABORT_ERROR("Attempted to create a Decomp with name {} but a Decomp of "
                  "that name already exists",
                  Name);
   }

   // create a new decomp on the heap and put it in a map of
   // unique_ptrs, which will manage its lifetime
   auto *NewDecomp =
       new Decomp(Name, Env, NParts, Method, HaloWidth, MeshFileName);
   AllDecomps.emplace(Name, NewDecomp);

   TimerFlag = Pacer::stop("Decomp create", 1) && TimerFlag;
   if (!TimerFlag)
      LOG_WARN("Decomp::create: Error encountered in timers");
   return NewDecomp;
} // end Decomp create

// Destructor
//------------------------------------------------------------------------------
// Destroys a decomposition and deallocates all arrays

Decomp::~Decomp() {

   // No operations needed, Kokkos arrays removed when no longer in scope

} // end decomp destructor

//------------------------------------------------------------------------------
// Removes a decomposition from list and destroys it

void Decomp::erase(std::string InName // [in] name of decomp to remove
) {

   AllDecomps.erase(InName); // removes the decomp from the list (map) and in
                             // the process, calls the destructor

} // end decomp erase

//------------------------------------------------------------------------------
// Removes all decompositions to clean up before exit

void Decomp::clear() {

   AllDecomps.clear();      // removes all decomps from the list (map) and in
                            // the process, calls the destructors for each
   DefaultDecomp = nullptr; // prevent dangling pointer

} // end decomp clear

// Retrieval functions
//------------------------------------------------------------------------------
// Get default decomposition
Decomp *Decomp::getDefault() { return Decomp::DefaultDecomp; }

//------------------------------------------------------------------------------
// Get decomposition by name
Decomp *Decomp::get(const std::string Name ///< [in] Name of environment
) {

   // look for an instance of this name
   auto it = AllDecomps.find(Name);

   // if not found abort with an error
   if (it == AllDecomps.end())
      ABORT_ERROR(
          "Decomp::get: Attempt to retrieve non-existent Decomposition: {}",
          Name);

   // return the decomposition pointer
   return it->second.get();

} // end get Decomposition

//------------------------------------------------------------------------------
// Trivially partition cells in the case of single task
// It sets the NCells sizes (NCellsOwned,
// NCellsHalo array, NCellsAll and NCellsSize) and the final CellID
// and CellLoc arrays
void Decomp::partCellsSingleTask() {

   NCellsOwned = NCellsGlobal;

   NCellsHaloH = HostArray1DI4("NCellsHalo", HaloWidth);
   for (int Halo = 0; Halo < HaloWidth; ++Halo) {
      NCellsHaloH(Halo) = NCellsGlobal;
   }

   NCellsAll  = NCellsGlobal;
   NCellsSize = NCellsGlobal + 1;

   CellIDH  = HostArray1DI4("CellID", NCellsSize);
   CellLocH = HostArray2DI4("CellLoc", NCellsSize, 2);
   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      CellIDH(Cell)     = Cell + 1;
      CellLocH(Cell, 0) = 0;
      CellLocH(Cell, 1) = Cell;
   }
} // end function partCellsSingleTask

//------------------------------------------------------------------------------
// Partition the cells using the serial Metis KWay method
// After this partitioning, the decomposition class member CellID and
// CellLocator arrays have been set as well as the class NCells size variables
// (owned, halo, all).

void Decomp::partCellsMetisKWay(
    const MachEnv *InEnv, // [in] input machine environment with MPI info
    const std::vector<I4> &CellsOnCellInit // [in] cell nbrs in linear distrb
) {

   // Retrieve some info on the MPI layout
   MPI_Comm Comm = InEnv->getComm();
   I4 NumTasks   = InEnv->getNumTasks();
   I4 MyTask     = InEnv->getMyTask();

   // Allocate adjacency arrays and a buffer for portions of the CellsOnCell
   // array.
   std::vector<idx_t> AdjAdd(NCellsGlobal + 1, 0);
   std::vector<idx_t> Adjacency(2 * NEdgesGlobal, 0);
   I4 NCellsChunk     = (NCellsGlobal - 1) / NumTasks + 1;
   I4 CellsOnCellSize = NCellsChunk * MaxEdges;
   std::vector<I4> CellsOnCellBuf(CellsOnCellSize, 0);

   // This is an address counter needed to keep track of the starting
   // address for each cell in the packed adjacency array.
   I4 Add = 0;

   // One at a time, each task broadcasts its portion of the CellsOnCell data
   // and then unpacks it into the adjacency array in the packed form needed
   // by METIS/ParMETIS
   bool TimerFlag = Pacer::start("Gather adjacency", 2);
   for (int ITask = 0; ITask < NumTasks; ++ITask) {

      // If it is this task's turn, pack up the CellsOnCell data and broadcast
      // to other tasks.
      if (MyTask == ITask) {
         for (int n = 0; n < CellsOnCellSize; ++n) {
            CellsOnCellBuf[n] = CellsOnCellInit[n];
         } // end loop CellsOnCell
      } // end if this is MyTask
      Broadcast(CellsOnCellBuf, InEnv, ITask);

      // Create the adjacency graph by aggregating the individual
      // chunks. Prune edges that don't have neighbors.
      for (int Cell = 0; Cell < NCellsChunk; ++Cell) {

         I4 CellGlob = ITask * NCellsChunk + Cell; // global cell ID

         // when chunks do not divide evenly this happens and we break out
         if (CellGlob >= NCellsGlobal)
            break;

         AdjAdd[CellGlob] = Add; // start add for cell in Adjacency array
         for (int Edge = 0; Edge < MaxEdges; ++Edge) {
            I4 BufAdd  = Cell * MaxEdges + Edge;
            I4 NbrCell = CellsOnCellBuf[BufAdd];
            // Skip edges with no neighbors
            if (validCellID(NbrCell)) {
               // switch to 0-based indx
               Adjacency[Add] = NbrCell - 1;
               ++Add; // increment address counter
            }
         }
      } // end cell loop for buffer
   } // end task loop

   AdjAdd[NCellsGlobal] = Add; // Add the ending address
   TimerFlag            = Pacer::stop("Gather adjacency", 2) && TimerFlag;

   // Set up remaining partitioning variables

   // NConstraints is the number of balancing constraints, mostly for
   // use when multiple vertex weights are assigned. Must be at least 1.
   // We do not yet support weighted partitions
   idx_t NConstraints = 1;

   // Arrays needed for weighted decompositions. If no weighting used
   // set pointers to null.
   idx_t *VrtxWgtPtr{nullptr};
   idx_t *EdgeWgtPtr{nullptr};
   idx_t *VrtxSize{nullptr};

   // Use default metis options
   idx_t *Options{nullptr};

   // These are for multi-constraint partitions where the vertex weight
   // has to be distributed among the multiple constraints.
   // We do not use them so set them to null
   real_t *TpWgts{nullptr};
   real_t *Ubvec{nullptr};

   // Results are stored in a new partition array which returns
   // the processor (partition) assigned to the cell (vrtx)
   // The routine also returns the number of edge cuts in the partition
   std::vector<idx_t> CellTask(NCellsGlobal);
   idx_t Edgecut = 0;

   // Convert to idx_t from Omega::I4, in case these aren't the same
   idx_t NCellsMetis   = NCellsGlobal;
   idx_t NumTasksMetis = NumTasks;

   // Call METIS routine to partition the mesh
   // METIS routines are C code that expect pointers, so we use the
   // idiom &Var[0] to extract the pointer to the data in std::vector
   TimerFlag    = Pacer::start("Metis partitioning", 2) && TimerFlag;
   int MetisErr = METIS_PartGraphKway(&NCellsMetis, &NConstraints, &AdjAdd[0],
                                      &Adjacency[0], VrtxWgtPtr, VrtxSize,
                                      EdgeWgtPtr, &NumTasksMetis, TpWgts, Ubvec,
                                      Options, &Edgecut, &CellTask[0]);

   if (MetisErr != METIS_OK)
      ABORT_ERROR("Decomp: Error in ParMETIS");
   TimerFlag = Pacer::stop("Metis partitioning", 2) && TimerFlag;

   // Determine the initial sizes needed by address arrays

   std::vector<I4> TaskCount(NumTasks, 0);
   for (int Cell = 0; Cell < NCellsGlobal; ++Cell) {
      I4 TaskLoc  = CellTask[Cell];
      I4 LocalAdd = TaskCount[TaskLoc];
      ++TaskCount[TaskLoc]; // increment number of cells assigned to task
   }
   NCellsOwned = TaskCount[MyTask];

   // Assign the full address (TaskID, local index) for each cell. These
   // are implicity sorted in CellID order.
   // During this process we also create an ordered list of all local
   // (owned+halo) cells using std::set for later use in halo setup

   std::vector<I4> CellLocAll(2 * NCellsGlobal);
   std::vector<I4> CellIDTmp(NCellsOwned,
                             NCellsGlobal + 1);    // initial size only
   std::vector<I4> CellLocTmp(2 * NCellsOwned, 0); // will grow when halo added
   std::set<I4> CellsInList; // list of unique cells in local owned, halo

   // reset local address counter
   for (int ITask = 0; ITask < NumTasks; ++ITask)
      TaskCount[ITask] = 0;

   for (int Cell = 0; Cell < NCellsGlobal; ++Cell) {
      I4 TaskLoc  = CellTask[Cell];
      I4 LocalAdd = TaskCount[TaskLoc];
      ++TaskCount[TaskLoc]; // increment number of cells assigned to task
      CellLocAll[2 * Cell]     = TaskLoc;  // Task location
      CellLocAll[2 * Cell + 1] = LocalAdd; // local address within task
      // If this cell is on the local task, store as a local cell
      if (TaskLoc == MyTask) {
         CellIDTmp[LocalAdd] = Cell + 1; // IDs are 1-based
         CellsInList.insert(Cell + 1);
         CellLocTmp[2 * LocalAdd]     = TaskLoc;
         CellLocTmp[2 * LocalAdd + 1] = LocalAdd;
      } // end if my task
   } // end loop over all cells

   // Find and add the halo cells to the cell list. Here we use the
   // adjacency array to find the active neighbor cells and store if they
   // are not already owned by the task. We use the std::set container
   // to automatically sort each halo layer by cellID.
   I4 CellLocStart = 0;
   I4 CellLocEnd   = NCellsOwned - 1;
   I4 CurSize      = NCellsOwned;
   HostArray1DI4 NCellsHaloTmp("NCellsHalo", HaloWidth);
   std::set<I4> HaloList;
   // Loop over each halo layer
   for (int Halo = 0; Halo < HaloWidth; ++Halo) {
      // For each of the local cells in the previous level, retrieve
      // the GlobalID
      HaloList.clear(); // reset list for this halo layer
      I4 NbrGlob;       // use for global address of each neighbor
      I4 NbrID;         // global cell ID of each neighbor
      for (int CellLoc = CellLocStart; CellLoc <= CellLocEnd; ++CellLoc) {
         I4 CellGlob = CellIDTmp[CellLoc] - 1;
         // Loop through the neighbor cells (stored in the adjacency arrays)
         // to see if they are remote or local
         I4 NbrStart = AdjAdd[CellGlob];     // start add in adjacency array
         I4 NbrEnd   = AdjAdd[CellGlob + 1]; // end+1 add in adjacency array
         for (int NbrCell = NbrStart; NbrCell < NbrEnd; ++NbrCell) {
            // Get global ID for each neighbor cell
            NbrGlob = Adjacency[NbrCell];
            NbrID   = NbrGlob + 1; // Cell IDs are 1-based

            // Check to see if this cell is on a remote task and has not
            // already been added
            I4 NbrTask = CellLocAll[2 * NbrGlob];
            if (NbrTask != MyTask) {
               // only add to the list if the cell is not already listed
               // (the find function will return the end iterator if not found)
               if (CellsInList.find(NbrID) == CellsInList.end()) {
                  HaloList.insert(NbrID);
                  CellsInList.insert(NbrID);
               } // end search for existing entry
            } // end if not on task

         } // end loop over neighbors

      } // end loop over previous layer

      // Extract the values from the Halo list into the ID and location
      // vectors. Extend size of ID, Loc arrays.
      I4 HaloAdd = CellLocEnd;
      CurSize += HaloList.size();
      NCellsHaloTmp(Halo) = CurSize;
      CellIDTmp.resize(CurSize);
      CellLocTmp.resize(2 * CurSize);

      for (auto IHalo = HaloList.begin(); IHalo != HaloList.end(); IHalo++) {
         ++HaloAdd;
         NbrID              = *IHalo; // the neighbor cell id from the HaloList
         NbrGlob            = NbrID - 1; // global address of this neighbor
         CellIDTmp[HaloAdd] = NbrID;
         CellLocTmp[2 * HaloAdd]     = CellLocAll[2 * NbrGlob];
         CellLocTmp[2 * HaloAdd + 1] = CellLocAll[2 * NbrGlob + 1];
      }

      // Reset for next halo layer
      CellLocStart = CellLocEnd + 1;
      CellLocEnd   = NCellsHaloTmp(Halo) - 1;
   }
   NCellsAll  = NCellsHaloTmp(HaloWidth - 1);
   NCellsSize = NCellsAll + 1; // extra entry to store boundary/undefined value

   // The cell decomposition is now complete, copy the information
   // into the final locations as class members on host (copy to device later)

   NCellsHaloH = NCellsHaloTmp;

   // Copy global ID for each cell, both owned and halo.
   // Copy cell location (task, local add) for each cell, both owned and halo

   HostArray1DI4 CellIDHTmp("CellID", NCellsSize);
   HostArray2DI4 CellLocHTmp("CellLoc", NCellsSize, 2);
   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      CellIDHTmp(Cell)     = CellIDTmp[Cell];
      CellLocHTmp(Cell, 0) = CellLocTmp[2 * Cell];     // task owning this cell
      CellLocHTmp(Cell, 1) = CellLocTmp[2 * Cell + 1]; // local address on task
   }
   CellIDH  = CellIDHTmp;
   CellLocH = CellLocHTmp;

   // All done
   if (!TimerFlag)
      LOG_WARN("Decomp::partCellsKWay: Error in timers");
   return;

} // end function partCellsMetisKWay

//------------------------------------------------------------------------------
// Partition the cells using the parallel ParMetis KWay method
// After this partitioning, the decomposition class member CellID and
// CellLocator arrays have been set as well as the class NCells size variables
// (owned, halo, all).

void Decomp::partCellsParMetisKWay(
    const MachEnv *InEnv, // [in] input machine environment with MPI info
    const std::vector<I4> &CellsOnCellInit // [in] cell nbrs in linear distrb
) {

   // Retrieve some info on the MPI layout
   MPI_Comm Comm = InEnv->getComm();
   I4 NumTasks   = InEnv->getNumTasks();
   I4 MyTask     = InEnv->getMyTask();

   // On entry, the cells are distributed linearly across tasks so the metis
   // array describing the distribution is just the beginning of each chunk
   // (Note: We use cell centers as the graph vertices in Metis)
   I4 NCellsChunk = (NCellsGlobal - 1) / NumTasks + 1;
   std::vector<idx_t> VrtxDist(NumTasks + 1, 0);
   for (int IChunk = 0; IChunk < NumTasks; ++IChunk) {
      VrtxDist[IChunk] = IChunk * NCellsChunk;
   }
   VrtxDist[NumTasks] = NCellsGlobal; // last address
   I4 NCellsLocal     = NCellsChunk;  // most processors have full chunk
   if (MyTask == NumTasks - 1) { // if uneven distribution, must adjust last
      I4 NCellsStart = (NumTasks - 1) * NCellsChunk;
      NCellsLocal    = NCellsGlobal - NCellsStart;
   }

   // Allocate adjacency arrays
   I4 AdjSize = NCellsLocal * MaxEdges;
   std::vector<idx_t> AdjAdd(NCellsLocal + 1, 0);
   std::vector<idx_t> Adjacency(AdjSize, 0);

   // This is an address counter needed to keep track of the starting
   // address for each cell in the packed adjacency array.
   I4 Add = 0;

   // Now we pack the local adjacency array into the form needed by ParMetis

   // Create the adjacency graph by extracting the cells on cell info
   // and packing into the adjacency array and pruning special entries
   for (int Cell = 0; Cell < NCellsLocal; ++Cell) {

      AdjAdd[Cell] = Add; // starting address for neighbors of this local cell

      for (int Edge = 0; Edge < MaxEdges; ++Edge) {
         I4 VecAdd  = Cell * MaxEdges + Edge; // vector add into CellsOnCell
         I4 NbrCell = CellsOnCellInit[VecAdd];
         // Only include valid neighbors (eg non-existent edges or edges along
         // boundaries with no neighbors)
         if (NbrCell > 0 && NbrCell <= NCellsGlobal) {
            // switch to 0-based indx
            Adjacency[Add] = NbrCell - 1;
            ++Add; // increment address counter
         }
      } // end edge loop
   } // end cell loop
   // Add last end address
   AdjAdd[NCellsLocal] = Add;
   Adjacency.resize(Add);

   // Set up remaining partitioning variables

   // Options array
   // Options[0] = 0 for default (1 otherwise), Options[1] = 1 for timing info
   // Options[2] = random number seed
   std::vector<idx_t> Options(3, 0); // default options

   // Arrays needed for weighted decompositions. If no weighting used
   // set pointers to null.
   idx_t WgtFlag = 0; // no weights (1 for edge wgts, 2 for vrtx, 3 for both)
   idx_t *VrtxWgtPtr{nullptr};
   idx_t *EdgeWgtPtr{nullptr};

   // NConstraints is the number of balancing constraints, mostly for
   // use when multiple vertex weights are assigned. Must be at least 1.
   // We do not yet support weighted partitions
   idx_t NConstraints = 1;

   // These are for multi-constraint partitions where the vertex weight
   // has to be distributed among the multiple constraints.
   // We do not use them so set them to an even distribution of weight
   std::vector<real_t> TpWgts(NConstraints * NumTasks, 1.0 / NumTasks);
   std::vector<real_t> Ubvec(NConstraints, 1.05);
   real_t *TpWgtsPtr = TpWgts.data();
   real_t *UbvecPtr  = Ubvec.data();

   // Results are stored in a new partition array which returns
   // the processor (partition) assigned to the cell (vrtx)
   // The routine also returns the number of edge cuts in the partition
   std::vector<idx_t> CellTask(NCellsChunk);
   idx_t Edgecut = 0;

   // Additional info for Metis
   idx_t NumFlagMetis  = 0;        // C-style base-0 indexing
   idx_t NumPartsMetis = NumTasks; // number of partitions = num MPI tasks

   // Call ParMETIS routine to partition the mesh
   // METIS routines are C code that expect pointers, so we use the
   // idiom &Var[0] to extract the pointer to any data in std::vector
   bool TimerFlag = Pacer::start("Metis partitioning", 2);
   int MetisErr   = ParMETIS_V3_PartKway(
       &VrtxDist[0], &AdjAdd[0], &Adjacency[0], VrtxWgtPtr, EdgeWgtPtr,
       &WgtFlag, &NumFlagMetis, &NConstraints, &NumPartsMetis, TpWgtsPtr,
       UbvecPtr, &Options[0], &Edgecut, &CellTask[0], &Comm);

   TimerFlag = Pacer::stop("Metis partitioning", 2);
   if (MetisErr != METIS_OK)
      ABORT_ERROR("Decomp: Error in ParMETIS");

   // The location of each cell is now distributed across tasks, so we
   // need to communicate that information to their final locations. We do
   // this by having each task broadcast its piece of the distribution
   // information and gathering the cellIDs owned by this task. We also
   // gather some connectivity info to create the first level of ghost
   // cells.

   // Create message buffer for cell location and the nbr cell IDs
   // Create temporary vector to hold CellIDs, locators and list of Halo cells
   TimerFlag  = Pacer::start("partCellsGatherOwned", 2);
   I4 BufSize = NCellsChunk * (1 + MaxEdges);
   std::vector<I4> MsgBuf(BufSize);
   std::vector<I4> CellIDTmp;    // Global ID for all cells needed by this task
   std::vector<I4> CellTaskTmp;  // Task owner for remote halo cells
   std::vector<I4> CellLocalTmp; // Local index for remote halo cells
   std::set<I4> HaloList;        // List of halo/nbr points for next pass
   std::vector<I4> TaskCounter(NumTasks + 1, 0);
   I4 CellCount = 0; // counter for number of local owned and halo cells

   // First make a pass to determine owned cells and collect global IDs for
   // the first halo layer
   for (int ITask = 0; ITask < NumTasks; ++ITask) {

      // If it is this tasks turn, fill the buffer with cell loc and nbrs
      if (MyTask == ITask) {
         for (int I = 0; I < NCellsLocal; ++I) {
            int BufAdd     = I * (1 + MaxEdges);
            MsgBuf[BufAdd] = CellTask[I];
            for (int J = 0; J < MaxEdges; ++J) {
               int VecAdd             = I * MaxEdges + J;
               MsgBuf[BufAdd + 1 + J] = CellsOnCellInit[VecAdd];
            }
         }
         // handle case where NCells does not divide evenly over tasks
         for (int I = NCellsLocal; I < NCellsChunk; ++I) {
            int BufAdd     = I * (1 + MaxEdges);
            MsgBuf[BufAdd] = NumTasks;
         }
      }

      // broadcast the buffer to all tasks
      TimerFlag = Pacer::start("partCellsOwnerBcast", 3);
      Broadcast(MsgBuf, InEnv, ITask);
      TimerFlag = Pacer::stop("partCellsOwnerBcast", 3);

      // for each cell in the buffer, determine its location (task, local indx)
      // if it's on this task, store as an owned cell and gather nbr info for
      // the first halo level
      for (int BufCell = 0; BufCell < NCellsChunk; ++BufCell) {
         int BufAdd    = BufCell * (1 + MaxEdges); // location of task owner
         int TaskLoc   = MsgBuf[BufAdd];
         int LocalIndx = TaskCounter[TaskLoc];
         if (TaskLoc == MyTask) { // this cell owned by this task
            int GlobalID = ITask * NCellsChunk + BufCell + 1; // IDs are 1-based
            CellIDTmp.push_back(GlobalID);
            CellTaskTmp.push_back(TaskLoc);
            CellLocalTmp.push_back(LocalIndx);
            ++CellCount;
            // collect nbr ids for this cell and create the halo list
            for (int Nbr = 0; Nbr < MaxEdges; ++Nbr) {
               int NbrID = MsgBuf[BufAdd + Nbr + 1];
               if (validCellID(NbrID)) // active neighbor ID
                  HaloList.insert(NbrID);
            }
         }
         // increment the task counter for next pass
         TaskCounter[TaskLoc] = LocalIndx + 1;
      }
   }

   // Set the number of owned cells and add halo points to CellID list
   NCellsOwned = CellCount;
   for (auto Iter = HaloList.begin(); Iter != HaloList.end(); ++Iter) {
      int HaloID = *Iter;
      // add Halo cell if not already an owned cell
      if (std::find(CellIDTmp.begin(), CellIDTmp.end(), HaloID) ==
          CellIDTmp.end()) { // not found
         CellIDTmp.push_back(HaloID);
         ++CellCount;
      }
   }
   // Add number of owned+halo points for this first halo level
   HostArray1DI4 NCellsHaloTmp("NCellsHalo", HaloWidth);
   NCellsHaloTmp(0) = CellCount;

   TimerFlag = Pacer::stop("partCellsGatherOwned", 2);

   // For each halo layer, we find the location of the current layer and
   // gather the global address of the next layer. The process is similar to
   // the owned process above, broadcasting and extracting needed info
   TimerFlag      = Pacer::start("partCellsGatherHalos", 2);
   I4 SrchIDBegin = NCellsOwned; // initialize search range for ID search
   I4 SrchIDEnd   = CellCount;
   for (int HaloLvl = 0; HaloLvl < HaloWidth; ++HaloLvl) {

      // reset task counter and halo list
      for (int ITask = 0; ITask <= NumTasks; ++ITask)
         TaskCounter[ITask] = 0;
      HaloList.clear();

      // Loop over tasks to broadcast the partition info on each
      for (int ITask = 0; ITask < NumTasks; ++ITask) {

         // If it is this tasks turn, fill the buffer with cell loc and nbrs
         if (MyTask == ITask) {
            for (int I = 0; I < NCellsLocal; ++I) {
               int BufAdd     = I * (1 + MaxEdges);
               MsgBuf[BufAdd] = CellTask[I];
               for (int J = 0; J < MaxEdges; ++J) {
                  int VecAdd             = I * MaxEdges + J;
                  MsgBuf[BufAdd + J + 1] = CellsOnCellInit[VecAdd];
               }
            }
            // handle case where NCells does not divide evenly over tasks
            for (int I = NCellsLocal; I < NCellsChunk; ++I) {
               int BufAdd     = I * (1 + MaxEdges);
               MsgBuf[BufAdd] = NumTasks; // fill with invalid task ID
            }
         }

         // broadcast the buffer to all tasks
         TimerFlag = Pacer::start("partCellsHaloBCast", 3);
         Broadcast(MsgBuf, InEnv, ITask);
         TimerFlag = Pacer::stop("partCellsHaloBCast", 3);

         // for each cell in the buffer, determine location (task, local indx)
         // if it's needed by the halo, store the info and gather the nbrs for
         // the next halo level
         TimerFlag = Pacer::start("partCellsHaloSearch", 3);
         for (int BufCell = 0; BufCell < NCellsChunk; ++BufCell) {
            int GlobalID         = ITask * NCellsChunk + BufCell + 1;
            int BufAdd           = BufCell * (1 + MaxEdges);
            int TaskLoc          = MsgBuf[BufAdd];
            int LocalIndx        = TaskCounter[TaskLoc];
            TaskCounter[TaskLoc] = LocalIndx + 1;
            // If tasks did not divide mesh evenly, the last group is
            // padded with a task outside the range, so skip these
            if (TaskLoc < 0 or TaskLoc >= NumTasks)
               continue;
            // search this halo level to see if this cell is needed
            for (int SrchIndx = SrchIDBegin; SrchIndx < SrchIDEnd; ++SrchIndx) {
               if (CellIDTmp[SrchIndx] == GlobalID) { // found
                  CellTaskTmp.push_back(TaskLoc);
                  CellLocalTmp.push_back(LocalIndx);
                  // collect nbr ids for this cell and create the halo list
                  for (int Nbr = 0; Nbr < MaxEdges; ++Nbr) {
                     int NbrID = MsgBuf[BufAdd + Nbr + 1];
                     if (validCellID(NbrID)) // active neighbor ID
                        HaloList.insert(NbrID);
                  }
                  break;
               }
            }
         }
         TimerFlag = Pacer::stop("partCellsHaloSearch", 3);
      } // end task loop

      // Add new halo points to CellID vector (except last pass) and update
      // count indices
      if (HaloLvl < HaloWidth - 1) {
         for (auto Iter = HaloList.begin(); Iter != HaloList.end(); ++Iter) {
            int HaloID = *Iter;
            // add Halo cell if not already an owned cell
            if (std::find(CellIDTmp.begin(), CellIDTmp.end(), HaloID) ==
                CellIDTmp.end()) { // not found
               CellIDTmp.push_back(HaloID);
               ++CellCount;
            }
         } // end halo list loop
         // Set new halo level index and reset search range for next halo lvl
         SrchIDBegin                = SrchIDEnd;
         SrchIDEnd                  = CellCount;
         NCellsHaloTmp(HaloLvl + 1) = CellCount;
      } // endif halolvl

   } // end loop halo levels
   TimerFlag = Pacer::stop("partCellsGatherHalos", 2);

   // The cell decomposition is now complete, copy the information
   // into the final locations as class members on host (copy to device later)

   NCellsAll   = CellCount;
   NCellsSize  = NCellsAll + 1;
   NCellsHaloH = NCellsHaloTmp;

   // Copy global ID for each cell, both owned and halo.
   // Copy cell location (task, local add) for each cell, both owned and halo

   HostArray1DI4 CellIDHTmp("CellID", NCellsSize);
   HostArray2DI4 CellLocHTmp("CellLoc", NCellsSize, 2);
   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      CellIDHTmp(Cell)     = CellIDTmp[Cell];
      CellLocHTmp(Cell, 0) = CellTaskTmp[Cell];  // task owning this cell
      CellLocHTmp(Cell, 1) = CellLocalTmp[Cell]; // local address on task
   }
   CellIDH  = CellIDHTmp;
   CellLocH = CellLocHTmp;

   // All done
   return;

} // end function partCellsParMetisKWay

//------------------------------------------------------------------------------
// Partition the edges based on the cell decomposition. The first cell ID in
// the CellsOnEdge array for a given edge is assigned ownership of the edge.

void Decomp::partEdges(
    const MachEnv *InEnv, ///< [in] input machine environment with MPI info
    const std::vector<I4> &CellsOnEdgeInit // [in] cell nbrs for each edge
) {

   // Retrieve some info on the MPI layout
   MPI_Comm Comm = InEnv->getComm();
   I4 NumTasks   = InEnv->getNumTasks();
   I4 MyTask     = InEnv->getMyTask();

   // Calculate some quantities associated with the initial linear
   // distribution
   I4 NCellsChunk = (NCellsGlobal - 1) / NumTasks + 1;
   I4 NEdgesChunk = (NEdgesGlobal - 1) / NumTasks + 1;
   I4 NEdgesLocal = NEdgesChunk;

   // If cells/edges do not divide evenly over processors, the last processor
   // only has the remaining cells and not a full block (chunk)
   if (MyTask == NumTasks - 1) {
      I4 StartAdd = NCellsChunk * (NumTasks - 1);
      StartAdd    = NEdgesChunk * (NumTasks - 1);
      NEdgesLocal = NEdgesGlobal - StartAdd;
   }

   // For the tasks below, it is useful to keep a sorted list of various cells
   // and edges using the std::set container and related search/sort functions.
   // We create sets for local owned cells, all local edges,
   // local owned edges and local halo edges.
   bool TimerFlag = Pacer::start("partEdgesOwned", 3);
   std::set<I4> CellsOwned;
   std::set<I4> CellsAll;
   std::set<I4> EdgesAll;
   std::set<I4> EdgesOwned;
   std::set<I4> EdgesOwnedHalo1;

   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      CellsAll.insert(CellIDH(Cell));
      if (Cell < NCellsOwned)
         CellsOwned.insert(CellIDH(Cell));
   }

   // From the EdgesOnCell array, we count and create a sorted list of
   // all edges that are needed locally. Edges that surround owned
   // cells will either be owned edges or in the first halo so we track
   // those as well for later sorting.

   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      for (int Edge = 0; Edge < MaxEdges; ++Edge) {
         I4 EdgeGlob = EdgesOnCellH(Cell, Edge);
         if (validEdgeID(EdgeGlob)) {
            EdgesAll.insert(EdgeGlob);
            if (Cell < NCellsOwned)
               EdgesOwnedHalo1.insert(EdgeGlob);
         }
      }
   }
   NEdgesAll  = EdgesAll.size();
   NEdgesSize = NEdgesAll + 1;

   // To determine whether the edge is owned by this task, we first
   // determine ownership as the first valid cell in the CellsOnEdge
   // array. CellsOnEdge is only in the initial linear decomposition so
   // we compute it there and broadcast the data to other MPI tasks.

   std::vector<I4> EdgeOwnerInit(NEdgesChunk, NCellsGlobal + 1);
   for (int Edge = 0; Edge < NEdgesLocal; ++Edge) {
      for (int Cell = 0; Cell < MaxCellsOnEdge; ++Cell) {
         I4 CellGlob = CellsOnEdgeInit[Edge * MaxCellsOnEdge + Cell];
         if (validCellID(CellGlob)) {
            EdgeOwnerInit[Edge] = CellGlob;
            break;
         }
      }
   }
   TimerFlag = Pacer::stop("partEdgesOwned", 3) && TimerFlag;

   // Broadcast the edge ownership to all tasks, one chunk at a time
   // and create a sorted list of all owned edges.

   std::vector<I4> EdgeBuf(NEdgesChunk);
   for (int ITask = 0; ITask < NumTasks; ++ITask) {

      TimerFlag = Pacer::start("partEdgesOwnerBcast", 3) && TimerFlag;
      // if it is this task's turn, fill the buffer with the owner info
      if (ITask == MyTask) {
         for (int Edge = 0; Edge < NEdgesChunk; ++Edge) {
            EdgeBuf[Edge] = EdgeOwnerInit[Edge];
         }
      }
      // Broadcast this buffer
      Broadcast(EdgeBuf, InEnv, ITask);
      TimerFlag = Pacer::stop("partEdgesOwnerBcast", 3) && TimerFlag;

      // For each edge in the buffer, check to see if the task owns
      // the cell. If so, add the edge ID to the owned edges list.
      TimerFlag = Pacer::start("partEdgesOwnerSearch", 3) && TimerFlag;
      for (int Edge = 0; Edge < NEdgesChunk; ++Edge) {

         I4 EdgeGlob =
             ITask * NEdgesChunk + Edge + 1; // Global ID in initial distrb
         // when number of edges doesn't divide evenly, the last
         // task has some invalid edges at the end so break out of loop
         if (EdgeGlob > NEdgesGlobal)
            break;

         I4 CellOwner = EdgeBuf[Edge];
         if (CellsOwned.find(CellOwner) != CellsOwned.end()) {
            // this task owns cell and therefore the edge
            EdgesOwned.insert(EdgeGlob);
         }
      }
      TimerFlag = Pacer::stop("partEdgesOwnerSearch", 3) && TimerFlag;

   } // end task loop

   // For compatibility with the previous MPAS model, we sort the
   // edges based on the order encounted in EdgesOnCell. The first halo
   // level is actually stored in reverse order from the end inward. We sort
   // edge IDs, locations and CellsOnEdge with this ordering.

   TimerFlag   = Pacer::start("partEdgesHalo", 3) && TimerFlag;
   NEdgesOwned = EdgesOwned.size();
   HostArray1DI4 NEdgesHaloTmp("NEdgesHalo", HaloWidth);
   I4 HaloCount     = EdgesOwnedHalo1.size();
   NEdgesHaloTmp(0) = HaloCount;

   HostArray1DI4 EdgeIDTmp("EdgeID", NEdgesSize);
   deepCopy(EdgeIDTmp, NEdgesGlobal + 1);

   // The owned and first halo of edges comes from the edges around
   // the owned cells, so start with these.
   I4 EdgeCount = 0;
   HaloCount--; // initialize for first halo entry in reverse order
   for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
      for (int CellEdge = 0; CellEdge < MaxEdges; ++CellEdge) {
         I4 EdgeGlob = EdgesOnCellH(Cell, CellEdge);
         // if this is a valid edge and the edge has not already been
         // processed, add the edge to the edge list and remove from
         // the relevant search arrays
         if (validEdgeID(EdgeGlob) &&
             EdgesAll.find(EdgeGlob) != EdgesAll.end()) { // edge in list
            // determine whether location is in owned or halo list
            if (EdgesOwned.find(EdgeGlob) != EdgesOwned.end()) {
               EdgeIDTmp(EdgeCount) = EdgeGlob;
               ++EdgeCount;
               EdgesAll.erase(EdgeGlob); // remove from search list
            } else {
               EdgeIDTmp(HaloCount) = EdgeGlob;
               --HaloCount;
               EdgesAll.erase(EdgeGlob); // remove from search list
            }
         }
      }
   }

   // Fill in the remaining halo regions
   I4 CellStart = NCellsOwned;
   HaloCount    = NEdgesHaloTmp(0); // reset to end of halo 1

   // The edge halos contain edges needed for the n-1 cell halo layer
   for (int Halo = 0; Halo < HaloWidth; ++Halo) {
      I4 CellEnd = NCellsHaloH(Halo);
      for (int Cell = CellStart; Cell < CellEnd; ++Cell) {
         for (int CellEdge = 0; CellEdge < MaxEdges; ++CellEdge) {
            I4 EdgeGlob = EdgesOnCellH(Cell, CellEdge);
            // if this is a valid edge and the edge has not already been
            // processed, add the edge to the halo and remove from
            // the relevant search arrays
            if (validEdgeID(EdgeGlob) &&
                EdgesAll.find(EdgeGlob) != EdgesAll.end()) { // edge in list
               EdgeIDTmp(HaloCount) = EdgeGlob;
               ++HaloCount;
               EdgesAll.erase(EdgeGlob);
            } // end if valid edge
         } // end loop over cell edges
      } // end cell loop
      // reset address range for next halo and set NEdgesHalo
      CellStart = CellEnd;
      if ((Halo + 1) < HaloWidth)
         NEdgesHaloTmp(Halo + 1) = HaloCount;
   } // end halo loop
   TimerFlag = Pacer::stop("partEdgesHalo", 3) && TimerFlag;

   // Now that we have the local lists, update the final location
   // (task, local edge address) of each of the local edges. This
   // requires one more round of communication.
   // Resize the buffer to make sure we have enough room - the distribution
   // may be less even than the original chunk size.

   TimerFlag = Pacer::start("partEdgesFinalLoc", 3) && TimerFlag;
   HostArray2DI4 EdgeLocTmp("EdgeLoc", NEdgesSize, 2);
   EdgeBuf.resize(2 * NEdgesChunk);

   // For local owned cells, the location is obvious
   // For others, we initialize to NEdgesAll and perform a search below
   std::vector<bool> EdgeFound(NEdgesSize, false);
   std::vector<I4> RemoteID;
   for (int Edge = 0; Edge < NEdgesOwned; ++Edge) {
      EdgeLocTmp(Edge, 0) = MyTask;
      EdgeLocTmp(Edge, 1) = Edge;
      EdgeFound[Edge]     = true;
   }
   for (int Edge = NEdgesOwned; Edge < NEdgesSize; ++Edge) {
      EdgeLocTmp(Edge, 0) = MyTask;
      EdgeLocTmp(Edge, 1) = NEdgesAll;
   }

   // Determine remote locations by having each task broadcast its list
   // of owned edges. Then each task searches the list for the edges it
   // needs and stores the remote address
   for (int ITask = 0; ITask < NumTasks; ++ITask) {

      TimerFlag = Pacer::start("partEdgesFinalBcast", 3) && TimerFlag;
      // fill broadcast buffer with the list of owned edges. The
      // first entry in the vector is the number of edges owned by
      // this task.
      if (ITask == MyTask) {
         EdgeBuf[0] = NEdgesOwned;
         for (int BufEdge = 0; BufEdge < NEdgesOwned; ++BufEdge) {
            EdgeBuf[BufEdge + 1] = EdgeIDTmp(BufEdge);
         }
      }
      // Broadcast the list of edges owned by this task
      Broadcast(EdgeBuf, InEnv, ITask);
      TimerFlag = Pacer::stop("partEdgesFinalBcast", 3) && TimerFlag;

      // Extract the buffer into a local search vector
      TimerFlag   = Pacer::start("partEdgesFinalSearch", 3) && TimerFlag;
      I4 BufOwned = EdgeBuf[0];
      RemoteID.resize(BufOwned);
      for (int Edge = 0; Edge < BufOwned; ++Edge)
         RemoteID[Edge] = EdgeBuf[Edge + 1];

      // For each halo point that hasn't yet been found, search the
      // current vector and store the remote address if it's found
      for (int Edge = NEdgesOwned; Edge < NEdgesAll; ++Edge) {
         if (!EdgeFound[Edge]) {
            I4 GlobID = EdgeIDTmp(Edge);
            I4 BufLoc = srchVector(RemoteID, GlobID);
            if (BufLoc < BufOwned) {
               EdgeLocTmp(Edge, 0) = ITask;  // Task that owns edge
               EdgeLocTmp(Edge, 1) = BufLoc; // Local address on task
               EdgeFound[Edge]     = true;   // Mark as found
            }
         }
      }
      TimerFlag = Pacer::stop("partEdgesFinalSearch", 3) && TimerFlag;
   }
   TimerFlag = Pacer::stop("partEdgesFinalLoc", 3) && TimerFlag;

   // Copy ID and location arrays into permanent storage
   EdgeIDH     = EdgeIDTmp;
   EdgeLocH    = EdgeLocTmp;
   NEdgesHaloH = NEdgesHaloTmp;

   if (!TimerFlag)
      LOG_WARN("Decomp::partEdges: Error encountered in timers");
   return;

} // end function partEdges

//------------------------------------------------------------------------------
// Partition the vertices based on the cell decomposition. The first cell ID
// in the CellsOnVertex array for a given vertex is assigned ownership of the
// vertex.

void Decomp::partVertices(
    const MachEnv *InEnv, ///< [in] input machine environment with MPI info
    const std::vector<I4> &CellsOnVertexInit // [in] cell nbrs for each vrtx
) {

   // Retrieve some info on the MPI layout
   MPI_Comm Comm = InEnv->getComm();
   I4 NumTasks   = InEnv->getNumTasks();
   I4 MyTask     = InEnv->getMyTask();

   // Calculate some quantities associated with the initial linear
   // distribution
   I4 NCellsChunk    = (NCellsGlobal - 1) / NumTasks + 1;
   I4 NVerticesChunk = (NVerticesGlobal - 1) / NumTasks + 1;
   I4 NVerticesLocal = NVerticesChunk;

   // If cells/edges do not divide evenly over processors, the last processor
   // only has the remaining cells and not a full block (chunk)
   if (MyTask == NumTasks - 1) {
      I4 StartAdd    = NCellsChunk * (NumTasks - 1);
      StartAdd       = NVerticesChunk * (NumTasks - 1);
      NVerticesLocal = NVerticesGlobal - StartAdd;
   }

   // For the tasks below, it is useful to keep a sorted list of various cells
   // and vertices using the std::set container and related search/sort
   // functions. We create sets for local owned cells, all local vertices,
   // local owned vertices and local halo vertices.
   bool TimerFlag = Pacer::start("partVerticesOwned", 3);
   std::set<I4> CellsOwned;
   std::set<I4> CellsAll;
   std::set<I4> VerticesAll;
   std::set<I4> VerticesOwned;
   std::set<I4> VerticesOwnedHalo1;

   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      CellsAll.insert(CellIDH(Cell));
      if (Cell < NCellsOwned)
         CellsOwned.insert(CellIDH(Cell));
   }

   // From the VerticesOnCell array, we count and create a sorted list of
   // all vertices that are needed locally. Vertices that surround owned
   // cells will either be owned edges or in the first halo so we track
   // those as well for later sorting.

   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      for (int Vrtx = 0; Vrtx < MaxEdges; ++Vrtx) {
         I4 VrtxGlob = VerticesOnCellH(Cell, Vrtx);
         if (validVertexID(VrtxGlob)) {
            VerticesAll.insert(VrtxGlob);
            if (Cell < NCellsOwned)
               VerticesOwnedHalo1.insert(VrtxGlob);
         }
      }
   }
   NVerticesAll  = VerticesAll.size();
   NVerticesSize = NVerticesAll + 1;

   // To determine whether the vertex is owned by this task, we first
   // determine ownership as the first valid cell in the CellsOnVertex
   // array. CellsOnVertex is only in the initial linear decomposition so
   // we compute it there and broadcast the data to other MPI tasks.

   std::vector<I4> VrtxOwnerInit(NVerticesChunk, NCellsGlobal + 1);
   for (int Vrtx = 0; Vrtx < NVerticesLocal; ++Vrtx) {
      for (int Cell = 0; Cell < VertexDegree; ++Cell) {
         I4 CellGlob = CellsOnVertexInit[Vrtx * VertexDegree + Cell];
         if (validCellID(CellGlob)) {
            VrtxOwnerInit[Vrtx] = CellGlob;
            break;
         }
      }
   }
   TimerFlag = Pacer::stop("partVerticesOwned", 3) && TimerFlag;

   // Broadcast the vertex ownership to all tasks, one chunk at a time
   // and create a sorted list of all owned vertices.

   std::vector<I4> VrtxBuf(NVerticesChunk);
   for (int ITask = 0; ITask < NumTasks; ++ITask) {

      TimerFlag = Pacer::start("partVerticesOwnedBcast", 3) && TimerFlag;
      // if it is this task's turn, fill the buffer with the owner info
      if (ITask == MyTask) {
         for (int Vrtx = 0; Vrtx < NVerticesChunk; ++Vrtx) {
            VrtxBuf[Vrtx] = VrtxOwnerInit[Vrtx];
         }
      }
      // Broadcast this buffer
      Broadcast(VrtxBuf, InEnv, ITask);
      TimerFlag = Pacer::stop("partVerticesOwnedBcast", 3) && TimerFlag;

      // For each vertex in the buffer, check to see if the task owns
      // the cell. If so, add the vertex ID to the owned vertices list.
      TimerFlag = Pacer::start("partVerticesOwnedSearch", 3) && TimerFlag;
      for (int Vrtx = 0; Vrtx < NVerticesChunk; ++Vrtx) {

         I4 VrtxGlob =
             ITask * NVerticesChunk + Vrtx + 1; // Global ID init distrb
         // when number of vertices doesn't divide evenly, the last
         // task has some invalid vertices at the end so break out of loop
         if (VrtxGlob > NVerticesGlobal)
            break;

         I4 CellOwner = VrtxBuf[Vrtx];
         if (CellsOwned.find(CellOwner) != CellsOwned.end()) {
            // this task owns cell and therefore the edge
            VerticesOwned.insert(VrtxGlob);
         }
      }
      TimerFlag = Pacer::stop("partVerticesOwnedSearch", 3) && TimerFlag;

   } // end task loop

   // For compatibility with the previous MPAS model, we sort the
   // vertices based on the order encounted in VerticesOnCell. The first halo
   // level is actually stored in reverse order from the end inward. We sort
   // vertex IDs, locations and CellsOnVertex with this ordering.

   TimerFlag      = Pacer::start("partVerticesHalo", 3) && TimerFlag;
   NVerticesOwned = VerticesOwned.size();
   HostArray1DI4 NVerticesHaloTmp("NVerticesHalo", HaloWidth);
   I4 HaloCount        = VerticesOwnedHalo1.size();
   NVerticesHaloTmp(0) = HaloCount;

   HostArray1DI4 VertexIDTmp("VertexID", NVerticesSize);
   deepCopy(VertexIDTmp, NVerticesGlobal + 1);

   // The owned and first halo of vertices comes from the vertices around
   // the owned cells, so start with these.
   I4 VrtxCount = 0;
   HaloCount--; // initialize for first halo entry in reverse order
   for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
      for (int CellVrtx = 0; CellVrtx < MaxEdges; ++CellVrtx) {
         I4 VrtxGlob = VerticesOnCellH(Cell, CellVrtx);
         // if this is a valid vertex and the vertex has not already been
         // processed, add the vertex to the vertex list and remove from
         // the relevant search arrays
         if (validVertexID(VrtxGlob) &&
             VerticesAll.find(VrtxGlob) !=
                 VerticesAll.end()) { // vertex in list
            // determine whether location is in owned or halo list
            if (VerticesOwned.find(VrtxGlob) != VerticesOwned.end()) {
               VertexIDTmp(VrtxCount) = VrtxGlob;
               ++VrtxCount;
               VerticesAll.erase(VrtxGlob); // remove from search list
            } else {
               VertexIDTmp(HaloCount) = VrtxGlob;
               --HaloCount;
               VerticesAll.erase(VrtxGlob); // remove from search list
            }
         }
      }
   }

   // Fill in the remaining halo regions
   I4 CellStart = NCellsOwned;
   HaloCount    = NVerticesHaloTmp(0); // reset to end of halo 1

   // The vertex halos contain vertices needed for the n-1 cell halo layer
   for (int Halo = 0; Halo < HaloWidth; ++Halo) {
      I4 CellEnd = NCellsHaloH(Halo);
      for (int Cell = CellStart; Cell < CellEnd; ++Cell) {
         for (int CellVrtx = 0; CellVrtx < MaxEdges; ++CellVrtx) {
            I4 VrtxGlob = VerticesOnCellH(Cell, CellVrtx);
            // if this is a valid vertex and the vertex has not already been
            // processed, add the vertex to the halo and remove from
            // the relevant search arrays
            if (validVertexID(VrtxGlob) &&
                VerticesAll.find(VrtxGlob) != VerticesAll.end()) { // in list
               VertexIDTmp(HaloCount) = VrtxGlob;
               ++HaloCount;
               VerticesAll.erase(VrtxGlob);
            }
         }
      }
      // reset address range for next halo and set NVerticesHalo
      CellStart = CellEnd;
      if ((Halo + 1) < HaloWidth)
         NVerticesHaloTmp(Halo + 1) = HaloCount;
   } // end halo loop
   TimerFlag = Pacer::stop("partVerticesHalo", 3) && TimerFlag;

   // Now that we have the local lists, update the final location
   // (task, local edge address) of each of the local vertices. This
   // requires one more round of communication.
   // Resize the buffer to make sure we have enough room - the distribution
   // may be less even than the original chunk size.

   TimerFlag = Pacer::start("partVerticesFinalLoc", 3) && TimerFlag;
   HostArray2DI4 VertexLocTmp("VertexLoc", NVerticesSize, 2);
   VrtxBuf.resize(2 * NVerticesChunk);

   // For local owned vertices, the location is obvious
   // For halo vertices, we initialize to NVerticesAll and perform a search
   std::vector<bool> VertexFound(NVerticesSize, false);
   std::vector<I4> RemoteID;
   for (int Vrtx = 0; Vrtx < NVerticesOwned; ++Vrtx) {
      VertexLocTmp(Vrtx, 0) = MyTask;
      VertexLocTmp(Vrtx, 1) = Vrtx;
      VertexFound[Vrtx]     = true;
   }
   for (int Vrtx = NVerticesOwned; Vrtx < NVerticesSize; ++Vrtx) {
      VertexLocTmp(Vrtx, 0) = MyTask;
      VertexLocTmp(Vrtx, 1) = NVerticesAll;
   }

   // Determine remote locations by having each task broadcast its list
   // of owned edges. Then each task searches the list for the edges it
   // needs and stores the remote address
   for (int ITask = 0; ITask < NumTasks; ++ITask) {

      // fill broadcast buffer with the list of owned vertices. The
      // first entry in the vector is the number of vertices owned by
      // this task.
      TimerFlag = Pacer::start("partVerticesFinalBcast", 3) && TimerFlag;
      if (ITask == MyTask) {
         VrtxBuf[0] = NVerticesOwned;
         for (int BufVrtx = 0; BufVrtx < NVerticesOwned; ++BufVrtx) {
            VrtxBuf[BufVrtx + 1] = VertexIDTmp(BufVrtx);
         }
      }
      // Broadcast the list of vertices owned by this task
      Broadcast(VrtxBuf, InEnv, ITask);
      TimerFlag = Pacer::stop("partVerticesFinalBcast", 3) && TimerFlag;

      // Extract the buffer into a local search vector
      TimerFlag   = Pacer::start("partVerticesFinalSearch", 3) && TimerFlag;
      I4 BufOwned = VrtxBuf[0];
      RemoteID.resize(BufOwned);
      for (int Vrtx = 0; Vrtx < BufOwned; ++Vrtx)
         RemoteID[Vrtx] = VrtxBuf[Vrtx + 1];

      // For each halo point that hasn't yet been found, search the
      // current vector and store the remote address if it's found
      for (int Vrtx = NVerticesOwned; Vrtx < NVerticesAll; ++Vrtx) {
         if (!VertexFound[Vrtx]) {
            I4 GlobID = VertexIDTmp(Vrtx);
            I4 BufLoc = srchVector(RemoteID, GlobID);
            if (BufLoc < BufOwned) {
               VertexLocTmp(Vrtx, 0) = ITask;  // Task that owns edge
               VertexLocTmp(Vrtx, 1) = BufLoc; // Local address on task
               VertexFound[Vrtx]     = true;   // Mark as found
            }
         }
      }
      TimerFlag = Pacer::stop("partVerticesFinalSearch", 3) && TimerFlag;
   }
   TimerFlag = Pacer::stop("partVerticesFinalLoc", 3) && TimerFlag;

   // Copy ID and location arrays into permanent storage
   VertexIDH      = VertexIDTmp;
   VertexLocH     = VertexLocTmp;
   NVerticesHaloH = NVerticesHaloTmp;

   if (!TimerFlag)
      LOG_WARN("Decomp::partVertices: error in timers");
   return;

} // end function partVertices

//------------------------------------------------------------------------------
// Redistribute the various XxOnCell index arrays to the final cell
// decomposition. The inputs are the various XxOnCell arrays in the
// initial linear distribution. On exit, all the XxOnCell arrays are
// in the correct final domain decomposition.

void Decomp::rearrangeCellArrays(
    const MachEnv *InEnv, // input machine environment for MPI layout
    const std::vector<I4> &CellsOnCellInit,   //< [in] cell nbrs on each edge
    const std::vector<I4> &EdgesOnCellInit,   //< [in] edges around each cell
    const std::vector<I4> &VerticesOnCellInit //< [in] vertices around cell
) {

   bool TimerFlag = Pacer::start("rearrangeCellArrays", 3);

   // Extract some MPI information
   MPI_Comm Comm = InEnv->getComm();
   I4 NumTasks   = InEnv->getNumTasks();
   I4 MyTask     = InEnv->getMyTask();

   // Define the chunk sizes for the initial linear distribution
   I4 NCellsChunk = (NCellsGlobal - 1) / NumTasks + 1;

   // Create a buffer for sending all the cell information
   // Each of the 3 arrays are (NCellsChunk,MaxEdges)
   I4 SizePerCell = MaxEdges * 3;
   I4 BufSize     = NCellsChunk * SizePerCell;
   std::vector<I4> CellBuf(BufSize);

   // Create temporary arrays for holding the XxOnCell results
   // and initialize to NXxGlobal+1 to denote a non-existent entry
   HostArray2DI4 CellsOnCellTmp("CellsOnCell", NCellsSize, MaxEdges);
   HostArray2DI4 EdgesOnCellTmp("EdgesOnCell", NCellsSize, MaxEdges);
   HostArray2DI4 VerticesOnCellTmp("VerticesOnCell", NCellsSize, MaxEdges);
   HostArray1DI4 NEdgesOnCellTmp("NEdgesOnCell", NCellsSize);
   deepCopy(CellsOnCellTmp, NCellsGlobal + 1);
   deepCopy(EdgesOnCellTmp, NEdgesGlobal + 1);
   deepCopy(VerticesOnCellTmp, NVerticesGlobal + 1);
   deepCopy(NEdgesOnCellTmp, 0);

   // For each cell this task owns, determine the task and local address
   // of the cell in the initial linear distribution
   std::vector<I4> TaskInit(NCellsSize, NumTasks);
   std::vector<I4> AddInit(NCellsSize, NCellsGlobal + 1);
   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      int GlobalID = CellIDH(Cell); // one-based
      if (validCellID(GlobalID)) {
         TaskInit[Cell] = (GlobalID - 1) / NCellsChunk;
         AddInit[Cell]  = GlobalID - 1 - TaskInit[Cell] * NCellsChunk;
      }
   }

   // Each task will broadcast the cells it owns in the initial linear
   // distribution and all tasks will search their list of cells to determine
   // entries it needs.
   for (int ITask = 0; ITask < NumTasks; ++ITask) {

      TimerFlag = Pacer::start("rearrangeCellsBcast", 3) && TimerFlag;
      // If it is this task's turn to send, fill the buffer with the local
      // chunk of all three arrays.
      if (MyTask == ITask) { // Fill buffer with local chunk
         for (int Cell = 0; Cell < NCellsChunk; ++Cell) {
            for (int Edge = 0; Edge < MaxEdges; ++Edge) {
               I4 BufAdd           = Cell * SizePerCell + Edge * 3;
               I4 ArrayAdd         = Cell * MaxEdges + Edge;
               CellBuf[BufAdd]     = CellsOnCellInit[ArrayAdd];
               CellBuf[BufAdd + 1] = VerticesOnCellInit[ArrayAdd];
               CellBuf[BufAdd + 2] = EdgesOnCellInit[ArrayAdd];
            }
         }
      }
      Broadcast(CellBuf, InEnv, ITask);
      TimerFlag = Pacer::stop("rearrangeCellsBcast", 3) && TimerFlag;

      // For each cell needed locally, we can compute the task and address
      // of the cell in the linear distribution and extract from the buffer
      // if needed

      TimerFlag = Pacer::start("rearrangeCellsSearch", 3) && TimerFlag;
      for (int Cell = 0; Cell < NCellsAll; ++Cell) {
         if (TaskInit[Cell] == ITask) {
            if (AddInit[Cell] < NCellsGlobal + 1) {

               // Local cell needs the info so extract from the buffer
               // into the local address. For edges, we only store the
               // active edges and maintain a count of the edges.
               I4 EdgeCount = 0;
               for (int Edge = 0; Edge < MaxEdges; ++Edge) {
                  I4 BufAdd  = AddInit[Cell] * SizePerCell + Edge * 3;
                  I4 NbrCell = CellBuf[BufAdd];
                  I4 NbrVrtx = CellBuf[BufAdd + 1];
                  I4 NbrEdge = CellBuf[BufAdd + 2];
                  if (validCellID(NbrCell)) {
                     CellsOnCellTmp(Cell, Edge) = NbrCell;
                  } else {
                     CellsOnCellTmp(Cell, Edge) = NCellsGlobal + 1;
                  }
                  if (validVertexID(NbrVrtx)) {
                     VerticesOnCellTmp(Cell, Edge) = NbrVrtx;
                  } else {
                     VerticesOnCellTmp(Cell, Edge) = NVerticesGlobal + 1;
                  }
                  if (validEdgeID(NbrEdge)) {
                     EdgesOnCellTmp(Cell, EdgeCount) = NbrEdge;
                     EdgeCount++;
                  }
               }
               NEdgesOnCellTmp(Cell) = EdgeCount;
            }
         }
      }
      TimerFlag = Pacer::stop("rearrangeCellsSearch", 3) && TimerFlag;
   } // end loop over MPI tasks

   // Copy to final location on host - wait to create device copies until
   // the entries are translated to local addresses rather than global IDs
   CellsOnCellH    = CellsOnCellTmp;
   EdgesOnCellH    = EdgesOnCellTmp;
   VerticesOnCellH = VerticesOnCellTmp;
   NEdgesOnCellH   = NEdgesOnCellTmp;

   // All done
   TimerFlag = Pacer::stop("rearrangeCellArrays", 3) && TimerFlag;
   return;

} // end function rearrangeCellArrays

//------------------------------------------------------------------------------
// Redistribute the various XxOnEdge index arrays to the final edge
// decomposition. The inputs are the various XxOnEdge arrays in the
// initial linear distribution. On exit, all the XxOnEdge arrays are
// in the correct final domain decomposition.

void Decomp::rearrangeEdgeArrays(
    const MachEnv *InEnv, // input machine environment for MPI layout
    const std::vector<I4> &CellsOnEdgeInit,   //< [in] cell nbrs on each edge
    const std::vector<I4> &EdgesOnEdgeInit,   //< [in] edges around nbr cells
    const std::vector<I4> &VerticesOnEdgeInit //< [in] vertices at edge end
) {

   bool TimerFlag = true; // timer return value

   // Extract some MPI information
   MPI_Comm Comm = InEnv->getComm();
   I4 NumTasks   = InEnv->getNumTasks();
   I4 MyTask     = InEnv->getMyTask();

   // Define the chunk sizes for the initial linear distribution
   I4 NEdgesChunk = (NEdgesGlobal - 1) / NumTasks + 1;

   // Create a buffer for sending all the edge array information
   // Each of the arrays are either NEdgesChunk*2 (cell, vertex)
   // or NEdgesChunk*MaxEdges*2 (edge)
   I4 SizePerEdge = 2 * MaxEdges + MaxCellsOnEdge + 2;
   I4 BufSize     = NEdgesChunk * SizePerEdge;
   std::vector<I4> EdgeBuf(BufSize);

   // Create temporary arrays for holding the XxOnEdge results
   // and initialize to NXxGlobal+1 to denote a non-existent entry
   HostArray2DI4 CellsOnEdgeTmp("CellsOnEdge", NEdgesSize, MaxCellsOnEdge);
   HostArray2DI4 EdgesOnEdgeTmp("EdgesOnEdge", NEdgesSize, 2 * MaxEdges);
   HostArray2DI4 VerticesOnEdgeTmp("VerticesOnEdge", NEdgesSize, 2);
   HostArray1DI4 NEdgesOnEdgeTmp("NEdgesOnEdge", NEdgesSize);
   deepCopy(CellsOnEdgeTmp, NCellsGlobal + 1);
   deepCopy(EdgesOnEdgeTmp, NEdgesGlobal + 1);
   deepCopy(VerticesOnEdgeTmp, NVerticesGlobal + 1);
   deepCopy(NEdgesOnEdgeTmp, 0);

   // Compute task and local address for needed edges in initial linear
   // distribution. Initialize to invalid task, addresses.
   std::vector<I4> TaskInit(NEdgesAll, NumTasks);
   std::vector<I4> AddInit(NEdgesAll, NEdgesGlobal + 1);
   for (int Edge = 0; Edge < NEdgesAll; ++Edge) {
      I4 GlobalID = EdgeIDH(Edge); // one-based ID
      if (validEdgeID(GlobalID)) {
         TaskInit[Edge] = (GlobalID - 1) / NEdgesChunk;
         AddInit[Edge]  = GlobalID - 1 - TaskInit[Edge] * NEdgesChunk;
      }
   }

   // Each task will broadcast the array chunks it owns in the initial linear
   // distribution and all tasks will search that list and extract the
   // entries it owns.
   for (int ITask = 0; ITask < NumTasks; ++ITask) {

      TimerFlag = Pacer::start("rearrangeEdgeArraysBcast", 3) && TimerFlag;
      // If it is this task's turn to send, fill the buffer with the local
      // chunk of all three arrays.
      if (MyTask == ITask) { // Fill buffer with local chunk
         for (int Edge = 0; Edge < NEdgesChunk; ++Edge) {
            I4 BufAdd = Edge * SizePerEdge;
            for (int Cell = 0; Cell < MaxCellsOnEdge; ++Cell) {
               I4 ArrayAdd     = Edge * MaxCellsOnEdge + Cell;
               EdgeBuf[BufAdd] = CellsOnEdgeInit[ArrayAdd];
               ++BufAdd;
            }
            for (int Vrtx = 0; Vrtx < 2; ++Vrtx) {
               I4 ArrayAdd     = Edge * 2 + Vrtx;
               EdgeBuf[BufAdd] = VerticesOnEdgeInit[ArrayAdd];
               ++BufAdd;
            }
            for (int NbrEdge = 0; NbrEdge < 2 * MaxEdges; ++NbrEdge) {
               I4 ArrayAdd     = Edge * 2 * MaxEdges + NbrEdge;
               EdgeBuf[BufAdd] = EdgesOnEdgeInit[ArrayAdd];
               ++BufAdd;
            }
         }
      }
      Broadcast(EdgeBuf, InEnv, ITask);
      TimerFlag = Pacer::stop("rearrangeEdgeArraysBcast", 3) && TimerFlag;

      // If the local Edge array has points in this buffer, extract the
      // array information into the proper location
      TimerFlag = Pacer::start("rearrangeEdgeArraysSearch", 3) && TimerFlag;
      for (int Edge = 0; Edge < NEdgesAll; ++Edge) {
         if (TaskInit[Edge] == ITask) {
            if (AddInit[Edge] < NEdgesGlobal + 1) {

               // Local edge exists in the buffer, compute starting address
               // and extract CellsOnEdge, VerticesOnEdge
               I4 BufAdd = AddInit[Edge] * SizePerEdge;
               for (int Cell = 0; Cell < MaxCellsOnEdge; ++Cell) {
                  CellsOnEdgeTmp(Edge, Cell) = EdgeBuf[BufAdd];
                  ++BufAdd;
               }
               for (int Vrtx = 0; Vrtx < 2; ++Vrtx) {
                  VerticesOnEdgeTmp(Edge, Vrtx) = EdgeBuf[BufAdd];
                  ++BufAdd;
               }
               // In the EdgeOnEdge array, a zero entry must be kept in
               // place but assigned the boundary value NEdgesGlobal+1
               I4 EdgeCount = 0;
               for (int NbrEdge = 0; NbrEdge < 2 * MaxEdges; ++NbrEdge) {
                  I4 EdgeID = EdgeBuf[BufAdd];
                  ++BufAdd;
                  if (EdgeID == 0) {
                     EdgesOnEdgeTmp(Edge, EdgeCount) = NEdgesGlobal + 1;
                     EdgeCount++;
                  } else if (validEdgeID(EdgeID)) {
                     EdgesOnEdgeTmp(Edge, EdgeCount) = EdgeID;
                     EdgeCount++;
                  }
               }
               NEdgesOnEdgeTmp(Edge) = EdgeCount;
            } // end if address in buffer
         } // end if address on this task
      } // end loop over local edges
      TimerFlag = Pacer::stop("rearrangeEdgeArraysSearch", 3) && TimerFlag;
   } // end loop over MPI tasks

   // Copy to final location on host - wait to create device copies until
   // the entries are translated to local addresses rather than global IDs
   CellsOnEdgeH    = CellsOnEdgeTmp;
   EdgesOnEdgeH    = EdgesOnEdgeTmp;
   VerticesOnEdgeH = VerticesOnEdgeTmp;
   NEdgesOnEdgeH   = NEdgesOnEdgeTmp;

   // All done
   if (!TimerFlag)
      LOG_WARN("Decomp::rearrangeEdgeArrays: Error in timers");
   return;

} // end function rearrangeEdgeArrays

//------------------------------------------------------------------------------
// Redistribute the various XxOnVertex index arrays to the final vertex
// decomposition. The inputs are the various XxOnVertex arrays in the
// initial linear distribution. On exit, all the XxOnVertex arrays are
// in the correct final domain decomposition.

void Decomp::rearrangeVertexArrays(
    const MachEnv *InEnv, // input machine environment for MPI layout
    const std::vector<I4> &CellsOnVertexInit, //< [in] cells at each vrtx
    const std::vector<I4> &EdgesOnVertexInit  //< [in] edges joined at vrtx
) {

   bool TimerFlag = true; // timer return code

   // Extract some MPI information
   MPI_Comm Comm = InEnv->getComm();
   I4 NumTasks   = InEnv->getNumTasks();
   I4 MyTask     = InEnv->getMyTask();

   // Define the chunk sizes for the initial linear distribution
   I4 NVerticesChunk = (NVerticesGlobal - 1) / NumTasks + 1;

   // Create a buffer for sending all the vertex array information
   // Both of the arrays should be of size NVerticesChunk*VertexDegree
   // so the full buffer is twice that.
   I4 SizePerVrtx = 2 * VertexDegree;
   I4 BufSize     = NVerticesChunk * SizePerVrtx;
   std::vector<I4> VrtxBuf(BufSize);

   // Create temporary arrays for holding the XxOnVertex results
   // and initialize to NXxGlobal+1 to denote a non-existent entry
   HostArray2DI4 CellsOnVertexTmp("CellsOnVertex", NVerticesSize, VertexDegree);
   HostArray2DI4 EdgesOnVertexTmp("EdgesOnVertex", NVerticesSize, VertexDegree);
   deepCopy(CellsOnVertexTmp, NCellsGlobal + 1);
   deepCopy(EdgesOnVertexTmp, NEdgesGlobal + 1);

   // Compute the task and local address of each needed local vertex in the
   // original linear decomposition
   std::vector<I4> TaskInit(NVerticesAll, NumTasks);
   std::vector<I4> AddInit(NVerticesAll, NVerticesGlobal + 1);
   for (int Vrtx = 0; Vrtx < NVerticesAll; ++Vrtx) {
      I4 GlobalID = VertexIDH(Vrtx); // one-based ID
      if (validVertexID(GlobalID)) {
         TaskInit[Vrtx] = (GlobalID - 1) / NVerticesChunk;
         AddInit[Vrtx]  = GlobalID - 1 - TaskInit[Vrtx] * NVerticesChunk;
      }
   }

   // Each task will broadcast the array chunks it owns in the initial linear
   // distribution and all tasks will search that list and extract the
   // entries it owns.
   for (int ITask = 0; ITask < NumTasks; ++ITask) {

      // If it is this task's turn to send, fill the buffer with the local
      // chunk of both arrays.
      TimerFlag = Pacer::start("rearrangeVertexArraysBcast", 3) && TimerFlag;
      if (MyTask == ITask) { // Fill buffer with local chunk
         for (int Vrtx = 0; Vrtx < NVerticesChunk; ++Vrtx) {
            I4 BufAdd = Vrtx * SizePerVrtx;
            for (int Cell = 0; Cell < VertexDegree; ++Cell) {
               I4 ArrayAdd     = Vrtx * VertexDegree + Cell;
               VrtxBuf[BufAdd] = CellsOnVertexInit[ArrayAdd];
               ++BufAdd;
            }
            for (int Edge = 0; Edge < VertexDegree; ++Edge) {
               I4 ArrayAdd     = Vrtx * VertexDegree + Edge;
               VrtxBuf[BufAdd] = EdgesOnVertexInit[ArrayAdd];
               ++BufAdd;
            }
         }
      }
      Broadcast(VrtxBuf, InEnv, ITask);
      TimerFlag = Pacer::stop("rearrangeVertexArraysBcast", 3) && TimerFlag;

      // For each local vertex in the distribution, determine whether the
      // buffer contains the vertex and extract the arrays
      TimerFlag = Pacer::start("rearrangeVertexArraysSearch", 3) && TimerFlag;
      for (int Vrtx = 0; Vrtx < NVerticesAll; ++Vrtx) {
         if (TaskInit[Vrtx] == ITask) {
            if (AddInit[Vrtx] < NVerticesGlobal + 1) {

               // Local vertex is in the buffer so compute a buffer starting
               // address
               I4 BufAdd = AddInit[Vrtx] * SizePerVrtx;
               for (int Cell = 0; Cell < VertexDegree; ++Cell) {
                  CellsOnVertexTmp(Vrtx, Cell) = VrtxBuf[BufAdd];
                  ++BufAdd;
               }
               for (int Edge = 0; Edge < VertexDegree; ++Edge) {
                  EdgesOnVertexTmp(Vrtx, Edge) = VrtxBuf[BufAdd];
                  ++BufAdd;
               }
            } // end if valid vertex
         } // end if task has the vertex
      } // end loop over local vertices
      TimerFlag = Pacer::stop("rearrangeVertexArraysSearch", 3) && TimerFlag;
   } // end loop over MPI tasks

   // Copy to final location on host - wait to create device copies until
   // the entries are translated to local addresses rather than global IDs
   CellsOnVertexH = CellsOnVertexTmp;
   EdgesOnVertexH = EdgesOnVertexTmp;

   // All done
   if (!TimerFlag)
      LOG_WARN("Decomp::rearrangeVertexArrays: Error in timers");
   return;

} // end function rearrangeVertexArrays

//------------------------------------------------------------------------------
// Utility routine to convert a partition method string into PartMethod enum

PartMethod getPartMethodFromStr(const std::string &InMethod) {

   // convert string to lower case for easier equivalence checking
   std::string MethodComp = InMethod;
   std::transform(MethodComp.begin(), MethodComp.end(), MethodComp.begin(),
                  [](unsigned char c) { return std::tolower(c); });

   // Check supported methods and return appropriate enum
   if (MethodComp == "metiskway") { // serial Metis KWay
      return PartMethodMetisKWay;

   } else if (MethodComp == "parmetiskway") { // parallel ParMetis KWay
      return PartMethodParMetisKWay;

   } else {
      return PartMethodUnknown;

   } // end branch on method string

} // End getPartMethodFromStr

//------------------------------------------------------------------------------
// end Decomp methods

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
