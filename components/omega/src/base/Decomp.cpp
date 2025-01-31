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
#include "Config.h"
#include "DataTypes.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
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

int readMesh(const int MeshFileID, // file ID for open mesh file
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

   int Err = 0;

   // Retrieve some info on the MPI layout
   MPI_Comm Comm = InEnv->getComm();
   I4 NumTasks   = InEnv->getNumTasks();
   I4 MyTask     = InEnv->getMyTask();
   I4 MasterTask = InEnv->getMasterTask();
   bool IsMaster = InEnv->isMasterTask();

   // Read in mesh size information - these are dimension lengths in
   // the input mesh file. Check both the name under Omega name conventions
   // and the older MPAS name.
   std::string DimName    = "NCells";
   std::string DimNameOld = "nCells";
   I4 NCellsID;
   Err = IO::getDimFromFile(MeshFileID, DimName, NCellsID, NCellsGlobal);
   if (Err != 0) { // dim not found, try again with older MPAS name
      Err = IO::getDimFromFile(MeshFileID, DimNameOld, NCellsID, NCellsGlobal);
      if (Err != 0 or NCellsGlobal <= 0)
         LOG_CRITICAL("Decomp: error reading nCells");
   }

   DimName    = "NEdges";
   DimNameOld = "nEdges";
   I4 NEdgesID;
   Err = IO::getDimFromFile(MeshFileID, DimName, NEdgesID, NEdgesGlobal);
   if (Err != 0) { // dim not found, try again with older MPAS name
      Err = IO::getDimFromFile(MeshFileID, DimNameOld, NEdgesID, NEdgesGlobal);
      if (Err != 0 or NEdgesGlobal <= 0)
         LOG_CRITICAL("Decomp: error reading NEdges");
   }

   DimName    = "NVertices";
   DimNameOld = "nVertices";
   I4 NVerticesID;
   Err = IO::getDimFromFile(MeshFileID, DimName, NVerticesID, NVerticesGlobal);
   if (Err != 0) { // dim not found, try again with older MPAS name
      Err = IO::getDimFromFile(MeshFileID, DimNameOld, NVerticesID,
                               NVerticesGlobal);
      if (Err != 0 or NVerticesGlobal <= 0)
         LOG_CRITICAL("Decomp: error reading NVertices");
   }

   DimName    = "MaxEdges";
   DimNameOld = "maxEdges";
   I4 MaxEdgesID;
   Err = IO::getDimFromFile(MeshFileID, DimName, MaxEdgesID, MaxEdges);
   if (Err != 0) { // dim not found, try again with older MPAS name
      Err = IO::getDimFromFile(MeshFileID, DimNameOld, MaxEdgesID, MaxEdges);
      if (Err != 0 or MaxEdges <= 0)
         LOG_CRITICAL("Decomp: error reading MaxEdges");
   }

   DimName    = "VertexDegree";
   DimNameOld = "vertexDegree";
   I4 VertexDegreeID;
   Err = IO::getDimFromFile(MeshFileID, DimName, VertexDegreeID, VertexDegree);
   if (Err != 0) { // dim not found, try again with older MPAS name
      Err = IO::getDimFromFile(MeshFileID, DimNameOld, VertexDegreeID,
                               VertexDegree);
      if (Err != 0 or VertexDegree <= 0)
         LOG_CRITICAL("Decomp: error reading VertexDegree");
   }

   DimName    = "MaxCellsOnEdge";
   DimNameOld = "TWO";
   I4 MaxCellsOnEdgeID;
   Err = IO::getDimFromFile(MeshFileID, DimNameOld, MaxCellsOnEdgeID,
                            MaxCellsOnEdge);
   if (Err != 0) { // dim not found, try again with older MPAS name
      Err = IO::getDimFromFile(MeshFileID, DimName, MaxCellsOnEdgeID,
                               MaxCellsOnEdge);
      if (Err != 0 or MaxCellsOnEdge <= 0)
         LOG_CRITICAL("Decomp: error reading MaxCellsOnEdge");
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
   I4 OnCellDecomp;
   I4 OnEdgeDecomp;
   I4 OnEdgeDecomp2;
   I4 OnVertexDecomp;
   Err = IO::createDecomp(OnCellDecomp, IO::IOTypeI4, NDims, OnCellDims,
                          OnCellSize, OnCellOffset, Rearr);
   if (Err != 0)
      LOG_CRITICAL("Decomp: error creating OnCell IO decomposition");
   Err = IO::createDecomp(OnEdgeDecomp, IO::IOTypeI4, NDims, OnEdgeDims,
                          OnEdgeSize, OnEdgeOffset, Rearr);
   if (Err != 0)
      LOG_CRITICAL("Decomp: error creating OnEdge IO decomposition");
   Err = IO::createDecomp(OnEdgeDecomp2, IO::IOTypeI4, NDims, OnEdgeDims2,
                          OnEdgeSize2, OnEdgeOffset2, Rearr);
   if (Err != 0)
      LOG_CRITICAL("Decomp: error creating OnEdg2 IO decomposition");
   Err = IO::createDecomp(OnVertexDecomp, IO::IOTypeI4, NDims, OnVertexDims,
                          OnVertexSize, OnVertexOffset, Rearr);
   if (Err != 0)
      LOG_CRITICAL("Decomp: error creating Vertex IO decomposition");

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
   if (Err != 0) { // not found, try again under older name
      Err = IO::readArray(&CellsOnCellInit[0], OnCellSize, VarNameOld,
                          MeshFileID, OnCellDecomp, CellsOnCellID);
      if (Err != 0)
         LOG_CRITICAL("Decomp: error reading CellsOnCell");
   }

   VarName    = "EdgesOnCell";
   VarNameOld = "edgesOnCell";
   int EdgesOnCellID;
   Err = IO::readArray(&EdgesOnCellInit[0], OnCellSize, VarName, MeshFileID,
                       OnCellDecomp, EdgesOnCellID);
   if (Err != 0) { // not found, try again under older name
      Err = IO::readArray(&EdgesOnCellInit[0], OnCellSize, VarNameOld,
                          MeshFileID, OnCellDecomp, EdgesOnCellID);
      if (Err != 0)
         LOG_CRITICAL("Decomp: error reading EdgesOnCell");
   }

   VarName    = "VerticesOnCell";
   VarNameOld = "verticesOnCell";
   int VerticesOnCellID;
   Err = IO::readArray(&VerticesOnCellInit[0], OnCellSize, VarName, MeshFileID,
                       OnCellDecomp, VerticesOnCellID);
   if (Err != 0) { // not found, try again under older name
      Err = IO::readArray(&VerticesOnCellInit[0], OnCellSize, VarNameOld,
                          MeshFileID, OnCellDecomp, VerticesOnCellID);
      if (Err != 0)
         LOG_CRITICAL("Decomp: error reading VerticesOnCell");
   }

   VarName    = "CellsOnEdge";
   VarNameOld = "cellsOnEdge";
   int CellsOnEdgeID;
   Err = IO::readArray(&CellsOnEdgeInit[0], OnEdgeSize, VarName, MeshFileID,
                       OnEdgeDecomp, CellsOnEdgeID);
   if (Err != 0) { // not found, try again under older name
      Err = IO::readArray(&CellsOnEdgeInit[0], OnEdgeSize, VarNameOld,
                          MeshFileID, OnEdgeDecomp, CellsOnEdgeID);
      if (Err != 0)
         LOG_CRITICAL("Decomp: error reading CellsOnEdge");
   }

   VarName    = "EdgesOnEdge";
   VarNameOld = "edgesOnEdge";
   int EdgesOnEdgeID;
   Err = IO::readArray(&EdgesOnEdgeInit[0], OnEdgeSize2, VarName, MeshFileID,
                       OnEdgeDecomp2, EdgesOnEdgeID);
   if (Err != 0) { // not found, try again under older name
      Err = IO::readArray(&EdgesOnEdgeInit[0], OnEdgeSize2, VarNameOld,
                          MeshFileID, OnEdgeDecomp2, EdgesOnEdgeID);
      if (Err != 0)
         LOG_CRITICAL("Decomp: error reading EdgesOnEdge");
   }

   VarName    = "VerticesOnEdge";
   VarNameOld = "verticesOnEdge";
   int VerticesOnEdgeID;
   Err = IO::readArray(&VerticesOnEdgeInit[0], OnEdgeSize, VarName, MeshFileID,
                       OnEdgeDecomp, VerticesOnEdgeID);
   if (Err != 0) { // not found, try again under older name
      Err = IO::readArray(&VerticesOnEdgeInit[0], OnEdgeSize, VarNameOld,
                          MeshFileID, OnEdgeDecomp, VerticesOnEdgeID);
      if (Err != 0)
         LOG_CRITICAL("Decomp: error reading VerticesOnEdge");
   }

   VarName    = "CellsOnVertex";
   VarNameOld = "cellsOnVertex";
   int CellsOnVertexID;
   Err = IO::readArray(&CellsOnVertexInit[0], OnVertexSize, VarName, MeshFileID,
                       OnVertexDecomp, CellsOnVertexID);
   if (Err != 0) { // not found, try again under older name
      Err = IO::readArray(&CellsOnVertexInit[0], OnVertexSize, VarNameOld,
                          MeshFileID, OnVertexDecomp, CellsOnVertexID);
      if (Err != 0)
         LOG_CRITICAL("Decomp: error reading CellsOnVertex");
   }

   VarName    = "EdgesOnVertex";
   VarNameOld = "edgesOnVertex";
   int EdgesOnVertexID;
   Err = IO::readArray(&EdgesOnVertexInit[0], OnVertexSize, VarName, MeshFileID,
                       OnVertexDecomp, EdgesOnVertexID);
   if (Err != 0) { // not found, try again under older name
      Err = IO::readArray(&EdgesOnVertexInit[0], OnVertexSize, VarNameOld,
                          MeshFileID, OnVertexDecomp, EdgesOnVertexID);
      if (Err != 0)
         LOG_CRITICAL("Decomp: error reading EdgesOnVertex");
   }

   // Initial decompositions are no longer needed so remove them now
   Err = IO::destroyDecomp(OnCellDecomp);
   if (Err != 0)
      LOG_ERROR("Decomp: error destroying OnCell decomposition");
   Err = IO::destroyDecomp(OnEdgeDecomp);
   if (Err != 0)
      LOG_ERROR("Decomp: error destroying OnEdge decomposition");
   Err = IO::destroyDecomp(OnEdgeDecomp2);
   if (Err != 0)
      LOG_ERROR("Decomp: error destroying OnEdge2 decomposition");
   Err = IO::destroyDecomp(OnVertexDecomp);
   if (Err != 0)
      LOG_ERROR("Decomp: error destroying OnVertex decomposition");

   return Err;

} // end readMesh

//------------------------------------------------------------------------------
// Initialize the decomposition and create the default decomposition with
// (currently) one partition per MPI task using a ParMetis KWay method.

int Decomp::init(const std::string &MeshFileName) {

   int Err = 0; // default successful return code

   I4 InHaloWidth;
   std::string DecompMethodStr;

   // Retrieve options from Config
   Config *OmegaConfig = Config::getOmegaConfig();

   Config DecompConfig("Decomp");
   Err = OmegaConfig->get(DecompConfig);
   if (Err != 0) {
      LOG_CRITICAL("Decomp: Decomp group not found in Config");
      return Err;
   }

   Err = DecompConfig.get("HaloWidth", InHaloWidth);
   if (Err != 0) {
      LOG_CRITICAL("Decomp: HaloWidth not found in Decomp Config");
      return Err;
   }

   Err = DecompConfig.get("DecompMethod", DecompMethodStr);
   if (Err != 0) {
      LOG_CRITICAL("Decomp: DecompMethod not found in Decomp Config");
      return Err;
   }

   PartMethod Method = getPartMethodFromStr(DecompMethodStr);

   // Retrieve the default machine environment
   MachEnv *DefEnv = MachEnv::getDefault();

   // Use one partition per MPI task as the default
   I4 NParts = DefEnv->getNumTasks();

   // Create the default decomposition and set pointer to it
   Decomp::DefaultDecomp = Decomp::create("Default", DefEnv, NParts, Method,
                                          InHaloWidth, MeshFileName);

   return Err;

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

   int Err = 0; // internal error code

   // Retrieve some info on the MPI layout
   MPI_Comm Comm = InEnv->getComm();
   I4 NumTasks   = InEnv->getNumTasks();
   I4 MyTask     = InEnv->getMyTask();
   I4 MasterTask = InEnv->getMasterTask();
   bool IsMaster = InEnv->isMasterTask();

   // Open the mesh file for reading (assume IO has already been initialized)
   int FileID;
   MeshFileName = MeshFileName_;
   Err          = IO::openFile(FileID, MeshFileName, IO::ModeRead);
   if (Err != 0)
      LOG_CRITICAL("Decomp: error opening mesh file");

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

   Err = readMesh(FileID, InEnv, NCellsGlobal, NEdgesGlobal, NVerticesGlobal,
                  MaxEdges, MaxCellsOnEdge, VertexDegree, CellsOnCellInit,
                  EdgesOnCellInit, VerticesOnCellInit, CellsOnEdgeInit,
                  EdgesOnEdgeInit, VerticesOnEdgeInit, CellsOnVertexInit,
                  EdgesOnVertexInit);
   if (Err != 0)
      LOG_CRITICAL("Decomp: Error reading mesh connectivity");

   // Close file
   Err = IO::closeFile(FileID);

   // In the case of single task avoid calling a full partitioning routine and
   // just set the needed variables directly. This is done because some METIS
   // functions can raise SIGFPE when numparts == 1 due to division by zero
   // See: https://github.com/KarypisLab/METIS/issues/67
   if (NumTasks == 1) {
      partCellsSingleTask();
   } else {
      // Use the mesh adjacency information to create a partition of cells
      switch (Method) { // branch depending on method chosen

      //---------------------------------------------------------------------------
      // ParMetis KWay method
      case PartMethodMetisKWay: {

         Err = partCellsKWay(InEnv, CellsOnCellInit);
         if (Err != 0) {
            LOG_CRITICAL("Decomp: Error partitioning cells KWay");
            return;
         }
         break;
      } // end case MethodKWay

         //---------------------------------------------------------------------------
         // Unknown partitioning method

      default:
         LOG_CRITICAL("Decomp: Unknown or unsupported decomposition method");
         return;

      } // End switch on Method
   }

   //---------------------------------------------------------------------------

   // Cell partitioning complete. Redistribute the initial XXOnCell arrays
   // to their final locations.
   Err = rearrangeCellArrays(InEnv, CellsOnCellInit, EdgesOnCellInit,
                             VerticesOnCellInit);
   if (Err != 0) {
      LOG_CRITICAL("Decomp: Error rearranging XxOnCell arrays");
      return;
   }

   // Partition the edges
   Err = partEdges(InEnv, CellsOnEdgeInit);
   if (Err != 0) {
      LOG_CRITICAL("Decomp: Error partitioning edges");
      return;
   }

   // Edge partitioning complete. Redistribute the initial XXOnEdge arrays
   // to their final locations.
   Err = rearrangeEdgeArrays(InEnv, CellsOnEdgeInit, EdgesOnEdgeInit,
                             VerticesOnEdgeInit);
   if (Err != 0) {
      LOG_CRITICAL("Decomp: Error rearranging XxOnEdge arrays");
      return;
   }

   // Partition the vertices
   Err = partVertices(InEnv, CellsOnVertexInit);
   if (Err != 0) {
      LOG_CRITICAL("Decomp: Error partitioning vertices");
      return;
   }

   // Vertex partitioning complete. Redistribute the initial XXOnVertex arrays
   // to their final locations.
   Err = rearrangeVertexArrays(InEnv, CellsOnVertexInit, EdgesOnVertexInit);
   if (Err != 0) {
      LOG_CRITICAL("Decomp: Error rearranging XxOnVertex arrays");
      return;
   }

   // Convert global addresses to local addresses. Create the global to
   // local address ordered maps to simplify and optimize searches.
   // Invalid/non-existent edges have been assigned NXxGlobal+1 and we want
   // to map that to the local NXxAll+1 (NXxSize). We insert that value
   // first in the map so that later attempts to change will be ignored.

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

   // Create device copies of all arrays

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
} // end decomposition constructor

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

   // Check to see if a decomposition of the same name already exists and
   // if so, exit with an error
   if (AllDecomps.find(Name) != AllDecomps.end()) {
      LOG_ERROR("Attempted to create a Decomp with name {} but a Decomp of "
                "that name already exists",
                Name);
      return nullptr;
   }

   // create a new decomp on the heap and put it in a map of
   // unique_ptrs, which will manage its lifetime
   auto *NewDecomp =
       new Decomp(Name, Env, NParts, Method, HaloWidth, MeshFileName);
   AllDecomps.emplace(Name, NewDecomp);

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

   AllDecomps.clear(); // removes all decomps from the list (map) and in
                       // the process, calls the destructors for each

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

   // if found, return the decomposition pointer
   if (it != AllDecomps.end()) {
      return it->second.get();

      // otherwise print an error and return a null pointer
   } else {
      LOG_ERROR("Decomp::get: Attempt to retrieve non-existent Decomposition:");
      LOG_ERROR(" {} has not been defined or has been removed", Name);
      return nullptr;
   }

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
// Partition the cells using the Metis/ParMetis KWay method
// After this partitioning, the decomposition class member CellID and
// CellLocator arrays have been set as well as the class NCells size variables
// (owned, halo, all).

int Decomp::partCellsKWay(
    const MachEnv *InEnv, // [in] input machine environment with MPI info
    const std::vector<I4> &CellsOnCellInit // [in] cell nbrs in linear distrb
) {

   int Err = 0; // initialize return code

   // Retrieve some info on the MPI layout
   MPI_Comm Comm = InEnv->getComm();
   I4 NumTasks   = InEnv->getNumTasks();
   I4 MyTask     = InEnv->getMyTask();
   I4 MasterTask = InEnv->getMasterTask();
   bool IsMaster = InEnv->isMasterTask();

   // TEMPORARY:
   // Due to difficulties with ParMetis, we use serial Metis for now with
   // each task calling the serial form with the global adjacency data.
   // This requires us to communicate the CellsOnCell data and pack it
   // into the Metis structure.

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
   for (int Task = 0; Task < NumTasks; ++Task) {

      // If it is this task's turn, pack up the CellsOnCell data and broadcast
      // to other tasks.
      if (MyTask == Task) {
         for (int n = 0; n < CellsOnCellSize; ++n) {
            CellsOnCellBuf[n] = CellsOnCellInit[n];
         } // end loop CellsOnCell
      } // end if this is MyTask
      Err = MPI_Bcast(&CellsOnCellBuf[0], CellsOnCellSize, MPI_INT32_T, Task,
                      Comm);
      if (Err != 0) {
         LOG_CRITICAL("Decomp: Error communicating CellsOnCell info");
         return Err;
      }

      // Create the adjacency graph by aggregating the individual
      // chunks. Prune edges that don't have neighbors.
      for (int Cell = 0; Cell < NCellsChunk; ++Cell) {

         I4 CellGlob = Task * NCellsChunk + Cell; // global cell ID

         // when chunks do not divide evenly this happens and we break out
         if (CellGlob >= NCellsGlobal)
            break;

         AdjAdd[CellGlob] = Add; // start add for cell in Adjacency array
         for (int Edge = 0; Edge < MaxEdges; ++Edge) {
            I4 BufAdd  = Cell * MaxEdges + Edge;
            I4 NbrCell = CellsOnCellBuf[BufAdd];
            // Skip edges with no neighbors
            if (NbrCell > 0 && NbrCell <= NCellsGlobal) {
               // switch to 0-based indx
               Adjacency[Add] = NbrCell - 1;
               ++Add; // increment address counter
            }
         }
      } // end cell loop for buffer
   } // end task loop

   AdjAdd[NCellsGlobal] = Add; // Add the ending address

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
   int MetisErr = METIS_PartGraphKway(&NCellsMetis, &NConstraints, &AdjAdd[0],
                                      &Adjacency[0], VrtxWgtPtr, VrtxSize,
                                      EdgeWgtPtr, &NumTasksMetis, TpWgts, Ubvec,
                                      Options, &Edgecut, &CellTask[0]);

   if (MetisErr != METIS_OK) {
      LOG_CRITICAL("Decomp: Error in ParMETIS");
      Err = -1;
      return Err;
   }

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
   for (int Task = 0; Task < NumTasks; ++Task)
      TaskCount[Task] = 0;

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
   return Err;

} // end function partCellsKWay

//------------------------------------------------------------------------------
// Partition the edges based on the cell decomposition. The first cell ID in
// the CellsOnEdge array for a given edge is assigned ownership of the edge.

int Decomp::partEdges(
    const MachEnv *InEnv, ///< [in] input machine environment with MPI info
    const std::vector<I4> &CellsOnEdgeInit // [in] cell nbrs for each edge
) {

   I4 Err = 0; // default error code

   // Retrieve some info on the MPI layout
   MPI_Comm Comm = InEnv->getComm();
   I4 NumTasks   = InEnv->getNumTasks();
   I4 MyTask     = InEnv->getMyTask();
   I4 MasterTask = InEnv->getMasterTask();
   bool IsMaster = InEnv->isMasterTask();

   // Calculate some quantities associated with the initial linear
   // distribution
   I4 NCellsChunk = (NCellsGlobal - 1) / NumTasks + 1;
   I4 NEdgesChunk = (NEdgesGlobal - 1) / NumTasks + 1;
   I4 NCellsLocal = NCellsChunk;
   I4 NEdgesLocal = NEdgesChunk;

   // If cells/edges do not divide evenly over processors, the last processor
   // only has the remaining cells and not a full block (chunk)
   if (MyTask == NumTasks - 1) {
      I4 StartAdd = NCellsChunk * (NumTasks - 1);
      NCellsLocal = NCellsGlobal - StartAdd;
      StartAdd    = NEdgesChunk * (NumTasks - 1);
      NEdgesLocal = NEdgesGlobal - StartAdd;
   }

   // For the tasks below, it is useful to keep a sorted list of various cells
   // and edges using the std::set container and related search/sort functions.
   // We create sets for local owned cells, all local edges,
   // local owned edges and local halo edges.
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
      I4 EdgeGlob = MyTask * NEdgesChunk + Edge + 1;
      for (int Cell = 0; Cell < MaxCellsOnEdge; ++Cell) {
         I4 CellGlob = CellsOnEdgeInit[Edge * MaxCellsOnEdge + Cell];
         if (validCellID(CellGlob)) {
            EdgeOwnerInit[Edge] = CellGlob;
            break;
         }
      }
   }

   // Broadcast the edge ownership to all tasks, one chunk at a time
   // and create a sorted list of all owned edges.

   std::vector<I4> EdgeBuf(NEdgesChunk);
   for (int Task = 0; Task < NumTasks; ++Task) {

      // if it is this task's turn, fill the buffer with the owner info
      if (Task == MyTask) {
         for (int Edge = 0; Edge < NEdgesChunk; ++Edge) {
            EdgeBuf[Edge] = EdgeOwnerInit[Edge];
         }
      }
      // Broadcast this buffer
      Err = MPI_Bcast(&EdgeBuf[0], NEdgesChunk, MPI_INT32_T, Task, Comm);

      // For each edge in the buffer, check to see if the task owns
      // the cell. If so, add the edge ID to the owned edges list.
      for (int Edge = 0; Edge < NEdgesChunk; ++Edge) {

         I4 EdgeGlob =
             Task * NEdgesChunk + Edge + 1; // Global ID in initial distrb
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

   } // end task loop

   // For compatibility with the previous MPAS model, we sort the
   // edges based on the order encounted in EdgesOnCell. The first halo
   // level is actually stored in reverse order from the end inward. We sort
   // edge IDs, locations and CellsOnEdge with this ordering.

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

   // Now that we have the local lists, update the final location
   // (task, local edge address) of each of the local edges. This
   // requires one more round of communication.
   // Resize the buffer to make sure we have enough room - the distribution
   // may be less even than the original chunk size.

   HostArray2DI4 EdgeLocTmp("EdgeLoc", NEdgesSize, 2);
   EdgeBuf.resize(2 * NEdgesChunk);

   for (int Edge = 0; Edge < NEdgesSize; ++Edge) {
      EdgeLocTmp(Edge, 0) = MyTask;
      EdgeLocTmp(Edge, 1) = NEdgesAll;
   }

   for (int Task = 0; Task < NumTasks; ++Task) {

      // fill broadcast buffer with the list of owned edges. The
      // first entry in the vector is the number of edges owned by
      // this task.
      if (Task == MyTask) {
         EdgeBuf[0] = NEdgesOwned;
         for (int BufEdge = 0; BufEdge < NEdgesOwned; ++BufEdge) {
            EdgeBuf[BufEdge + 1] = EdgeIDTmp(BufEdge);
         }
      }
      // Broadcast the list of edges owned by this task
      Err = MPI_Bcast(&EdgeBuf[0], 2 * NEdgesChunk, MPI_INT32_T, Task, Comm);

      // For each edge in the buffer, search the full list of edges on
      // this task and store the location
      I4 BufOwned = EdgeBuf[0];
      for (int BufEdge = 0; BufEdge < BufOwned; ++BufEdge) {
         I4 GlobID = EdgeBuf[BufEdge + 1];
         I4 Edge   = srchVector(EdgeIDTmp, GlobID);
         if (Edge < NEdgesAll) {
            EdgeLocTmp(Edge, 0) = Task;    // Task that owns edge
            EdgeLocTmp(Edge, 1) = BufEdge; // Local address on task
         }
      }
   }

   // Copy ID and location arrays into permanent storage
   EdgeIDH     = EdgeIDTmp;
   EdgeLocH    = EdgeLocTmp;
   NEdgesHaloH = NEdgesHaloTmp;

   return Err;

} // end function partEdges

//------------------------------------------------------------------------------
// Partition the vertices based on the cell decomposition. The first cell ID
// in the CellsOnVertex array for a given vertex is assigned ownership of the
// vertex.

int Decomp::partVertices(
    const MachEnv *InEnv, ///< [in] input machine environment with MPI info
    const std::vector<I4> &CellsOnVertexInit // [in] cell nbrs for each vrtx
) {

   I4 Err = 0; // default error code

   // Retrieve some info on the MPI layout
   MPI_Comm Comm = InEnv->getComm();
   I4 NumTasks   = InEnv->getNumTasks();
   I4 MyTask     = InEnv->getMyTask();
   I4 MasterTask = InEnv->getMasterTask();
   bool IsMaster = InEnv->isMasterTask();

   // Calculate some quantities associated with the initial linear
   // distribution
   I4 NCellsChunk    = (NCellsGlobal - 1) / NumTasks + 1;
   I4 NCellsLocal    = NCellsChunk;
   I4 NVerticesChunk = (NVerticesGlobal - 1) / NumTasks + 1;
   I4 NVerticesLocal = NVerticesChunk;

   // If cells/edges do not divide evenly over processors, the last processor
   // only has the remaining cells and not a full block (chunk)
   if (MyTask == NumTasks - 1) {
      I4 StartAdd    = NCellsChunk * (NumTasks - 1);
      NCellsLocal    = NCellsGlobal - StartAdd;
      StartAdd       = NVerticesChunk * (NumTasks - 1);
      NVerticesLocal = NVerticesGlobal - StartAdd;
   }

   // For the tasks below, it is useful to keep a sorted list of various cells
   // and vertices using the std::set container and related search/sort
   // functions. We create sets for local owned cells, all local vertices,
   // local owned vertices and local halo vertices.
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
      I4 VrtxGlob = MyTask * NVerticesChunk + Vrtx + 1;
      for (int Cell = 0; Cell < VertexDegree; ++Cell) {
         I4 CellGlob = CellsOnVertexInit[Vrtx * VertexDegree + Cell];
         if (validCellID(CellGlob)) {
            VrtxOwnerInit[Vrtx] = CellGlob;
            break;
         }
      }
   }

   // Broadcast the vertex ownership to all tasks, one chunk at a time
   // and create a sorted list of all owned vertices.

   std::vector<I4> VrtxBuf(NVerticesChunk);
   for (int Task = 0; Task < NumTasks; ++Task) {

      // if it is this task's turn, fill the buffer with the owner info
      if (Task == MyTask) {
         for (int Vrtx = 0; Vrtx < NVerticesChunk; ++Vrtx) {
            VrtxBuf[Vrtx] = VrtxOwnerInit[Vrtx];
         }
      }
      // Broadcast this buffer
      Err = MPI_Bcast(&VrtxBuf[0], NVerticesChunk, MPI_INT32_T, Task, Comm);

      // For each vertex in the buffer, check to see if the task owns
      // the cell. If so, add the vertex ID to the owned vertices list.
      for (int Vrtx = 0; Vrtx < NVerticesChunk; ++Vrtx) {

         I4 VrtxGlob =
             Task * NVerticesChunk + Vrtx + 1; // Global ID init distrb
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

   } // end task loop

   // For compatibility with the previous MPAS model, we sort the
   // vertices based on the order encounted in VerticesOnCell. The first halo
   // level is actually stored in reverse order from the end inward. We sort
   // vertex IDs, locations and CellsOnVertex with this ordering.

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

   // Now that we have the local lists, update the final location
   // (task, local edge address) of each of the local vertices. This
   // requires one more round of communication.
   // Resize the buffer to make sure we have enough room - the distribution
   // may be less even than the original chunk size.

   HostArray2DI4 VertexLocTmp("VertexLoc", NVerticesSize, 2);
   VrtxBuf.resize(2 * NVerticesChunk);

   for (int Vrtx = 0; Vrtx < NVerticesSize; ++Vrtx) {
      VertexLocTmp(Vrtx, 0) = MyTask;
      VertexLocTmp(Vrtx, 1) = NVerticesAll;
   }

   for (int Task = 0; Task < NumTasks; ++Task) {

      // fill broadcast buffer with the list of owned vertices. The
      // first entry in the vector is the number of vertices owned by
      // this task.
      if (Task == MyTask) {
         VrtxBuf[0] = NVerticesOwned;
         for (int BufVrtx = 0; BufVrtx < NVerticesOwned; ++BufVrtx) {
            VrtxBuf[BufVrtx + 1] = VertexIDTmp(BufVrtx);
         }
      }
      // Broadcast the list of edges owned by this task
      Err = MPI_Bcast(&VrtxBuf[0], 2 * NVerticesChunk, MPI_INT32_T, Task, Comm);

      // For each vertex in the buffer, search the full list of vertices on
      // this task and store the location
      I4 BufOwned = VrtxBuf[0];
      for (int BufVrtx = 0; BufVrtx < BufOwned; ++BufVrtx) {
         I4 GlobID = VrtxBuf[BufVrtx + 1];
         I4 Vrtx   = srchVector(VertexIDTmp, GlobID);
         if (Vrtx < NVerticesAll) {
            VertexLocTmp(Vrtx, 0) = Task;    // Task that owns vertex
            VertexLocTmp(Vrtx, 1) = BufVrtx; // Local address on task
         }
      }
   }

   // Copy ID and location arrays into permanent storage
   VertexIDH      = VertexIDTmp;
   VertexLocH     = VertexLocTmp;
   NVerticesHaloH = NVerticesHaloTmp;

   return Err;

} // end function partVertices

//------------------------------------------------------------------------------
// Redistribute the various XxOnCell index arrays to the final cell
// decomposition. The inputs are the various XxOnCell arrays in the
// initial linear distribution. On exit, all the XxOnCell arrays are
// in the correct final domain decomposition.

int Decomp::rearrangeCellArrays(
    const MachEnv *InEnv, // input machine environment for MPI layout
    const std::vector<I4> &CellsOnCellInit,   //< [in] cell nbrs on each edge
    const std::vector<I4> &EdgesOnCellInit,   //< [in] edges around each cell
    const std::vector<I4> &VerticesOnCellInit //< [in] vertices around cell
) {

   int Err = 0; // default return code

   // Extract some MPI information
   MPI_Comm Comm = InEnv->getComm();
   I4 NumTasks   = InEnv->getNumTasks();
   I4 MyTask     = InEnv->getMyTask();
   I4 MasterTask = InEnv->getMasterTask();
   bool IsMaster = InEnv->isMasterTask();

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

   // Each task will broadcast the cells it owns in the initial linear
   // distribution and all tasks will search that list and extract the
   // entries it owns.
   for (int Task = 0; Task < NumTasks; ++Task) {

      // If it is this task's turn to send, fill the buffer with the local
      // chunk of all three arrays.
      if (MyTask == Task) { // Fill buffer with local chunk
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
      Err = MPI_Bcast(&CellBuf[0], BufSize, MPI_INT32_T, Task, Comm);
      if (Err != 0) {
         LOG_CRITICAL("rearrangeCellArrays: Error broadcasting cell buffer");
         return Err;
      }

      // For each cell in the message buffer, look through the local
      // cell IDs and if this task has the cell in either the owned or
      // halo entries, fill the cell arrays. For edges, we prune the
      // non-active edges and track the number of edges.
      for (int Cell = 0; Cell < NCellsChunk; ++Cell) {

         I4 GlobalID = Task * NCellsChunk + Cell + 1; // IDs are 1-based
         // if NCellsGlobal did not divide evenly, the last task will
         // have global IDs out of range so break out of loop in that case
         if (GlobalID > NCellsGlobal)
            break;

         // Search through local cell list to determine if a match with
         // the global ID
         I4 LocCell = srchVector(CellIDH, GlobalID);
         if (LocCell < NCellsAll) {

            // Local cell needs the info so extract from the buffer
            // into the local address. For edges, we only store the
            // active edges and maintain a count of the edges.
            I4 EdgeCount = 0;
            for (int Edge = 0; Edge < MaxEdges; ++Edge) {
               I4 BufAdd  = Cell * SizePerCell + Edge * 3;
               I4 NbrCell = CellBuf[BufAdd];
               I4 NbrVrtx = CellBuf[BufAdd + 1];
               I4 NbrEdge = CellBuf[BufAdd + 2];
               if (validCellID(NbrCell)) {
                  CellsOnCellTmp(LocCell, Edge) = NbrCell;
               } else {
                  CellsOnCellTmp(LocCell, Edge) = NCellsGlobal + 1;
               }
               if (validVertexID(NbrVrtx)) {
                  VerticesOnCellTmp(LocCell, Edge) = NbrVrtx;
               } else {
                  VerticesOnCellTmp(LocCell, Edge) = NVerticesGlobal + 1;
               }
               if (validEdgeID(NbrEdge)) {
                  EdgesOnCellTmp(LocCell, EdgeCount) = NbrEdge;
                  EdgeCount++;
               }
            }
            NEdgesOnCellTmp(LocCell) = EdgeCount;
         } // end if local cell
      } // end loop over chunk of global cells
   } // end loop over MPI tasks

   // Copy to final location on host - wait to create device copies until
   // the entries are translated to local addresses rather than global IDs
   CellsOnCellH    = CellsOnCellTmp;
   EdgesOnCellH    = EdgesOnCellTmp;
   VerticesOnCellH = VerticesOnCellTmp;
   NEdgesOnCellH   = NEdgesOnCellTmp;

   // All done
   return Err;

} // end function rearrangeCellArrays

//------------------------------------------------------------------------------
// Redistribute the various XxOnEdge index arrays to the final edge
// decomposition. The inputs are the various XxOnEdge arrays in the
// initial linear distribution. On exit, all the XxOnEdge arrays are
// in the correct final domain decomposition.

int Decomp::rearrangeEdgeArrays(
    const MachEnv *InEnv, // input machine environment for MPI layout
    const std::vector<I4> &CellsOnEdgeInit,   //< [in] cell nbrs on each edge
    const std::vector<I4> &EdgesOnEdgeInit,   //< [in] edges around nbr cells
    const std::vector<I4> &VerticesOnEdgeInit //< [in] vertices at edge end
) {

   int Err = 0; // default return code

   // Extract some MPI information
   MPI_Comm Comm = InEnv->getComm();
   I4 NumTasks   = InEnv->getNumTasks();
   I4 MyTask     = InEnv->getMyTask();
   I4 MasterTask = InEnv->getMasterTask();
   bool IsMaster = InEnv->isMasterTask();

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

   // Each task will broadcast the array chunks it owns in the initial linear
   // distribution and all tasks will search that list and extract the
   // entries it owns.
   for (int Task = 0; Task < NumTasks; ++Task) {

      // If it is this task's turn to send, fill the buffer with the local
      // chunk of all three arrays.
      if (MyTask == Task) { // Fill buffer with local chunk
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
      Err = MPI_Bcast(&EdgeBuf[0], BufSize, MPI_INT32_T, Task, Comm);
      if (Err != 0) {
         LOG_CRITICAL("rearrangeEdgeArrays: Error broadcasting edge buffer");
         return Err;
      }

      // For each edge in the message buffer, look through the local
      // edge IDs and if this task has the edge in either the owned or
      // halo entries, fill the edge arrays. For EdgesOnEdge, we prune the
      // non-active edges and track the number of active entries.
      for (int Edge = 0; Edge < NEdgesChunk; ++Edge) {

         I4 GlobalID = Task * NEdgesChunk + Edge + 1; // IDs are 1-based
         // if NEdgesGlobal did not divide evenly, the last task will
         // have global IDs out of range so break out of loop in that case
         if (GlobalID > NEdgesGlobal)
            break;

         // Search through local edge list to determine if a match with
         // the global ID
         I4 LocEdge = srchVector(EdgeIDH, GlobalID);
         if (LocEdge < NEdgesAll) {

            // Local task owns this edge so extract the array info
            // into the local address. For edges, we only store the
            // active edges and maintain a count of the edges.
            I4 BufAdd = Edge * SizePerEdge;
            for (int Cell = 0; Cell < MaxCellsOnEdge; ++Cell) {
               CellsOnEdgeTmp(LocEdge, Cell) = EdgeBuf[BufAdd];
               ++BufAdd;
            }
            for (int Vrtx = 0; Vrtx < 2; ++Vrtx) {
               VerticesOnEdgeTmp(LocEdge, Vrtx) = EdgeBuf[BufAdd];
               ++BufAdd;
            }
            // In the EdgeOnEdge array, a zero entry must be kept in
            // place but assigned the boundary value NEdgesGlobal+1
            I4 EdgeCount = 0;
            for (int NbrEdge = 0; NbrEdge < 2 * MaxEdges; ++NbrEdge) {
               I4 EdgeID = EdgeBuf[BufAdd];
               ++BufAdd;
               if (EdgeID == 0) {
                  EdgesOnEdgeTmp(LocEdge, EdgeCount) = NEdgesGlobal + 1;
                  EdgeCount++;
               } else if (validEdgeID(EdgeID)) {
                  EdgesOnEdgeTmp(LocEdge, EdgeCount) = EdgeID;
                  EdgeCount++;
               }
            }
            NEdgesOnEdgeTmp(LocEdge) = EdgeCount;
         } // end if local cell
      } // end loop over chunk of global cells
   } // end loop over MPI tasks

   // Copy to final location on host - wait to create device copies until
   // the entries are translated to local addresses rather than global IDs
   CellsOnEdgeH    = CellsOnEdgeTmp;
   EdgesOnEdgeH    = EdgesOnEdgeTmp;
   VerticesOnEdgeH = VerticesOnEdgeTmp;
   NEdgesOnEdgeH   = NEdgesOnEdgeTmp;

   // All done
   return Err;

} // end function rearrangeEdgeArrays

//------------------------------------------------------------------------------
// Redistribute the various XxOnVertex index arrays to the final vertex
// decomposition. The inputs are the various XxOnVertex arrays in the
// initial linear distribution. On exit, all the XxOnVertex arrays are
// in the correct final domain decomposition.

int Decomp::rearrangeVertexArrays(
    const MachEnv *InEnv, // input machine environment for MPI layout
    const std::vector<I4> &CellsOnVertexInit, //< [in] cells at each vrtx
    const std::vector<I4> &EdgesOnVertexInit  //< [in] edges joined at vrtx
) {

   int Err = 0; // default return code

   // Extract some MPI information
   MPI_Comm Comm = InEnv->getComm();
   I4 NumTasks   = InEnv->getNumTasks();
   I4 MyTask     = InEnv->getMyTask();
   I4 MasterTask = InEnv->getMasterTask();
   bool IsMaster = InEnv->isMasterTask();

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

   // Each task will broadcast the array chunks it owns in the initial linear
   // distribution and all tasks will search that list and extract the
   // entries it owns.
   for (int Task = 0; Task < NumTasks; ++Task) {

      // If it is this task's turn to send, fill the buffer with the local
      // chunk of both arrays.
      if (MyTask == Task) { // Fill buffer with local chunk
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
      Err = MPI_Bcast(&VrtxBuf[0], BufSize, MPI_INT32_T, Task, Comm);
      if (Err != 0) {
         LOG_CRITICAL("rearrangeVertexArrays: Error broadcasting buffer");
         return Err;
      }

      // For each vertex in the message buffer, look through the local
      // vertex IDs and if this task has the vertex in either the owned or
      // halo entries, fill the vertex arrays.
      for (int Vrtx = 0; Vrtx < NVerticesChunk; ++Vrtx) {

         I4 GlobalID = Task * NVerticesChunk + Vrtx + 1; // IDs are 1-based
         // if NVerticesGlobal did not divide evenly, the last task will
         // have global IDs out of range so break out of loop in that case
         if (GlobalID > NVerticesGlobal)
            break;

         // Search through local vertex list to determine if a match with
         // the global ID
         I4 LocVrtx = srchVector(VertexIDH, GlobalID);
         if (LocVrtx < NVerticesAll) {

            // Local task owns this vertex so extract the array info
            // into the local address.
            I4 BufAdd = Vrtx * SizePerVrtx;
            for (int Cell = 0; Cell < VertexDegree; ++Cell) {
               CellsOnVertexTmp(LocVrtx, Cell) = VrtxBuf[BufAdd];
               ++BufAdd;
            }
            for (int Edge = 0; Edge < VertexDegree; ++Edge) {
               EdgesOnVertexTmp(LocVrtx, Edge) = VrtxBuf[BufAdd];
               ++BufAdd;
            }
         } // end if local cell
      } // end loop over chunk of global cells
   } // end loop over MPI tasks

   // Copy to final location on host - wait to create device copies until
   // the entries are translated to local addresses rather than global IDs
   CellsOnVertexH = CellsOnVertexTmp;
   EdgesOnVertexH = EdgesOnVertexTmp;

   // All done
   return Err;

} // end function rearrangeVertexArrays

//------------------------------------------------------------------------------
// Utility routine to convert a partition method string into PartMethod enum

PartMethod getPartMethodFromStr(const std::string &InMethod) {

   // convert string to lower case for easier equivalence checking
   std::string MethodComp = InMethod;
   std::transform(MethodComp.begin(), MethodComp.end(), MethodComp.begin(),
                  [](unsigned char c) { return std::tolower(c); });

   // Check supported methods and return appropriate enum
   // Currently, only the METIS/ParMETIS KWay option is supported
   if (MethodComp == "metiskway") {
      return PartMethodMetisKWay;

   } else {
      return PartMethodUnknown;

   } // end branch on method string

} // End getPartMethodFromStr

//------------------------------------------------------------------------------
// end Decomp methods

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
