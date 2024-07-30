#ifndef OMEGA_DECOMP_H
#define OMEGA_DECOMP_H
//===-- base/Decomp.h - domain decomposition --------------------*- C++ -*-===//
//
/// \file
/// \brief Decomposes/partitions an OMEGA horizontal domain
///
/// The decomposition (Decomp) class partitions an OMEGA horizontal domain into
/// a number of sub-domains that can be distributed across a parallel
/// environment. Currently, this relies on the (Par)Metis partitioning package
/// to partition a set of nodes (cells) based on the adjacency graph created
/// during the OMEGA mesh generators. A default decomposition is created
/// from the default Machine Env that specifies the MPI details and number
/// of MPI tasks. Other decompositions can be created with the same mesh
/// on any of the subset environments that are possible in MachEnv.
/// The Decomp class stores a number of index-space arrays that describe
/// the partition, neighbor information, global IDs for cell, edge and
/// vertex points in an Omega mesh.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "MachEnv.h"
#include "mpi.h"
#include "parmetis.h"

#include <memory>
#include <string>

namespace OMEGA {

/// A class for the input adjacency graph file contains information
/// related to the file itself as well as the contents and size of
/// the adjacency graph.

/// Supported partitioning methods
enum PartMethod {
   PartMethodUnknown,   ///< Unknown or undefined method
   PartMethodMetisKWay, ///< Metis K-way partitioning (default)
   PartMethodMetisRB    ///< Metis recursive bisection (not yet supported)
};

/// Translates an input string for partition method option to the
/// enum for later use
PartMethod getPartMethodFromStr(
    const std::string &InMethod ///< [in] choice of partition method
);

/// The Decomp class creates and maintains most of the information related
/// to the mesh index space and its distribution across partitions or processors
/// in a parallel domain decomposition. This information includes the location
/// of every cell/edge/vertex point and the adjacency or neighbor information
/// for all mesh locations. It relies on the ParMetis package for performing
/// the actual partitions.
class Decomp {

 private:
   /// The default decomposition describes the index space decomposition
   /// used for most of the model. Because it is used most often,
   /// we store the extra pointer here for easier retrieval.
   static Decomp *DefaultDecomp;

   /// All decompositions are tracked/stored within the class as a
   /// map paired with a name for later retrieval.
   static std::map<std::string, std::unique_ptr<Decomp>> AllDecomps;

   /// Partition cells by calling the METIS/ParMETIS KWay routine
   /// It starts with the CellsOnCell array from the input mesh file
   /// distributed across tasks in linear contiguous chunks
   /// On output, it has defined all the NCells sizes (NCellsOwned,
   /// NCellsHalo array, NCellsAll and NCellsSize) and the final CellID
   /// and CellLoc arrays
   int partCellsKWay(
       const MachEnv *InEnv,                  ///< [in] MachEnv with MPI info
       const std::vector<I4> &CellsOnCellInit ///< [in] cell nbrs in init dstrb
   );

   /// Partition the edges given the cell partition and edge connectivity
   /// The first cell ID associated with an edge in the CellsOnEdge array
   /// is assumed to own the edge. The inputs are the edge-cell connectivity
   /// arrays that have been initially partitioned across MPI tasks in a
   /// linear distribution. On return, this function populates all of the
   /// NEdge sizes, edge halo indices, and final edge-related connectivity
   /// arrays.
   int partEdges(
       const MachEnv *InEnv,                  ///< [in] MachEnv with MPI info
       const std::vector<I4> &CellsOnEdgeInit ///< [in] cell nbrs on each edge
   );

   /// Partition the vertices given the cell partition and vertex connectivity
   /// The first cell ID associated with an vertex in the CellsOnVertex array
   /// is assumed to own the vertex. The inputs are the vertex-cell connectivity
   /// arrays that have been initially partitioned across MPI tasks in a
   /// linear distribution. On return, this function populates all of the
   /// NVertex sizes, vertex halo indices, and final vertex-related connectivity
   /// arrays.
   int partVertices(
       const MachEnv *InEnv,                    ///< [in] MachEnv with MPI info
       const std::vector<I4> &CellsOnVertexInit ///< [in] cells at each vertex
   );

   /// Redistribute the various XxOnCell index arrays to the final cell
   /// decomposition. The inputs are the various XxOnCell arrays in the
   /// initial linear distribution. On exit, all the XxOnCell arrays are
   /// in the correct final domain decomposition.
   int rearrangeCellArrays(
       const MachEnv *InEnv, ///< [in] MachEnv for the new partition
       const std::vector<I4> &CellsOnCellInit, ///< [in] cell nbrs on each edge
       const std::vector<I4> &EdgesOnCellInit, ///< [in] edges around each cell
       const std::vector<I4> &VerticesOnCellInit ///< [in] vertices around cell
   );

   /// Redistribute the various XxOnEdge index arrays to the final edge
   /// decomposition. The inputs are the various XxOnEdge arrays in the
   /// initial linear distribution. On exit, all the XxOnEdge arrays are
   /// in the correct final domain decomposition.
   int rearrangeEdgeArrays(
       const MachEnv *InEnv, ///< [in] MachEnv for the new partition
       const std::vector<I4> &CellsOnEdgeInit, ///< [in] cell nbrs on each edge
       const std::vector<I4> &EdgesOnEdgeInit, ///< [in] edges around each cell
       const std::vector<I4> &VerticesOnEdgeInit ///< [in] vertices around cell
   );

   /// Redistribute the various XxOnVertex index arrays to the final vertex
   /// decomposition. The inputs are the various XxOnVertex arrays in the
   /// initial linear distribution. On exit, all the XxOnVertex arrays are
   /// in the correct final domain decomposition.
   int rearrangeVertexArrays(
       const MachEnv *InEnv, ///< [in] MachEnv for the new partition
       const std::vector<I4> &CellsOnVertexInit, ///< [in] cells at each vertex
       const std::vector<I4> &EdgesOnVertexInit  ///< [in] edges at each vertex
   );

   /// Construct a new decomposition across an input MachEnv with
   /// NPart partitions of a mesh that is read from a mesh file.
   Decomp(const std::string &Name, ///< [in] Name for new decomposition
          const MachEnv *InEnv,    ///< [in] MachEnv for the new partition
          I4 NParts,               ///< [in] num of partitions for new decomp
          PartMethod Method,       ///< [in] method for partitioning
          I4 InHaloWidth,          ///< [in] width of halo in new decomp
          const std::string &MeshFileName_ ///< [in] name of file with mesh info
   );

   // forbid copy and move construction
   Decomp(const Decomp &) = delete;
   Decomp(Decomp &&)      = delete;

 public:
   // Variables
   // Since these are used frequently, we make them public to reduce the
   // number of retrievals required.

   std::string MeshFileName; ///< The name of the file with mesh info

   // Sizes and global IDs
   // Note that all sizes are actual counts (1-based) so that loop extents
   // should always use the 0:NCellsXX-1 form.

   I4 HaloWidth; ///< Number of halo layers for cell-based variables

   I4 NCellsGlobal; ///< Number of cells in the full global mesh
   I4 NCellsOwned;  ///< Number of cells owned by this task
   I4 NCellsAll;    ///< Total number of local cells (owned + all halo)
   I4 NCellsSize;   ///< Array size (incl padding, bndy cell) for cell arrays
   I4 MaxEdges;     ///< Max number of edges around a cell

   Array1DI4 NCellsHalo;      ///< num cells owned+halo for halo layer
   HostArray1DI4 NCellsHaloH; ///< num cells owned+halo for halo layer
   Array1DI4 CellID;          ///< global cell ID for each local cell
   HostArray1DI4 CellIDH;     ///< global cell ID for each local cell
   Array2DI4 CellLoc;      ///< location (task, local add) for local cells,halo
   HostArray2DI4 CellLocH; ///< location (task, local add) for local cells,halo

   I4 NEdgesGlobal;   ///< Number of edges in the full global mesh
   I4 NEdgesOwned;    ///< Number of edges owned by this task
   I4 NEdgesAll;      ///< Total number (owned+halo) of local edges
   I4 NEdgesSize;     ///< Array length (incl padding, bndy) for edge dim
   I4 MaxCellsOnEdge; ///< Max number of cells sharing an edge

   Array1DI4 NEdgesHalo;      ///< num cells owned+halo for halo layer
   HostArray1DI4 NEdgesHaloH; ///< num cells owned+halo for halo layer
   Array1DI4 EdgeID;          ///< global cell ID for each local cell
   HostArray1DI4 EdgeIDH;     ///< global cell ID for each local cell
   Array2DI4 EdgeLoc;      ///< location (task, local add) for local edges,halo
   HostArray2DI4 EdgeLocH; ///< location (task, local add) for local edges,halo

   I4 NVerticesGlobal; ///< Number of vertices in the full global mesh
   I4 NVerticesOwned;  ///< Number of vertices owned by this task
   I4 NVerticesAll;    ///< Total number (owned+halo) of local vertices
   I4 NVerticesSize;   ///< Array length (incl padding, bndy) for vrtx dim
   I4 VertexDegree;    ///< Number of cells that meet at each vertex

   Array1DI4 NVerticesHalo;      ///< num cells owned+halo for halo layer
   HostArray1DI4 NVerticesHaloH; ///< num cells owned+halo for halo layer
   Array1DI4 VertexID;           ///< global vertex ID for each local cell
   HostArray1DI4 VertexIDH;      ///< global vertex ID for each local cell
   Array2DI4 VertexLoc;      ///< location (task, local add) for local vrtx halo
   HostArray2DI4 VertexLocH; ///< location (task, local add) for local vrtx halo

   // Mesh connectivity

   Array2DI4 CellsOnCell;      ///< Indx of cells that neighbor each cell
   HostArray2DI4 CellsOnCellH; ///< Indx of cells that neighbor each cell

   Array2DI4 EdgesOnCell;      ///< Indx of edges that border each cell
   HostArray2DI4 EdgesOnCellH; ///< Indx of edges that border each cell

   Array1DI4 NEdgesOnCell;      ///< Num of active edges around each cell
   HostArray1DI4 NEdgesOnCellH; ///< Num of active edges around each cell

   Array2DI4 VerticesOnCell;      ///< Indx of vertices bordering each cell
   HostArray2DI4 VerticesOnCellH; ///< Indx of vertices bordering each cell

   Array2DI4 CellsOnEdge;      ///< Indx of cells straddling each edge
   HostArray2DI4 CellsOnEdgeH; ///< Indx of cells straddling each edge

   Array2DI4 EdgesOnEdge;      ///< Indx of edges around cells across each edge
   HostArray2DI4 EdgesOnEdgeH; ///< Indx of edges around cells across each edge

   Array1DI4 NEdgesOnEdge;      ///< Num of edges around the cells across edge
   HostArray1DI4 NEdgesOnEdgeH; ///< Num of edges around the cells across edge

   Array2DI4 VerticesOnEdge;      ///< Indx of vertices straddling each edge
   HostArray2DI4 VerticesOnEdgeH; ///< Indx of vertices straddling each edge

   Array2DI4 CellsOnVertex;      ///< Indx of cells that share a vertex
   HostArray2DI4 CellsOnVertexH; ///< Indx of cells that share a vertex

   Array2DI4 EdgesOnVertex;      ///< Indx of edges sharing vertex as endpoint
   HostArray2DI4 EdgesOnVertexH; ///< Indx of edges sharing vertex as endpoint

   // Methods

   /// Initializes Omega decomposition info and creates the default
   /// decomposition based on the default MachEnv and configuration
   /// options.
   static int init(const std::string &MeshFileName = "OmegaMesh.nc");

   // Creates a new decomposition using the constructor and puts it in the
   // AllDecomps map
   static Decomp *
   create(const std::string &Name, ///< [in] Name for new decomposition
          const MachEnv *Env,      ///< [in] MachEnv for the new partition
          I4 NParts,               ///< [in] num of partitions for new decomp
          PartMethod Method,       ///< [in] method for partitioning
          I4 HaloWidth,            ///< [in] width of halo in new decomp
          const std::string &MeshFileName ///< [in] name of file with mesh info
   );

   /// Destructor - deallocates all memory and deletes a Decomp.
   ~Decomp();

   /// Erase - removes a defined decomposition
   static void erase(std::string InName ///< [in] name of decomp to remove
   );

   /// Clear - removes all defined decompositions to clean up
   static void clear();

   // Retrieval functions

   /// Retrieves a pointer the default decomposition. The preference is
   /// to pass the decomp as an argument, but a retrieval is necessary
   /// for sharing the info between initialization and multiple run phases
   static Decomp *getDefault();

   /// Retrieve a decomposition by name.
   static Decomp *get(std::string name);

   /// Query functions

   /// Checks a global cell ID to make sure it is in the valid range
   bool validCellID(I4 InCellID ///< [in] a cell ID to check
   );

   /// Checks a global edge ID to make sure it is in the valid range
   bool validEdgeID(I4 InEdgeID ///< [in] a edge ID to check
   );

   /// Checks a global vertex ID to make sure it is in the valid range
   bool validVertexID(I4 InVertexID ///< [in] a vertex ID to check
   );

}; // end class Decomp

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_DECOMP_H
