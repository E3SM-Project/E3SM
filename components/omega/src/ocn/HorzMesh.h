#ifndef OMEGA_HORZMESH_H
#define OMEGA_HORZMESH_H
//===-- base/HorzMesh.h - horizontal mesh --------------------*- C++ -*-===//
//
/// \file
/// \brief Contains the mesh variables for an OMEGA horizontal sub-domain
///
/// The HorzMesh class contains the data to represent a sub-domain of the
/// global horizontal mesh.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "Decomp.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"

#include <memory>
#include <string>

namespace OMEGA {

/// A class for the horizontal mesh information

/// The HorzMesh class reads in the remaining mesh information that is
/// not used to create the domain decomposition from mesh file.
/// It handles computing any dependent mesh quantities and transfers
/// the relevant information to the device
class HorzMesh {

 private:
   void initParallelIO(Decomp *MeshDecomp);

   void finalizeParallelIO();

   void createDimensions(Decomp *MeshDecomp);

   void readVertexArray(HostArray1DReal &VertexArrayH,
                        const std::string &MPASName);
   void readEdgeArray(HostArray1DReal &EdgeArrayH, const std::string &MPASName);
   void readCellArray(HostArray1DReal &CellArrayH, const std::string &MPASName);

   void readCoordinates();

   void readBottomDepth();

   void readMeasurements();

   void readWeights();

   void readCoriolis();

   // void computeEdgeSign();

   void copyToDevice();

   // int computeMesh();
   I4 CellDecompR8;
   I4 EdgeDecompR8;
   I4 VertexDecompR8;
   I4 OnEdgeDecompR8;
   I4 OnVertexDecompR8;

   static HorzMesh *DefaultHorzMesh;

   static std::map<std::string, std::unique_ptr<HorzMesh>> AllHorzMeshes;

   /// Construct a new local mesh for a given decomposition
   HorzMesh(const std::string &Name, ///< [in] Name for mesh
            Decomp *Decomp,          ///< [in] Decomposition for mesh
            I4 InNVertLevels         ///< [in] num vertical levels
   );

   // Forbid copy and move construction
   HorzMesh(const HorzMesh &) = delete;
   HorzMesh(HorzMesh &&)      = delete;

 public:
   // KOKKOS_LAMBDA does not allow to have parallel_* functions inside of a
   // private function.
   void computeEdgeSign();

   void setMasks(int NVertLevels);

   void setMeshScaling();

   // Variables
   // Since these are used frequently, we make them public to reduce the
   // number of retrievals required.

   std::string MeshName;
   std::string MeshFileName;
   int MeshFileID;

   // Sizes and global IDs
   // Note that all sizes are actual counts (1-based) so that loop extents
   // should always use the 0:NCellsXX-1 form.

   I4 NVertLevels; ///< number of vertical levels

   Array1DI4 NCellsHalo;      ///< num cells owned+halo for halo layer
   HostArray1DI4 NCellsHaloH; ///< num cells owned+halo for halo layer
   I4 NCellsOwned;            ///< Number of cells owned by this task
   I4 NCellsAll;  ///< Total number of local cells (owned + all halo)
   I4 NCellsSize; ///< Array size (incl padding, bndy cell) for cell arrays

   Array1DI4 NEdgesHalo;      ///< num cells owned+halo for halo layer
   HostArray1DI4 NEdgesHaloH; ///< num cells owned+halo for halo layer
   I4 NEdgesOwned;            ///< Number of edges owned by this task
   I4 NEdgesAll;              ///< Total number (owned+halo) of local edges
   I4 NEdgesSize;     ///< Array length (incl padding, bndy) for edge dim
   I4 MaxCellsOnEdge; ///< Max number of cells sharing an edge
   I4 MaxEdges;       ///< Max number of edges around a cell
   I4 MaxEdges2;      ///< Max number of edges around a cell x2

   Array1DI4 NVerticesHalo;      ///< num cells owned+halo for halo layer
   HostArray1DI4 NVerticesHaloH; ///< num cells owned+halo for halo layer
   I4 NVerticesOwned;            ///< Number of vertices owned by this task
   I4 NVerticesAll;  ///< Total number (owned+halo) of local vertices
   I4 NVerticesSize; ///< Array length (incl padding, bndy) for vrtx dim
   I4 VertexDegree;  ///< Number of cells that meet at each vertex

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

   // Coordinates

   Array1DReal XCell;      ///< X Coordinates of cell centers (m)
   HostArray1DReal XCellH; ///< X Coordinates of cell centers (m)

   Array1DReal YCell;      ///< Y Coordinates of cell centers (m)
   HostArray1DReal YCellH; ///< Y Coordinates of cell centers (m)

   HostArray1DReal ZCellH;   ///< Z Coordinates of cell centers (m)
   HostArray1DReal LonCellH; ///< Longitude location of cell centers (radians)
   HostArray1DReal LatCellH; ///< Latitude location of cell centers (radians)

   Array1DReal XEdge;      ///< X Coordinate of edge midpoints (m)
   HostArray1DReal XEdgeH; ///< X Coordinate of edge midpoints (m)

   Array1DReal YEdge;      ///< Y Coordinate of edge midpoints (m)
   HostArray1DReal YEdgeH; ///< Y Coordinate of edge midpoints (m)

   HostArray1DReal ZEdgeH;   ///< Z Coordinate of edge midpoints (m)
   HostArray1DReal LonEdgeH; ///< Longitude location of edge midpoints (radians)
   HostArray1DReal LatEdgeH; ///< Latitude location of edge midpoints (radians)

   HostArray1DReal XVertexH;   ///< X Coordinate of vertices (m)
   HostArray1DReal YVertexH;   ///< Y Coordinate of vertices (m)
   HostArray1DReal ZVertexH;   ///< Z Coordinate of vertices (m)
   HostArray1DReal LonVertexH; ///< Longitude location of vertices (radians)
   HostArray1DReal LatVertexH; ///< Latitude location of vertices (radians)

   // Mesh measurements

   Array1DReal AreaCell;      ///< Area of each cell (m^2)
   HostArray1DReal AreaCellH; ///< Area of each cell (m^2)

   Array1DReal AreaTriangle; ///< Area of each triangle in the dual grid (m^2)
   HostArray1DReal
       AreaTriangleH; ///< Area of each triangle in the dual grid (m^2)

   Array2DReal
       KiteAreasOnVertex; ///< Area of the portions of each dual cell that
                          ///  are part of each cellsOnVertex (m^2)
   HostArray2DReal
       KiteAreasOnVertexH; ///< Area of the portions of each dual cell that
                           ///  are part of each cellsOnVertex (m^2)

   Array1DReal
       DvEdge; ///< Length of each edge, computed as the distance between
               ///  verticesOnEdge (m)
   HostArray1DReal DvEdgeH; ///< Length of each edge, computed as the distance
                            ///  between verticesOnEdge (m)

   Array1DReal
       DcEdge; ///< Length of each edge, computed as the distance between
               ///  CellsOnEdge (m)
   HostArray1DReal DcEdgeH; ///< Length of each edge, computed as the distance
                            ///  between CellsOnEdge (m)

   Array1DReal AngleEdge; ///< Angle the edge normal makes with local eastward
                          ///  direction (radians)
   HostArray1DReal AngleEdgeH; ///< Angle the edge normal makes with local
                               ///  eastward direction (radians)

   HostArray1DReal
       MeshDensityH; ///< Value of density function used to generate a
                     ///  particular mesh at cell centers

   // Weights

   Array2DReal
       WeightsOnEdge; ///< Reconstruction weights associated with each of
                      ///  the edgesOnEdge
   HostArray2DReal WeightsOnEdgeH; ///< Reconstruction weights associated with
                                   ///  each of the edgesOnEdge

   // Coriolis

   Array1DReal FEdge;      ///< Coriolis parameter at edges (radians s^-1)
   HostArray1DReal FEdgeH; ///< Coriolis parameter at edges (radians s^-1)

   Array1DReal FCell; ///< Coriolis parameter at cell centers (radians s^-1)
   HostArray1DReal
       FCellH; ///< Coriolis parameter at cell centers (radians s^-1)

   Array1DReal FVertex;      ///< Coriolis parameter at vertices (radians s^-1)
   HostArray1DReal FVertexH; ///< Coriolis parameter at vertices (radians s^-1)

   // Depth

   Array1DReal BottomDepth;      ///< Depth of the bottom of the ocean (m)
   HostArray1DReal BottomDepthH; ///< Depth of the bottom of the ocean (m)

   // Edge sign

   Array2DReal EdgeSignOnCell;      ///< Sign of vector connecting cells
   HostArray2DReal EdgeSignOnCellH; ///< Sign of vector connecting cells

   Array2DReal EdgeSignOnVertex;      ///< Sign of vector connecting vertices
   HostArray2DReal EdgeSignOnVertexH; ///< Sign of vector connecting vertices

   // Masks
   Array2DReal EdgeMask;      ///< Mask to determine if computations should be
                              ///  done on edge
   HostArray2DReal EdgeMaskH; ///< Mask to determine if computations should be
                              ///  done on edge

   // Mesh scaling
   Array1DReal MeshScalingDel2;      /// Coef to Laplacian mixing terms
   HostArray1DReal MeshScalingDel2H; /// Coef to Laplacian mixing terms

   Array1DReal MeshScalingDel4;      /// Coef to biharmonic mixing terms
   HostArray1DReal MeshScalingDel4H; /// Coef to biharmonic mixing terms

   // Methods

   /// Initialize Omega local mesh
   static int init();

   /// Creates a new mesh by calling the constructor and puts it in the
   /// AllHorzMeshes map
   static HorzMesh *create(const std::string &Name, ///< [in] Name for mesh
                           Decomp *Decomp,  ///< [in] Decomposition for mesh
                           I4 InNVertLevels ///< [in] num vertival levels
   );

   /// Destructor - deallocates all memory and deletes a HorzMesh
   ~HorzMesh();

   /// Deallocates arrays
   static void clear();

   /// Remove mesh by name
   static void erase(std::string InName ///< [in] name of mesh to remove
   );

   static HorzMesh *getDefault();

   static HorzMesh *get(std::string name);

}; // end class HorzMesh

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_HORZMESH_H
