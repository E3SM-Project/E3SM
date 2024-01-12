#ifndef OMEGA_HORZMESH_H
#define OMEGA_HORZMESH_H
//===-- base/HorzMesh.h - horizontal mesh --------------------*- C++ -*-===//
//
/// \file
/// \brief Contains the data for an OMEGA horizontal sub-domain
///
/// The HorzMesh class contains the data to represent a sub-domain of the 
/// global horizontal mesh.  
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "MachEnv.h"

#include <string>

namespace OMEGA {

/// A class for the horizontal mesh information

/// The HorzMesh class reads in the remaining mesh information that is
/// not used to create the domain decomposition from mesh file.
/// It handles computing any dependent mesh quantities and transfers 
/// the relevant information to the device
class HorzMesh {

 private:

    void readCoordinates();

    void readBottomDepth();

    //int readMesh();

    //int transferToDevice();

    //int computeMesh();

    I4 CellDecompR8;
    I4 EdgeDecompR8;
    I4 VertexDecompR8;

 public:
   // Variables
   // Since these are used frequently, we make them public to reduce the
   // number of retrievals required.

   std::string MeshFileName;
   int MeshFileID;

   // Sizes and global IDs
   // Note that all sizes are actual counts (1-based) so that loop extents
   // should always use the 0:NCellsXX-1 form.

   I4 NCellsOwned;  ///< Number of cells owned by this task
   I4 NCellsAll;    ///< Total number of local cells (owned + all halo)
   I4 NCellsSize;   ///< Array size (incl padding, bndy cell) for cell arrays

   I4 NEdgesOwned;    ///< Number of edges owned by this task
   I4 NEdgesAll;      ///< Total number (owned+halo) of local edges
   I4 NEdgesSize;     ///< Array length (incl padding, bndy) for edge dim
   I4 MaxCellsOnEdge; ///< Max number of cells sharing an edge
   I4 MaxEdges;       ///< Max number of edges around a cell

   I4 NVerticesOwned;  ///< Number of vertices owned by this task
   I4 NVerticesAll;    ///< Total number (owned+halo) of local vertices
   I4 NVerticesSize;   ///< Array length (incl padding, bndy) for vrtx dim
   I4 VertexDegree;    ///< Number of cells that meet at each vertex

   // Mesh connectivity

   Array2DI4 CellsOnCell;      ///< Indx of cells that neighbor each cell
   ArrayHost2DI4 CellsOnCellH; ///< Indx of cells that neighbor each cell

   Array2DI4 EdgesOnCell;      ///< Indx of edges that border each cell
   ArrayHost2DI4 EdgesOnCellH; ///< Indx of edges that border each cell

   Array1DI4 NEdgesOnCell;      ///< Num of active edges around each cell
   ArrayHost1DI4 NEdgesOnCellH; ///< Num of active edges around each cell

   Array2DI4 VerticesOnCell;      ///< Indx of vertices bordering each cell
   ArrayHost2DI4 VerticesOnCellH; ///< Indx of vertices bordering each cell

   Array2DI4 CellsOnEdge;      ///< Indx of cells straddling each edge
   ArrayHost2DI4 CellsOnEdgeH; ///< Indx of cells straddling each edge

   Array2DI4 EdgesOnEdge;      ///< Indx of edges around cells across each edge
   ArrayHost2DI4 EdgesOnEdgeH; ///< Indx of edges around cells across each edge

   Array1DI4 NEdgesOnEdge;      ///< Num of edges around the cells across edge
   ArrayHost1DI4 NEdgesOnEdgeH; ///< Num of edges around the cells across edge

   Array2DI4 VerticesOnEdge;      ///< Indx of vertices straddling each edge
   ArrayHost2DI4 VerticesOnEdgeH; ///< Indx of vertices straddling each edge

   Array2DI4 CellsOnVertex;      ///< Indx of cells that share a vertex
   ArrayHost2DI4 CellsOnVertexH; ///< Indx of cells that share a vertex

   Array2DI4 EdgesOnVertex;      ///< Indx of edges sharing vertex as endpoint
   ArrayHost2DI4 EdgesOnVertexH; ///< Indx of edges sharing vertex as endpoint

   // Coordinates

   ArrayHost1DR8 XCellH;     ///< X Coordinates of cell cetners (m)
   ArrayHost1DR8 YCellH;     ///< Y Coordinates of cell centers (m)
   ArrayHost1DR8 ZCellH;     ///< Z Coordinates of cell centers (m)
   ArrayHost1DR8 LonCellH;   ///< Longitude location of cell centers (radians)
   ArrayHost1DR8 LatCellH;   ///< Latitude location of cell centers (radians)

   ArrayHost1DR8 XEdgeH;     ///< X Coordinate of edge midpoints (m)
   ArrayHost1DR8 YEdgeH;     ///< Y Coordinate of edge midpoints (m)
   ArrayHost1DR8 ZEdgeH;     ///< Z Coordinate of edge midpoints (m)
   ArrayHost1DR8 LonEdgeH;   ///< Longitude location of edge midpoints (radians)
   ArrayHost1DR8 LatEdgeH;   ///< Latitude location of edge midpoints (radians)

   ArrayHost1DR8 XVertexH;   ///< X Coordinate of vertices (m)
   ArrayHost1DR8 YVertexH;   ///< Y Coordinate of vertices (m)
   ArrayHost1DR8 ZVertexH;   ///< Z Coordinate of vertices (m)
   ArrayHost1DR8 LonVertexH; ///< Longitude location of verticies (radians)
   ArrayHost1DR8 LatVertexH; ///< Latitude location of vertices (radians)

   // Mesh measurements

   Array1DR8 AreaCell;               ///< Area of each cell (m^2)
   ArrayHost1DR8 AreaCellH;          ///< Area of each cell (m^2)

   Array1DR8 AreaTriangle;           ///< Area of each triangle in the dual grid (m^2)
   ArrayHost1DR8 AreaTriangleH;      ///< Area of each triangle in the dual grid (m^2)

   Array2DR8 KiteAreasOnVertex;      ///< Area of the portions of each dual cell that are part of each cellsOnVertex (m^2)
   ArrayHost2DR8 KiteAreasOnVertexH; ///< Area of the portions of each dual cell that are part of each cellsOnVertex (m^2)

   Array1DR8 DvEdge;                 ///< Length of each edge, computed as the distance between verticesOnEdge (m)
   ArrayHost1DR8 DvEdgeH;            ///< Length of each edge, computed as the distance between verticesOnEdge (m)

   Array1DR8 DcEdge;                 ///< Length of each edge, computed as the distance between verticesOnEdge (m)
   ArrayHost1DR8 DcEdgeH;            ///< Length of each edge, computed as the distance between verticesOnEdge (m)

   Array1DR8 AngleEdge;              ///< Angle the edge normal makes with local eastward direction )radians)
   ArrayHost1DR8 AngleEdgeH;         ///< Angle the edge normal makes with local eastward direction )radians)

   ArrayHost1DR8 MeshDensityH;       ///< Value of density function used to generate a particular mesh at cell centers

   // Weights

   Array2DR8 WeightsOnEdge;      ///< Reconstruction weights associated with each of the edgesOnEdge
   ArrayHost2DR8 WeightsOnEdgeH; ///< Reconstruction weights associated with each of the edgesOnEdge
  
   // Coriolis 

   Array1DR8 FEdge;   ///< Coriolis parameter at edges (radians s^-1)
   ArrayHost1DR8 FEdgeH;   ///< Coriolis parameter at edges (radians s^-1)

   Array1DR8 FCell;   ///< Coriolis parameter at cell centers (radians s^-1)
   ArrayHost1DR8 FCellH;   ///< Coriolis parameter at cell centers (radians s^-1)

   Array1DR8 FVertex; ///< Coriolis parameter at vertices (radians s^-1)
   ArrayHost1DR8 FVertexH; ///< Coriolis parameter at vertices (radians s^-1)

   // Depth

   Array1DR8 BottomDepth;      ///< Depth of the bottom of the ocean (m)
   ArrayHost1DR8 BottomDepthH; ///< Depth of the bottom of the ocean (m)
  
   // Methods

   /// Construct a new local mesh for a given decomposition
   HorzMesh(Decomp *Decomp ///< [in] Decomposition for mesh
   );

   /// Destructor - deallocates all memory and deletes a HorzMesh
   ~HorzMesh();

}; // end class HorzMesh 

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_HORZMESH_H
