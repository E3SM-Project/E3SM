//===-- base/HorzMesh.cpp - horizontal mesh  methods ----------*- C++ -*-===//
//
//
//
//===----------------------------------------------------------------------===//

#include "Decomp.h"
#include "HorzMesh.h"
#include "DataTypes.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"

namespace OMEGA {

//------------------------------------------------------------------------------
// Construct a new local mesh given a decomposition

HorzMesh::HorzMesh(Decomp *MeshDecomp){

   MeshFileName = MeshDecomp->MeshFileName;
   
   
   // Retrieve mesh cell/edge/vertex totals from Decomp
   NCellsOwned = MeshDecomp->NCellsOwned;
   NCellsAll = MeshDecomp->NCellsAll;
   NCellsSize = MeshDecomp->NCellsSize;
   
   NEdgesOwned = MeshDecomp->NEdgesOwned;
   NEdgesAll = MeshDecomp->NEdgesAll;
   NEdgesSize = MeshDecomp->NEdgesSize;
   MaxCellsOnEdge = MeshDecomp->MaxCellsOnEdge;
   MaxEdges = MeshDecomp->MaxEdges;
   
   NVerticesOwned = MeshDecomp->NVerticesOwned;
   NVerticesAll = MeshDecomp->NVerticesAll;
   NVerticesSize = MeshDecomp->NVerticesSize;
   VertexDegree = MeshDecomp->VertexDegree;
   
   
   // Retrieve connectivity arrays from Decomp
   CellsOnCellH = MeshDecomp->CellsOnCellH;
   EdgesOnCellH = MeshDecomp->EdgesOnCellH;
   NEdgesOnCellH = MeshDecomp->NEdgesOnCellH;
   VerticesOnCellH = MeshDecomp->VerticesOnCellH;
   CellsOnEdgeH = MeshDecomp->CellsOnEdgeH;
   EdgesOnEdgeH = MeshDecomp->EdgesOnEdgeH;
   NEdgesOnEdgeH = MeshDecomp->NEdgesOnEdgeH;
   VerticesOnEdgeH = MeshDecomp->VerticesOnEdgeH;
   CellsOnVertexH = MeshDecomp->CellsOnVertexH;
   EdgesOnVertexH = MeshDecomp->EdgesOnVertexH;
   
   
   // Open the mesh file for reading (assume IO has already been initialized)
   I4 Err;
   Err = OMEGA::IOFileOpen(MeshFileID, MeshFileName, IOModeRead);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error opening mesh file");
   
   
   //Create parallel IO decomposition
   I4 NDims = 1;
   IORearranger Rearr = IORearrBox;
   
   std::vector<I4> CellDims{MeshDecomp->NCellsGlobal};
   std::vector<I4> CellID(NCellsAll);
   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
     CellID[Cell] = MeshDecomp->CellIDH(Cell) - 1;
   }
   
   Err = OMEGA::IOCreateDecomp(CellDecompR8, OMEGA::IOTypeR8, NDims, CellDims,
                               NCellsAll, CellID, Rearr);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error creating cell IO decomposition");
   
   std::vector<I4> EdgeDims{MeshDecomp->NEdgesGlobal};
   std::vector<I4> EdgeID(NEdgesAll);
   for (int Edge = 0; Edge < NEdgesAll; ++Edge) {
     EdgeID[Edge] = MeshDecomp->EdgeIDH(Edge) - 1;
   }
   
   Err = OMEGA::IOCreateDecomp(EdgeDecompR8, OMEGA::IOTypeR8, NDims, EdgeDims,
                               NEdgesAll, EdgeID, Rearr);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error creating edge IO decomposition");
   
   std::vector<I4> VertexDims{MeshDecomp->NVerticesGlobal};
   std::vector<I4> VertexID(NVerticesAll);
   for (int Vertex = 0; Vertex < NVerticesAll; ++Vertex) {
     VertexID[Vertex] = MeshDecomp->VertexIDH(Vertex) - 1;
   }
   
   Err = OMEGA::IOCreateDecomp(VertexDecompR8, OMEGA::IOTypeR8, NDims, VertexDims,
                               NVerticesAll, VertexID, Rearr);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error creating vertex IO decomposition");

   // Read mesh coordinates   
   readCoordinates();

   readBottomDepth();

   readMeasurements();

   readWeights();

   readCoriolis();


} // end constructor

HorzMesh::~HorzMesh() {

   // TODO: add deletes for all arrays and remove from AllDecomps map

} // end deconstructor

void HorzMesh::readCoordinates() {

   I4 Err;
   
   // Read mesh cell coordinates
   XCellH = ArrayHost1DR8("xCell", NCellsAll);
   Err = OMEGA::IOReadArray(XCellH.data(), NCellsAll, "xCell",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading xCell");
   
   YCellH = ArrayHost1DR8("yCell", NCellsAll);
   Err = OMEGA::IOReadArray(YCellH.data(), NCellsAll, "yCell",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading yCell");
   
   ZCellH = ArrayHost1DR8("zCell", NCellsAll);
   Err = OMEGA::IOReadArray(ZCellH.data(), NCellsAll, "zCell",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading zCell");
   
   
   LonCellH = ArrayHost1DR8("lonCell", NCellsAll);
   Err = OMEGA::IOReadArray(LonCellH.data(), NCellsAll, "lonCell",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading lonCell");
   
   LatCellH = ArrayHost1DR8("latCell", NCellsAll);
   Err = OMEGA::IOReadArray(LatCellH.data(), NCellsAll, "latCell",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading latCell");
   
   
   // Read mesh edge coordinates
   XEdgeH = ArrayHost1DR8("xEdge", NEdgesAll);
   Err = OMEGA::IOReadArray(XEdgeH.data(), NEdgesAll, "xEdge",
                            MeshFileID, EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading xEdge");
   
   YEdgeH = ArrayHost1DR8("yEdge", NEdgesAll);
   std::vector<R8> YEdgeTmp(NEdgesAll);
   Err = OMEGA::IOReadArray(YEdgeH.data(), NEdgesAll, "yEdge",
                            MeshFileID, EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading yEdge");
   
   ZEdgeH = ArrayHost1DR8("zEdge", NEdgesAll);
   Err = OMEGA::IOReadArray(ZEdgeH.data(), NEdgesAll, "zEdge",
                            MeshFileID, EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading zEdge");
   
   LonEdgeH = ArrayHost1DR8("lonEdge", NEdgesAll);
   Err = OMEGA::IOReadArray(LonEdgeH.data(), NEdgesAll, "lonEdge",
                            MeshFileID, EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading lonEdge");
   
   LatEdgeH = ArrayHost1DR8("latEdge", NEdgesAll);
   Err = OMEGA::IOReadArray(LatEdgeH.data(), NEdgesAll, "latEdge",
                            MeshFileID, EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading latEdge");

   
   // Read mesh vertex coordinates
   XVertexH = ArrayHost1DR8("xVertex", NVerticesAll);
   Err = OMEGA::IOReadArray(XVertexH.data(), NVerticesAll, "xVertex",
                            MeshFileID, VertexDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading xVertex");
   
   YVertexH = ArrayHost1DR8("yVertex", NVerticesAll);
   Err = OMEGA::IOReadArray(YVertexH.data(), NVerticesAll, "yVertex",
                            MeshFileID, VertexDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading yVertex");
   
   ZVertexH = ArrayHost1DR8("zVertex", NVerticesAll);
   Err = OMEGA::IOReadArray(ZVertexH.data(), NVerticesAll, "zVertex",
                            MeshFileID, VertexDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading zVertex");
   
   LonVertexH = ArrayHost1DR8("lonVertex", NVerticesAll);
   Err = OMEGA::IOReadArray(LonVertexH.data(), NVerticesAll, "lonVertex",
                            MeshFileID, VertexDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading lonVertex");
   
   LatVertexH = ArrayHost1DR8("latVertex", NVerticesAll);
   Err = OMEGA::IOReadArray(LatVertexH.data(), NVerticesAll, "latVertex",
                            MeshFileID, VertexDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading latVertex");

} // end readCoordinates

void HorzMesh::readBottomDepth() {

   I4 Err;

   BottomDepthH = ArrayHost1DR8("bottomDepth", NCellsAll);
   Err = OMEGA::IOReadArray(BottomDepthH.data(), NCellsAll, "bottomDepth",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading bottomDepth");
   
} // end readDepth

void HorzMesh::readMeasurements() {

   I4 Err;

   AreaCellH = ArrayHost1DR8("areaCell", NCellsAll);
   Err = OMEGA::IOReadArray(AreaCellH.data(), NCellsAll, "areaCell",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading areaCell");
   
   AreaTriangleH = ArrayHost1DR8("areaTriangle", NVerticesAll);
   Err = OMEGA::IOReadArray(AreaTriangleH.data(), NVerticesAll, "areaTriangle",
                            MeshFileID, VertexDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading areaTriangle");
   
   DvEdgeH = ArrayHost1DR8("dvEdge", NEdgesAll);
   Err = OMEGA::IOReadArray(DvEdgeH.data(), NEdgesAll, "dvEdge",
                            MeshFileID, EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading dvEdge");
   
   DcEdgeH = ArrayHost1DR8("dcEdge", NEdgesAll);
   Err = OMEGA::IOReadArray(DcEdgeH.data(), NEdgesAll, "dcEdge",
                            MeshFileID, EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading dcEdge");
   
   MeshDensityH = ArrayHost1DR8("meshDensity", NCellsAll);
   Err = OMEGA::IOReadArray(MeshDensityH.data(), NCellsAll, "meshDensity",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading meshDensity");

} // end readMeasurements

void HorzMesh::readWeights() {

} // end readWeights

void HorzMesh::readCoriolis() {

   int Err;

   FCellH = ArrayHost1DR8("fCell", NCellsAll);
   Err = OMEGA::IOReadArray(FCellH.data(), NCellsAll, "fCell",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading fCell");
   
   FVertexH = ArrayHost1DR8("fVertex", NVerticesAll);
   Err = OMEGA::IOReadArray(FVertexH.data(), NVerticesAll, "fVertex",
                            MeshFileID, VertexDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading fVertex");

   FEdgeH = ArrayHost1DR8("fEdge", NEdgesAll);
   Err = OMEGA::IOReadArray(FEdgeH.data(), NEdgesAll, "fEdge",
                            MeshFileID, EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading fEdge");

} // end readCoriolis



} // end namespace OMEGA

//===----------------------------------------------------------------------===//
