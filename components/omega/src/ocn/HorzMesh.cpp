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
   std::vector<R8> XCellTmp(NCellsAll);
   Err = OMEGA::IOReadArray(&XCellTmp[0], NCellsAll, "xCell",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading xCell");
   
   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      XCellH(Cell) = XCellTmp[Cell];
   }
   
   
   YCellH = ArrayHost1DR8("yCell", NCellsAll);
   std::vector<R8> YCellTmp(NCellsAll);
   Err = OMEGA::IOReadArray(&YCellTmp[0], NCellsAll, "yCell",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading yCell");
   
   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      YCellH(Cell) = YCellTmp[Cell];
   }
   
   
   ZCellH = ArrayHost1DR8("zCell", NCellsAll);
   std::vector<R8> ZCellTmp(NCellsAll);
   Err = OMEGA::IOReadArray(&ZCellTmp[0], NCellsAll, "zCell",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading zCell");
   
   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      ZCellH(Cell) = ZCellTmp[Cell];
   }
   
   
   LonCellH = ArrayHost1DR8("lonCell", NCellsAll);
   std::vector<R8> LonCellTmp(NCellsAll);
   Err = OMEGA::IOReadArray(&LonCellTmp[0], NCellsAll, "lonCell",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading lonCell");
   
   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      LonCellH(Cell) = LonCellTmp[Cell];
   }
   
   
   LatCellH = ArrayHost1DR8("latCell", NCellsAll);
   std::vector<R8> LatCellTmp(NCellsAll);
   Err = OMEGA::IOReadArray(&LatCellTmp[0], NCellsAll, "latCell",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading latCell");
   
   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      LatCellH(Cell) = LatCellTmp[Cell];
   }
   
   // Read mesh edge coordinates
   XEdgeH = ArrayHost1DR8("xEdge", NEdgesAll);
   std::vector<R8> XEdgeTmp(NEdgesAll);
   Err = OMEGA::IOReadArray(&XEdgeTmp[0], NEdgesAll, "xEdge",
                            MeshFileID, EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading xEdge");
   
   for (int Edge = 0; Edge < NEdgesAll; ++Edge) {
      XEdgeH(Edge) = XEdgeTmp[Edge];
   }
   
   
   YEdgeH = ArrayHost1DR8("yEdge", NEdgesAll);
   std::vector<R8> YEdgeTmp(NEdgesAll);
   Err = OMEGA::IOReadArray(&YEdgeTmp[0], NEdgesAll, "yEdge",
                            MeshFileID, EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading yEdge");
   
   for (int Edge = 0; Edge < NEdgesAll; ++Edge) {
      YEdgeH(Edge) = YEdgeTmp[Edge];
   }
   
   
   ZEdgeH = ArrayHost1DR8("zEdge", NEdgesAll);
   std::vector<R8> ZEdgeTmp(NEdgesAll);
   Err = OMEGA::IOReadArray(&ZEdgeTmp[0], NEdgesAll, "zEdge",
                            MeshFileID, EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading zEdge");
   
   for (int Edge = 0; Edge < NEdgesAll; ++Edge) {
      ZEdgeH(Edge) = ZEdgeTmp[Edge];
   }
   
   
   LonEdgeH = ArrayHost1DR8("lonEdge", NEdgesAll);
   std::vector<R8> LonEdgeTmp(NEdgesAll);
   Err = OMEGA::IOReadArray(&LonEdgeTmp[0], NEdgesAll, "lonEdge",
                            MeshFileID, EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading lonEdge");
   
   for (int Edge = 0; Edge < NEdgesAll; ++Edge) {
      LonEdgeH(Edge) = LonEdgeTmp[Edge];
   }
   
   
   LatEdgeH = ArrayHost1DR8("latEdge", NEdgesAll);
   std::vector<R8> LatEdgeTmp(NEdgesAll);
   Err = OMEGA::IOReadArray(&LatEdgeTmp[0], NEdgesAll, "latEdge",
                            MeshFileID, EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading latEdge");
   
   for (int Edge = 0; Edge < NEdgesAll; ++Edge) {
      LatEdgeH(Edge) = LatEdgeTmp[Edge];
   }
   
   // Read mesh vertex coordinates
   XVertexH = ArrayHost1DR8("xVertex", NVerticesAll);
   std::vector<R8> XVertexTmp(NVerticesAll);
   Err = OMEGA::IOReadArray(&XVertexTmp[0], NVerticesAll, "xVertex",
                            MeshFileID, VertexDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading xVertex");
   
   for (int Vertex = 0; Vertex < NVerticesAll; ++Vertex) {
      XVertexH(Vertex) = XVertexTmp[Vertex];
   }
   
   
   YVertexH = ArrayHost1DR8("yVertex", NVerticesAll);
   std::vector<R8> YVertexTmp(NVerticesAll);
   Err = OMEGA::IOReadArray(&YVertexTmp[0], NVerticesAll, "yVertex",
                            MeshFileID, VertexDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading yVertex");
   
   for (int Vertex = 0; Vertex < NVerticesAll; ++Vertex) {
      YVertexH(Vertex) = YVertexTmp[Vertex];
   }
   
   
   ZVertexH = ArrayHost1DR8("zVertex", NVerticesAll);
   std::vector<R8> ZVertexTmp(NVerticesAll);
   Err = OMEGA::IOReadArray(&ZVertexTmp[0], NVerticesAll, "zVertex",
                            MeshFileID, VertexDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading zVertex");
   
   for (int Vertex = 0; Vertex < NVerticesAll; ++Vertex) {
      ZVertexH(Vertex) = ZVertexTmp[Vertex];
   }
   
   
   LonVertexH = ArrayHost1DR8("lonVertex", NVerticesAll);
   std::vector<R8> LonVertexTmp(NVerticesAll);
   Err = OMEGA::IOReadArray(&LonVertexTmp[0], NVerticesAll, "lonVertex",
                            MeshFileID, VertexDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading lonVertex");
   
   for (int Vertex = 0; Vertex < NVerticesAll; ++Vertex) {
      LonVertexH(Vertex) = LonVertexTmp[Vertex];
   }
   
   
   LatVertexH = ArrayHost1DR8("latVertex", NVerticesAll);
   std::vector<R8> LatVertexTmp(NVerticesAll);
   Err = OMEGA::IOReadArray(&LatVertexTmp[0], NVerticesAll, "latVertex",
                            MeshFileID, VertexDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading latVertex");
   
   for (int Vertex = 0; Vertex < NVerticesAll; ++Vertex) {
      LatVertexH(Vertex) = LatVertexTmp[Vertex];
   }

} // end readCoordinates

void HorzMesh::readBottomDepth() {

   I4 Err;

   BottomDepthH = ArrayHost1DR8("bottomDepth", NCellsAll);
   std::vector<R8> BottomDepthTmp(NCellsAll);
   Err = OMEGA::IOReadArray(&BottomDepthTmp[0], NCellsAll, "bottomDepth",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading bottomDepth");
   
   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      BottomDepthH(Cell) = BottomDepthTmp[Cell];
   }
   
} // end readDepth

void HorzMesh::readMeasurements() {

   I4 Err;

   AreaCellH = ArrayHost1DR8("areaCell", NCellsAll);
   std::vector<R8> AreaCellTmp(NCellsAll);
   Err = OMEGA::IOReadArray(&AreaCellTmp[0], NCellsAll, "areaCell",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading areaCell");
   
   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      AreaCellH(Cell) = AreaCellTmp[Cell];
   }

   AreaTriangleH = ArrayHost1DR8("areaTriangle", NVerticesAll);
   std::vector<R8> AreaTriangleTmp(NVerticesAll);
   Err = OMEGA::IOReadArray(&AreaTriangleTmp[0], NVerticesAll, "areaTriangle",
                            MeshFileID, VertexDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading areaTriangle");
   
   for (int Vertex = 0; Vertex < NVerticesAll; ++Vertex) {
      AreaTriangleH(Vertex) = AreaTriangleTmp[Vertex];
   }

   DvEdgeH = ArrayHost1DR8("dvEdge", NEdgesAll);
   std::vector<R8> DvEdgeTmp(NEdgesAll);
   Err = OMEGA::IOReadArray(&DvEdgeTmp[0], NEdgesAll, "dvEdge",
                            MeshFileID, EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading dvEdge");
   
   for (int Edge = 0; Edge < NEdgesAll; ++Edge) {
      DvEdgeH(Edge) = DvEdgeTmp[Edge];
   }

   DcEdgeH = ArrayHost1DR8("dcEdge", NEdgesAll);
   std::vector<R8> DcEdgeTmp(NEdgesAll);
   Err = OMEGA::IOReadArray(&DcEdgeTmp[0], NEdgesAll, "dcEdge",
                            MeshFileID, EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading dcEdge");
   
   for (int Edge = 0; Edge < NEdgesAll; ++Edge) {
      DcEdgeH(Edge) = DcEdgeTmp[Edge];
   }

   MeshDensityH = ArrayHost1DR8("meshDensity", NCellsAll);
   std::vector<R8> MeshDensityTmp(NCellsAll);
   Err = OMEGA::IOReadArray(&MeshDensityTmp[0], NCellsAll, "meshDensity",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading meshDensity");
   
   for (int Cell = 0; Cell < NCellsAll; ++Cell) {
      MeshDensityH(Cell) = MeshDensityTmp[Cell];
   }

} // end readMeasurements

void HorzMesh::readWeights() {

} // end readWeights

void HorzMesh::readCoriolis() {

   //FCellH = ArrayHost1DR8("fCell", NCellsAll);
   //std::vector<R8> FCellTmp(NCellsAll);
   //Err = OMEGA::IOReadArray(&FCellTmp[0], NCellsAll, "fCell",
   //                         MeshFileID, CellDecompR8);
   //if (Err != 0)
   //   LOG_CRITICAL("HorzMesh: error reading fCell");
   //
   //for (int Cell = 0; Cell < NCellsAll; ++Cell) {
   //   FCellH(Cell) = FCellTmp[Cell];
   //}

   //FVertexH = ArrayHost1DR8("fVertex", NVerticesAll);
   //std::vector<R8> FVertexTmp(NVerticesAll);
   //Err = OMEGA::IOReadArray(&FVertexTmp[0], NVerticesAll, "fVertex",
   //                         MeshFileID, VertexDecompR8);
   //if (Err != 0)
   //   LOG_CRITICAL("HorzMesh: error reading fVertex");
   //
   //for (int Vertex = 0; Vertex < NVerticesAll; ++Vertex) {
   //   FVertexH(Vertex) = FVertexTmp[Vertex];
   //}

   //FEdgeH = ArrayHost1DR8("fEdge", NEdgesAll);
   //std::vector<R8> FEdgeTmp(NEdgesAll);
   //Err = OMEGA::IOReadArray(&FEdgeTmp[0], NEdgesAll, "fEdge",
   //                         MeshFileID, EdgeDecompR8);
   //if (Err != 0)
   //   LOG_CRITICAL("HorzMesh: error reading fEdge");
   //
   //for (int Edge = 0; Edge < NEdgesAll; ++Edge) {
   //   FEdgeH(Edge) = FEdgeTmp[Edge];
   //}

} // end readCoriolis



} // end namespace OMEGA

//===----------------------------------------------------------------------===//
