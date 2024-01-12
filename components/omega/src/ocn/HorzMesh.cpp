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
   std::vector<I4> CellID(NCellsOwned);
   for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
     CellID[Cell] = MeshDecomp->CellIDH(Cell) - 1;
   }
   
   Err = OMEGA::IOCreateDecomp(CellDecompR8, OMEGA::IOTypeR8, NDims, CellDims,
                               NCellsOwned, CellID, Rearr);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error creating cell IO decomposition");
   
   std::vector<I4> EdgeDims{MeshDecomp->NEdgesGlobal};
   std::vector<I4> EdgeID(NEdgesOwned);
   for (int Edge = 0; Edge < NEdgesOwned; ++Edge) {
     EdgeID[Edge] = MeshDecomp->EdgeIDH(Edge) - 1;
   }
   
   Err = OMEGA::IOCreateDecomp(EdgeDecompR8, OMEGA::IOTypeR8, NDims, EdgeDims,
                               NEdgesOwned, EdgeID, Rearr);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error creating edge IO decomposition");
   
   std::vector<I4> VertexDims{MeshDecomp->NVerticesGlobal};
   std::vector<I4> VertexID(NVerticesOwned);
   for (int Vertex = 0; Vertex < NVerticesOwned; ++Vertex) {
     VertexID[Vertex] = MeshDecomp->VertexIDH(Vertex) - 1;
   }
   
   Err = OMEGA::IOCreateDecomp(VertexDecompR8, OMEGA::IOTypeR8, NDims, VertexDims,
                               NVerticesOwned, VertexID, Rearr);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error creating vertex IO decomposition");
   
   readCoordinates();


} // end constructor

HorzMesh::~HorzMesh() {

   // TODO: add deletes for all arrays and remove from AllDecomps map

} // end deconstructor

void HorzMesh::readCoordinates() {

   I4 Err;
   
   // Read mesh cell coordinates
   XCellH = ArrayHost1DR8("xCell", NCellsOwned);
   std::vector<R8> XCellTmp(NCellsOwned);
   Err = OMEGA::IOReadArray(&XCellTmp[0], NCellsOwned, "xCell",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading xCell");
   
   for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
      XCellH(Cell) = XCellTmp[Cell];
   }
   
   
   YCellH = ArrayHost1DR8("yCell", NCellsOwned);
   std::vector<R8> YCellTmp(NCellsOwned);
   Err = OMEGA::IOReadArray(&YCellTmp[0], NCellsOwned, "yCell",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading yCell");
   
   for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
      YCellH(Cell) = YCellTmp[Cell];
   }
   
   
   ZCellH = ArrayHost1DR8("zCell", NCellsOwned);
   std::vector<R8> ZCellTmp(NCellsOwned);
   Err = OMEGA::IOReadArray(&ZCellTmp[0], NCellsOwned, "zCell",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading zCell");
   
   for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
      ZCellH(Cell) = ZCellTmp[Cell];
   }
   
   
   LonCellH = ArrayHost1DR8("lonCell", NCellsOwned);
   std::vector<R8> LonCellTmp(NCellsOwned);
   Err = OMEGA::IOReadArray(&LonCellTmp[0], NCellsOwned, "lonCell",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading lonCell");
   
   for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
      LonCellH(Cell) = LonCellTmp[Cell];
   }
   
   
   LatCellH = ArrayHost1DR8("latCell", NCellsOwned);
   std::vector<R8> LatCellTmp(NCellsOwned);
   Err = OMEGA::IOReadArray(&LatCellTmp[0], NCellsOwned, "latCell",
                            MeshFileID, CellDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading latCell");
   
   for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
      LatCellH(Cell) = LatCellTmp[Cell];
   }
   
   // Read mesh edge coordinates
   XEdgeH = ArrayHost1DR8("xEdge", NEdgesOwned);
   std::vector<R8> XEdgeTmp(NEdgesOwned);
   Err = OMEGA::IOReadArray(&XEdgeTmp[0], NEdgesOwned, "xEdge",
                            MeshFileID, EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading xEdge");
   
   for (int Edge = 0; Edge < NEdgesOwned; ++Edge) {
      XEdgeH(Edge) = XEdgeTmp[Edge];
   }
   
   
   YEdgeH = ArrayHost1DR8("yEdge", NEdgesOwned);
   std::vector<R8> YEdgeTmp(NEdgesOwned);
   Err = OMEGA::IOReadArray(&YEdgeTmp[0], NEdgesOwned, "yEdge",
                            MeshFileID, EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading yEdge");
   
   for (int Edge = 0; Edge < NEdgesOwned; ++Edge) {
      YEdgeH(Edge) = YEdgeTmp[Edge];
   }
   
   
   ZEdgeH = ArrayHost1DR8("zEdge", NEdgesOwned);
   std::vector<R8> ZEdgeTmp(NEdgesOwned);
   Err = OMEGA::IOReadArray(&ZEdgeTmp[0], NEdgesOwned, "zEdge",
                            MeshFileID, EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading zEdge");
   
   for (int Edge = 0; Edge < NEdgesOwned; ++Edge) {
      ZEdgeH(Edge) = ZEdgeTmp[Edge];
   }
   
   
   LonEdgeH = ArrayHost1DR8("lonEdge", NEdgesOwned);
   std::vector<R8> LonEdgeTmp(NEdgesOwned);
   Err = OMEGA::IOReadArray(&LonEdgeTmp[0], NEdgesOwned, "lonEdge",
                            MeshFileID, EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading lonEdge");
   
   for (int Edge = 0; Edge < NEdgesOwned; ++Edge) {
      LonEdgeH(Edge) = LonEdgeTmp[Edge];
   }
   
   
   LatEdgeH = ArrayHost1DR8("latEdge", NEdgesOwned);
   std::vector<R8> LatEdgeTmp(NEdgesOwned);
   Err = OMEGA::IOReadArray(&LatEdgeTmp[0], NEdgesOwned, "latEdge",
                            MeshFileID, EdgeDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading latEdge");
   
   for (int Edge = 0; Edge < NEdgesOwned; ++Edge) {
      LatEdgeH(Edge) = LatEdgeTmp[Edge];
   }
   
   // Read mesh vertex coordinates
   XVertexH = ArrayHost1DR8("xVertex", NVerticesOwned);
   std::vector<R8> XVertexTmp(NVerticesOwned);
   Err = OMEGA::IOReadArray(&XVertexTmp[0], NVerticesOwned, "xVertex",
                            MeshFileID, VertexDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading xVertex");
   
   for (int Vertex = 0; Vertex < NVerticesOwned; ++Vertex) {
      XVertexH(Vertex) = XVertexTmp[Vertex];
   }
   
   
   YVertexH = ArrayHost1DR8("yVertex", NVerticesOwned);
   std::vector<R8> YVertexTmp(NVerticesOwned);
   Err = OMEGA::IOReadArray(&YVertexTmp[0], NVerticesOwned, "yVertex",
                            MeshFileID, VertexDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading yVertex");
   
   for (int Vertex = 0; Vertex < NVerticesOwned; ++Vertex) {
      YVertexH(Vertex) = YVertexTmp[Vertex];
   }
   
   
   ZVertexH = ArrayHost1DR8("zVertex", NVerticesOwned);
   std::vector<R8> ZVertexTmp(NVerticesOwned);
   Err = OMEGA::IOReadArray(&ZVertexTmp[0], NVerticesOwned, "zVertex",
                            MeshFileID, VertexDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading zVertex");
   
   for (int Vertex = 0; Vertex < NVerticesOwned; ++Vertex) {
      ZVertexH(Vertex) = ZVertexTmp[Vertex];
   }
   
   
   LonVertexH = ArrayHost1DR8("lonVertex", NVerticesOwned);
   std::vector<R8> LonVertexTmp(NVerticesOwned);
   Err = OMEGA::IOReadArray(&LonVertexTmp[0], NVerticesOwned, "lonVertex",
                            MeshFileID, VertexDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading lonVertex");
   
   for (int Vertex = 0; Vertex < NVerticesOwned; ++Vertex) {
      LonVertexH(Vertex) = LonVertexTmp[Vertex];
   }
   
   
   LatVertexH = ArrayHost1DR8("latVertex", NVerticesOwned);
   std::vector<R8> LatVertexTmp(NVerticesOwned);
   Err = OMEGA::IOReadArray(&LatVertexTmp[0], NVerticesOwned, "latVertex",
                            MeshFileID, VertexDecompR8);
   if (Err != 0)
      LOG_CRITICAL("HorzMesh: error reading latVertex");
   
   for (int Vertex = 0; Vertex < NVerticesOwned; ++Vertex) {
      LatVertexH(Vertex) = LatVertexTmp[Vertex];
   }

} // end readCoordinates

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
