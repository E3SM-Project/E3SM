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

//std::cout << "Done with Decomp associations" << std::endl;

// Open the mesh file for reading (assume IO has already been initialized)
I4 FileID;
I4 Err;
Err = OMEGA::IOFileOpen(FileID, MeshFileName, IOModeRead);
if (Err != 0)
   LOG_CRITICAL("HorzMesh: error opening mesh file");

I4 NCellsGlobal;
NCellsGlobal = OMEGA::IOGetDimLength(FileID, "nCells");
std::cout << "Global nCells: " << NCellsGlobal << std::endl;

//Create parallel IO decomposition
I4 CellDecompI4;
I4 CellDecompR8;
I4 NDims = 1;
std::vector<I4> CellDims{MeshDecomp->NCellsGlobal};
std::vector<I4> CellID(NCellsOwned);
IORearranger Rearr = IORearrBox;
for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
  CellID[Cell] = MeshDecomp->CellIDH(Cell) - 1;
  //std::cout << CellID[Cell] << std::endl;  
}

Err = OMEGA::IOCreateDecomp(CellDecompI4, OMEGA::IOTypeI4, NDims, CellDims,
                            NCellsOwned, CellID, Rearr);
Err = OMEGA::IOCreateDecomp(CellDecompR8, OMEGA::IOTypeR8, NDims, CellDims,
                            NCellsOwned, CellID, Rearr);
if (Err != 0)
   LOG_CRITICAL("HorzMesh: error creating cell IO decomposition");


// Read additional mesh variables
ArrayHost1DI4 indexToCellID("indexToCellID", NCellsOwned);
std::vector<I4> indexToCellIDTmp(NCellsOwned);
Err = OMEGA::IOReadArray(&indexToCellIDTmp[0], NCellsOwned, "indexToCellID",
                         FileID, CellDecompI4);
if (Err != 0)
   LOG_CRITICAL("Decomp: error reading xCell");

for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
   indexToCellID(Cell) = indexToCellIDTmp[Cell];
   if ( indexToCellID(Cell)-1 != CellID[Cell] )
      LOG_CRITICAL("Global index read incorrectly");
}

// Read additional mesh variables
XCellH = ArrayHost1DR8("xCell", NCellsOwned);
std::vector<R8> XCellTmp(NCellsOwned);
Err = OMEGA::IOReadArray(&XCellTmp[0], NCellsOwned, "xCell",
                         FileID, CellDecompR8);
if (Err != 0)
   LOG_CRITICAL("HorzMesh: error reading xCell");

for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
   XCellH(Cell) = XCellTmp[Cell];
}


YCellH = ArrayHost1DR8("yCell", NCellsOwned);
std::vector<R8> YCellTmp(NCellsOwned);
Err = OMEGA::IOReadArray(&YCellTmp[0], NCellsOwned, "yCell",
                         FileID, CellDecompR8);
if (Err != 0)
   LOG_CRITICAL("HorzMesh: error reading yCell");

for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
   YCellH(Cell) = YCellTmp[Cell];
}


ZCellH = ArrayHost1DR8("zCell", NCellsOwned);
std::vector<R8> ZCellTmp(NCellsOwned);
Err = OMEGA::IOReadArray(&ZCellTmp[0], NCellsOwned, "zCell",
                         FileID, CellDecompR8);
if (Err != 0)
   LOG_CRITICAL("HorzMesh: error reading zCell");

for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
   ZCellH(Cell) = ZCellTmp[Cell];
}

} // End Constructor

HorzMesh::~HorzMesh() {

   // TODO: add deletes for all arrays and remove from AllDecomps map

} 


} // end namespace OMEGA

//===----------------------------------------------------------------------===//
