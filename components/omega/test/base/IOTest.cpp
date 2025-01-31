//===-- Test driver for OMEGA base IO ----------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA base IO
///
/// This driver tests the OMEGA lower-level IO routines. It writes several
/// distributed arrays and then reads the same file to verify the contents.
///
//
//===-----------------------------------------------------------------------===/

#include "IO.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Logging.h"
#include "MachEnv.h"
#include "mpi.h"

#include <iostream>

using namespace OMEGA;

//------------------------------------------------------------------------------
// The initialization routine for IO testing. It calls various
// init routines, including the creation of the default decomposition.

int initIOTest() {

   int Err = 0;

   // Initialize the Machine Environment class - this also creates
   // the default MachEnv. Then retrieve the default environment and
   // some needed data members.
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // Initialize the Logging system
   initLogging(DefEnv);

   // Open config file
   Config("Omega");
   Err = Config::readAll("omega.yml");
   if (Err != 0) {
      LOG_CRITICAL("IOTest: Error reading config file");
      return Err;
   }

   // Initialize the IO system
   Err = IO::init(DefComm);
   if (Err != 0)
      LOG_ERROR("IOTest: error initializing parallel IO");

   // Create the default decomposition (initializes the decomposition)
   Err = Decomp::init();
   if (Err != 0)
      LOG_ERROR("IOTest: error initializing default decomposition");

   return Err;
}

//------------------------------------------------------------------------------
// The test driver for IO. This creates several distributed arrays and
// associated metadata, writes that data to a file, then re-reads the data
// to ensure both reading and writing produce the same data.
//
int main(int argc, char *argv[]) {

   int RetVal = 0;

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   {
      // Call initialization routine to create the default decomposition
      // and initialize the parallel IO library
      int Err = initIOTest();
      if (Err != 0)
         LOG_CRITICAL("IOTest: Error initializing");

      // Get MPI vars if needed
      MachEnv *DefEnv = MachEnv::getDefault();
      MPI_Comm Comm   = DefEnv->getComm();
      I4 MyTask       = DefEnv->getMyTask();
      I4 NumTasks     = DefEnv->getNumTasks();
      bool IsMaster   = DefEnv->isMasterTask();

      // Retrieve the default decomposition
      Decomp *DefDecomp = Decomp::getDefault();
      if (DefDecomp) { // true if non-null ptr
         LOG_INFO("IOTest: Default decomp retrieval PASS");
      } else {
         LOG_INFO("IOTest: Default decomp retrieval FAIL");
         return -1;
      }

      // Create Kokkos arrays of each type and at various mesh locations
      I4 NCellsSize      = DefDecomp->NCellsSize;
      I4 NEdgesSize      = DefDecomp->NEdgesSize;
      I4 NVerticesSize   = DefDecomp->NVerticesSize;
      I4 NCellsOwned     = DefDecomp->NCellsOwned;
      I4 NEdgesOwned     = DefDecomp->NEdgesOwned;
      I4 NVerticesOwned  = DefDecomp->NVerticesOwned;
      I4 NCellsGlobal    = DefDecomp->NCellsGlobal;
      I4 NEdgesGlobal    = DefDecomp->NEdgesGlobal;
      I4 NVerticesGlobal = DefDecomp->NVerticesGlobal;
      I4 NVertLevels     = 128;

      HostArray1DI4 RefI4Vert("RefI4Vert", NVertLevels);
      HostArray1DI8 RefI8Vert("RefI8Vert", NVertLevels);
      HostArray1DR4 RefR4Vert("RefR4Vert", NVertLevels);
      HostArray1DR8 RefR8Vert("RefR8Vert", NVertLevels);

      HostArray2DI4 RefI4Cell("RefI4Cell", NCellsSize, NVertLevels);
      HostArray2DI8 RefI8Cell("RefI8Cell", NCellsSize, NVertLevels);
      HostArray2DR4 RefR4Cell("RefR4Cell", NCellsSize, NVertLevels);
      HostArray2DR8 RefR8Cell("RefR8Cell", NCellsSize, NVertLevels);

      HostArray2DI4 RefI4Edge("RefI4Edge", NEdgesSize, NVertLevels);
      HostArray2DI8 RefI8Edge("RefI8Edge", NEdgesSize, NVertLevels);
      HostArray2DR4 RefR4Edge("RefR4Edge", NEdgesSize, NVertLevels);
      HostArray2DR8 RefR8Edge("RefR8Edge", NEdgesSize, NVertLevels);

      HostArray2DI4 RefI4Vrtx("RefI4Vrtx", NVerticesSize, NVertLevels);
      HostArray2DI8 RefI8Vrtx("RefI8Vrtx", NVerticesSize, NVertLevels);
      HostArray2DR4 RefR4Vrtx("RefR4Vrtx", NVerticesSize, NVertLevels);
      HostArray2DR8 RefR8Vrtx("RefR8Vrtx", NVerticesSize, NVertLevels);

      HostArray1DI4 CellIDH = DefDecomp->CellIDH;
      HostArray1DI4 EdgeIDH = DefDecomp->EdgeIDH;
      HostArray1DI4 VrtxIDH = DefDecomp->VertexIDH;

      // Create local non-distributed arrays and scalars
      I4 RefI4Scalar = -1;
      I8 RefI8Scalar = -2;
      R4 RefR4Scalar = -3.1;
      R8 RefR8Scalar = -4.56789;
      for (int K = 0; K < NVertLevels; ++K) {
         RefI4Vert(K) = K;
         RefI8Vert(K) = K * 2;
         RefR4Vert(K) = K * 3.1;
         RefR8Vert(K) = K * 4.1234567;
      }

      // Offset arrays - initialize to -1, corresponding to entries
      // that should not be written;
      std::vector<int> OffsetCell(NCellsSize * NVertLevels, -1);
      std::vector<int> OffsetEdge(NEdgesSize * NVertLevels, -1);
      std::vector<int> OffsetVrtx(NVerticesSize * NVertLevels, -1);
      for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
         int GlobalCellAdd = CellIDH(Cell) - 1; // 0-based offset
         for (int k = 0; k < NVertLevels; ++k) {
            RefI4Cell(Cell, k)    = GlobalCellAdd * k;
            RefI8Cell(Cell, k)    = GlobalCellAdd * k * 1000000000;
            RefR4Cell(Cell, k)    = GlobalCellAdd * k * 123.45;
            RefR8Cell(Cell, k)    = GlobalCellAdd * k * 1.23456789;
            int VectorAdd         = Cell * NVertLevels + k;
            OffsetCell[VectorAdd] = GlobalCellAdd * NVertLevels + k;
         }
      }

      for (int Edge = 0; Edge < NEdgesOwned; ++Edge) {
         int GlobalEdgeAdd = EdgeIDH(Edge) - 1;
         for (int k = 0; k < NVertLevels; ++k) {
            RefI4Edge(Edge, k)    = GlobalEdgeAdd * k;
            RefI8Edge(Edge, k)    = GlobalEdgeAdd * k * 1000000000;
            RefR4Edge(Edge, k)    = GlobalEdgeAdd * k * 123.45;
            RefR8Edge(Edge, k)    = GlobalEdgeAdd * k * 1.23456789;
            int VectorAdd         = Edge * NVertLevels + k;
            OffsetEdge[VectorAdd] = GlobalEdgeAdd * NVertLevels + k;
         }
      }

      for (int Vrtx = 0; Vrtx < NVerticesOwned; ++Vrtx) {
         int GlobalVrtxAdd = VrtxIDH(Vrtx) - 1;
         for (int k = 0; k < NVertLevels; ++k) {
            RefI4Vrtx(Vrtx, k)    = GlobalVrtxAdd * k;
            RefI8Vrtx(Vrtx, k)    = GlobalVrtxAdd * k * 1000000000;
            RefR4Vrtx(Vrtx, k)    = GlobalVrtxAdd * k * 123.45;
            RefR8Vrtx(Vrtx, k)    = GlobalVrtxAdd * k * 1.23456789;
            int VectorAdd         = Vrtx * NVertLevels + k;
            OffsetVrtx[VectorAdd] = GlobalVrtxAdd * NVertLevels + k;
         }
      }

      // Create the needed decomposition offsets

      int DecompCellI4;
      int DecompCellI8;
      int DecompCellR4;
      int DecompCellR8;
      int DecompEdgeI4;
      int DecompEdgeI8;
      int DecompEdgeR4;
      int DecompEdgeR8;
      int DecompVrtxI4;
      int DecompVrtxI8;
      int DecompVrtxR4;
      int DecompVrtxR8;
      std::vector<int> CellDims{NCellsGlobal, NVertLevels};
      std::vector<int> EdgeDims{NEdgesGlobal, NVertLevels};
      std::vector<int> VrtxDims{NVerticesGlobal, NVertLevels};
      int CellArraySize = NCellsSize * NVertLevels;
      int EdgeArraySize = NEdgesSize * NVertLevels;
      int VrtxArraySize = NVerticesSize * NVertLevels;

      Err = IO::createDecomp(DecompCellI4, IO::IOTypeI4, 2, CellDims,
                             CellArraySize, OffsetCell, IO::DefaultRearr);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error creating cell decomp I4 FAIL");
      }
      Err = IO::createDecomp(DecompCellI8, IO::IOTypeI8, 2, CellDims,
                             CellArraySize, OffsetCell, IO::DefaultRearr);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error creating cell decomp I8 FAIL");
      }
      Err = IO::createDecomp(DecompCellR4, IO::IOTypeR4, 2, CellDims,
                             CellArraySize, OffsetCell, IO::DefaultRearr);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error creating cell decomp R4 FAIL");
      }
      Err = IO::createDecomp(DecompCellR8, IO::IOTypeR8, 2, CellDims,
                             CellArraySize, OffsetCell, IO::DefaultRearr);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error creating cell decomp R8 FAIL");
      }
      Err = IO::createDecomp(DecompEdgeI4, IO::IOTypeI4, 2, EdgeDims,
                             EdgeArraySize, OffsetEdge, IO::DefaultRearr);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error creating edge decomp I4 FAIL");
      }
      Err = IO::createDecomp(DecompEdgeI8, IO::IOTypeI8, 2, EdgeDims,
                             EdgeArraySize, OffsetEdge, IO::DefaultRearr);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error creating edge decomp I8 FAIL");
      }
      Err = IO::createDecomp(DecompEdgeR4, IO::IOTypeR4, 2, EdgeDims,
                             EdgeArraySize, OffsetEdge, IO::DefaultRearr);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error creating edge decomp R4 FAIL");
      }
      Err = IO::createDecomp(DecompEdgeR8, IO::IOTypeR8, 2, EdgeDims,
                             EdgeArraySize, OffsetEdge, IO::DefaultRearr);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error creating edge decomp R8 FAIL");
      }
      Err = IO::createDecomp(DecompVrtxI4, IO::IOTypeI4, 2, VrtxDims,
                             VrtxArraySize, OffsetVrtx, IO::DefaultRearr);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error creating vertex decomp I4 FAIL");
      }
      Err = IO::createDecomp(DecompVrtxI8, IO::IOTypeI8, 2, VrtxDims,
                             VrtxArraySize, OffsetVrtx, IO::DefaultRearr);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error creating vertex decomp I8 FAIL");
      }
      Err = IO::createDecomp(DecompVrtxR4, IO::IOTypeR4, 2, VrtxDims,
                             VrtxArraySize, OffsetVrtx, IO::DefaultRearr);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error creating vertex decomp R4 FAIL");
      }
      Err = IO::createDecomp(DecompVrtxR8, IO::IOTypeR8, 2, VrtxDims,
                             VrtxArraySize, OffsetVrtx, IO::DefaultRearr);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error creating vertex decomp R8 FAIL");
      }

      // Open a file for output
      int OutFileID;
      Err = IO::openFile(OutFileID, "IOTest.nc", IO::ModeWrite, IO::FmtDefault,
                         IO::IfExists::Replace);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error opening file for output FAIL");
      }
      // Define array dimensions
      int DimCellID;
      int DimEdgeID;
      int DimVrtxID;
      int DimVertID;
      Err = IO::defineDim(OutFileID, "NVertLevels", NVertLevels, DimVertID);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error defining vertical dimension FAIL");
      }
      Err = IO::defineDim(OutFileID, "NCells", NCellsGlobal, DimCellID);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error defining Cell dimension FAIL");
      }
      Err = IO::defineDim(OutFileID, "NEdges", NEdgesGlobal, DimEdgeID);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error defining Edge dimension FAIL");
      }
      Err = IO::defineDim(OutFileID, "NVertices", NVerticesGlobal, DimVrtxID);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error defining Vertex dimension FAIL");
      }

      // Write some global file metadata
      I4 FileMetaI4Ref          = 2;
      I8 FileMetaI8Ref          = 4;
      R4 FileMetaR4Ref          = 6.789;
      R8 FileMetaR8Ref          = 1.23456789;
      std::string FileMetaDescr = "OMEGA IO Unit test file";

      Err = IO::writeMeta("FileMetaI4", FileMetaI4Ref, OutFileID, IO::GlobalID);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing global I4 metadata FAIL");
      }
      Err = IO::writeMeta("FileMetaI8", FileMetaI8Ref, OutFileID, IO::GlobalID);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing global I8 metadata FAIL");
      }
      Err = IO::writeMeta("FileMetaR4", FileMetaR4Ref, OutFileID, IO::GlobalID);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing global R4 metadata FAIL");
      }
      Err = IO::writeMeta("FileMetaR8", FileMetaR8Ref, OutFileID, IO::GlobalID);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing global R8 metadata FAIL");
      }
      Err = IO::writeMeta("FileMetaDescr", FileMetaDescr, OutFileID,
                          IO::GlobalID);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing global char metadata FAIL");
      }
      Err = IO::writeMeta("StringLiteral", "MyString", OutFileID, IO::GlobalID);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing char literal metadata FAIL");
      }

      // Define variables/arrays
      int VarIDScalarI4;
      int VarIDScalarI8;
      int VarIDScalarR4;
      int VarIDScalarR8;
      int VarIDI4Vert;
      int VarIDI8Vert;
      int VarIDR4Vert;
      int VarIDR8Vert;
      int VarIDCellI4;
      int VarIDCellI8;
      int VarIDCellR4;
      int VarIDCellR8;
      int VarIDEdgeI4;
      int VarIDEdgeI8;
      int VarIDEdgeR4;
      int VarIDEdgeR8;
      int VarIDVrtxI4;
      int VarIDVrtxI8;
      int VarIDVrtxR4;
      int VarIDVrtxR8;

      int VertDimIDs[1] = {DimVertID};
      int CellDimIDs[2] = {DimCellID, DimVertID};
      int EdgeDimIDs[2] = {DimEdgeID, DimVertID};
      int VrtxDimIDs[2] = {DimVrtxID, DimVertID};

      Err = IO::defineVar(OutFileID, "ScalarI4", IO::IOTypeI4, 0, nullptr,
                          VarIDScalarI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: Error defining I4 scalar FAIL");
      }
      Err = IO::defineVar(OutFileID, "ScalarI8", IO::IOTypeI8, 0, nullptr,
                          VarIDScalarI8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: Error defining I8 scalar FAIL");
      }
      Err = IO::defineVar(OutFileID, "ScalarR4", IO::IOTypeR4, 0, nullptr,
                          VarIDScalarR4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: Error defining R4 scalar FAIL");
      }
      Err = IO::defineVar(OutFileID, "ScalarR8", IO::IOTypeR8, 0, nullptr,
                          VarIDScalarR8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: Error defining R8 scalar FAIL");
      }
      Err = IO::defineVar(OutFileID, "I4Vert", IO::IOTypeI4, 1, VertDimIDs,
                          VarIDI4Vert);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: Error defining I4Vert array FAIL");
      }
      Err = IO::defineVar(OutFileID, "I8Vert", IO::IOTypeI8, 1, VertDimIDs,
                          VarIDI8Vert);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: Error defining I8Vert array FAIL");
      }
      Err = IO::defineVar(OutFileID, "R4Vert", IO::IOTypeR4, 1, VertDimIDs,
                          VarIDR4Vert);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: Error defining R4Vert array FAIL");
      }
      Err = IO::defineVar(OutFileID, "R8Vert", IO::IOTypeR8, 1, VertDimIDs,
                          VarIDR8Vert);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: Error defining R8Vert array FAIL");
      }
      Err = IO::defineVar(OutFileID, "CellI4", IO::IOTypeI4, 2, CellDimIDs,
                          VarIDCellI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: Error defining CellI4 array FAIL");
      }
      Err = IO::defineVar(OutFileID, "CellI8", IO::IOTypeI8, 2, CellDimIDs,
                          VarIDCellI8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: Error defining CellI8 array FAIL");
      }
      Err = IO::defineVar(OutFileID, "CellR4", IO::IOTypeR4, 2, CellDimIDs,
                          VarIDCellR4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: Error defining CellR4 array FAIL");
      }
      Err = IO::defineVar(OutFileID, "CellR8", IO::IOTypeR8, 2, CellDimIDs,
                          VarIDCellR8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: Error defining CellR8 array FAIL");
      }

      Err = IO::defineVar(OutFileID, "EdgeI4", IO::IOTypeI4, 2, EdgeDimIDs,
                          VarIDEdgeI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: Error defining EdgeI4 array FAIL");
      }
      Err = IO::defineVar(OutFileID, "EdgeI8", IO::IOTypeI8, 2, EdgeDimIDs,
                          VarIDEdgeI8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: Error defining EdgeI8 array FAIL");
      }
      Err = IO::defineVar(OutFileID, "EdgeR4", IO::IOTypeR4, 2, EdgeDimIDs,
                          VarIDEdgeR4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: Error defining EdgeR4 array FAIL");
      }
      Err = IO::defineVar(OutFileID, "EdgeR8", IO::IOTypeR8, 2, EdgeDimIDs,
                          VarIDEdgeR8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: Error defining EdgeR8 array FAIL");
      }

      Err = IO::defineVar(OutFileID, "VrtxI4", IO::IOTypeI4, 2, VrtxDimIDs,
                          VarIDVrtxI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: Error defining VrtxI4 array FAIL");
      }
      Err = IO::defineVar(OutFileID, "VrtxI8", IO::IOTypeI8, 2, VrtxDimIDs,
                          VarIDVrtxI8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: Error defining VrtxI8 array FAIL");
      }
      Err = IO::defineVar(OutFileID, "VrtxR4", IO::IOTypeR4, 2, VrtxDimIDs,
                          VarIDVrtxR4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: Error defining VrtxR4 array FAIL");
      }
      Err = IO::defineVar(OutFileID, "VrtxR8", IO::IOTypeR8, 2, VrtxDimIDs,
                          VarIDVrtxR8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: Error defining VrtxR8 array FAIL");
      }

      // Add Variable metadata just for one array
      I4 VarMetaI4Ref             = 3;
      I8 VarMetaI8Ref             = 5;
      R4 VarMetaR4Ref             = 5.789;
      R8 VarMetaR8Ref             = 2.23456789;
      std::string VarMetaDescrRef = "Test array for I4 on Cells";

      Err = IO::writeMeta("VarMetaI4", VarMetaI4Ref, OutFileID, VarIDCellI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing var I4 metadata FAIL");
      }
      Err = IO::writeMeta("VarMetaI8", VarMetaI8Ref, OutFileID, VarIDCellI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing var I8 metadata FAIL");
      }
      Err = IO::writeMeta("VarMetaR4", VarMetaR4Ref, OutFileID, VarIDCellI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing var R4 metadata FAIL");
      }
      Err = IO::writeMeta("VarMetaR8", VarMetaR8Ref, OutFileID, VarIDCellI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing var R8 metadata FAIL");
      }
      Err = IO::writeMeta("VarMetaDescr", VarMetaDescrRef, OutFileID,
                          VarIDCellI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing var char metadata FAIL");
      }

      // Exit define mode
      Err = IO::endDefinePhase(OutFileID);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error ending define mode FAIL");
      }

      // Write variables
      I4 FillI4 = -999;
      I8 FillI8 = -999999;
      R4 FillR4 = -1.234e30;
      R8 FillR8 = -1.23456789e30;

      // Write non-distributed variables
      Err = IO::writeNDVar(&RefI4Scalar, OutFileID, VarIDScalarI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing I4 scalar FAIL");
      }

      Err = IO::writeNDVar(&RefI8Scalar, OutFileID, VarIDScalarI8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing I8 scalar FAIL");
      }

      Err = IO::writeNDVar(&RefR4Scalar, OutFileID, VarIDScalarR4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing R4 scalar FAIL");
      }

      Err = IO::writeNDVar(&RefR8Scalar, OutFileID, VarIDScalarR8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing R8 scalar FAIL");
      }

      Err = IO::writeNDVar(RefI4Vert.data(), OutFileID, VarIDI4Vert);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing I4Vert vector FAIL");
      }

      Err = IO::writeNDVar(RefI8Vert.data(), OutFileID, VarIDI8Vert);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing I8Vert vector FAIL");
      }

      Err = IO::writeNDVar(RefR4Vert.data(), OutFileID, VarIDR4Vert);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing R4Vert vector FAIL");
      }

      Err = IO::writeNDVar(RefR8Vert.data(), OutFileID, VarIDR8Vert);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing R8Vert vector FAIL");
      }

      // Write distributed arrays
      Err = IO::writeArray(RefI4Cell.data(), NCellsSize * NVertLevels, &FillI4,
                           OutFileID, DecompCellI4, VarIDCellI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing I4 array on cells FAIL");
      }
      Err = IO::writeArray(RefI8Cell.data(), NCellsSize * NVertLevels, &FillI8,
                           OutFileID, DecompCellI8, VarIDCellI8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing I8 array on cells FAIL");
      }
      Err = IO::writeArray(RefR4Cell.data(), NCellsSize * NVertLevels, &FillR4,
                           OutFileID, DecompCellR4, VarIDCellR4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing R4 array on cells FAIL");
      }
      Err = IO::writeArray(RefR8Cell.data(), NCellsSize * NVertLevels, &FillR8,
                           OutFileID, DecompCellR8, VarIDCellR8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing R8 array on cells FAIL");
      }

      Err = IO::writeArray(RefI4Edge.data(), NEdgesSize * NVertLevels, &FillI4,
                           OutFileID, DecompEdgeI4, VarIDEdgeI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing I4 array on Edges FAIL");
      }
      Err = IO::writeArray(RefI8Edge.data(), NEdgesSize * NVertLevels, &FillI8,
                           OutFileID, DecompEdgeI8, VarIDEdgeI8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing I8 array on Edges FAIL");
      }
      Err = IO::writeArray(RefR4Edge.data(), NEdgesSize * NVertLevels, &FillR4,
                           OutFileID, DecompEdgeR4, VarIDEdgeR4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing R4 array on Edges FAIL");
      }
      Err = IO::writeArray(RefR8Edge.data(), NEdgesSize * NVertLevels, &FillR8,
                           OutFileID, DecompEdgeR8, VarIDEdgeR8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing R8 array on Edges FAIL");
      }

      Err = IO::writeArray(RefI4Vrtx.data(), NVerticesSize * NVertLevels,
                           &FillI4, OutFileID, DecompVrtxI4, VarIDVrtxI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing I4 array on vertices FAIL");
      }
      Err = IO::writeArray(RefI8Vrtx.data(), NVerticesSize * NVertLevels,
                           &FillI8, OutFileID, DecompVrtxI8, VarIDVrtxI8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing I8 array on vertices FAIL");
      }
      Err = IO::writeArray(RefR4Vrtx.data(), NVerticesSize * NVertLevels,
                           &FillR4, OutFileID, DecompVrtxR4, VarIDVrtxR4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing R4 array on vertices FAIL");
      }
      Err = IO::writeArray(RefR8Vrtx.data(), NVerticesSize * NVertLevels,
                           &FillR8, OutFileID, DecompVrtxR8, VarIDVrtxR8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error writing R8 array on vertices FAIL");
      }

      // Finished writing, close file
      Err = IO::closeFile(OutFileID);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error closing output file FAIL");
      }

      // Open a file for reading to verify read/write
      int InFileID;
      Err = IO::openFile(InFileID, "IOTest.nc", IO::ModeRead);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error opening file for reading FAIL");
      }

      // Get dimension lengths to verify read/write of dimension info
      I4 NVertLevelsID;
      I4 NVertLevelsNew;
      Err = IO::getDimFromFile(InFileID, "NVertLevels", NVertLevelsID,
                               NVertLevelsNew);
      if (Err == 0 and NVertLevelsNew == NVertLevels) {
         LOG_INFO("IOTest: read/write vert dimension test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write vert dimension test FAIL");
      }

      I4 NCellsNewID;
      I4 NCellsNew;
      Err = IO::getDimFromFile(InFileID, "NCells", NCellsNewID, NCellsNew);
      if (Err == 0 and NCellsNew == NCellsGlobal) {
         LOG_INFO("IOTest: read/write cell dimension test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write cell dimension test FAIL");
      }

      I4 NEdgesNewID;
      I4 NEdgesNew;
      Err = IO::getDimFromFile(InFileID, "NEdges", NEdgesNewID, NEdgesNew);
      if (Err == 0 and NEdgesNew == NEdgesGlobal) {
         LOG_INFO("IOTest: read/write edge dimension test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write edge dimension test FAIL");
      }

      I4 NVerticesNewID;
      I4 NVerticesNew;
      Err = IO::getDimFromFile(InFileID, "NVertices", NVerticesNewID,
                               NVerticesNew);
      if (Err == 0 and NVerticesNew == NVerticesGlobal) {
         LOG_INFO("IOTest: read/write vertex dimension test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write vertex dimension test FAIL");
      }

      // Read global attributes
      I4 FileMetaI4New;
      I8 FileMetaI8New;
      R4 FileMetaR4New;
      R8 FileMetaR8New;
      std::string FileMetaDescrNew;

      Err = IO::readMeta("FileMetaI4", FileMetaI4New, InFileID, IO::GlobalID);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading file I4 metadata FAIL");
      }
      if (FileMetaI4New == FileMetaI4Ref) {
         LOG_INFO("IOTest: read/write file metadata I4 test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write file metadata I4 test FAIL");
      }

      Err = IO::readMeta("FileMetaI8", FileMetaI8New, InFileID, IO::GlobalID);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading file I8 metadata FAIL");
      }
      if (FileMetaI8New == FileMetaI8Ref) {
         LOG_INFO("IOTest: read/write file metadata I8 test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write file metadata I8 test FAIL");
      }

      Err = IO::readMeta("FileMetaR4", FileMetaR4New, InFileID, IO::GlobalID);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading file R4 metadata FAIL");
      }
      if (FileMetaR4New == FileMetaR4Ref) {
         LOG_INFO("IOTest: read/write file metadata R4 test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write file metadata R4 test FAIL");
      }

      Err = IO::readMeta("FileMetaR8", FileMetaR8New, InFileID, IO::GlobalID);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading file R8 metadata FAIL");
      }
      if (FileMetaR8New == FileMetaR8Ref) {
         LOG_INFO("IOTest: read/write file metadata R8 test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write file metadata R8 test FAIL");
      }

      Err = IO::readMeta("FileMetaDescr", FileMetaDescrNew, InFileID,
                         IO::GlobalID);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading file string metadata FAIL");
      }
      if (FileMetaDescrNew == FileMetaDescr) {
         LOG_INFO("IOTest: read/write file metadata string test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write file metadata string test FAIL");
      }

      std::string MyStringNew;
      Err = IO::readMeta("StringLiteral", MyStringNew, InFileID, IO::GlobalID);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading file string literal FAIL");
      }
      if (MyStringNew == "MyString") {
         LOG_INFO("IOTest: read/write file metadata string literal test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write file metadata string literal test FAIL");
      }
      // Read new variables

      I4 NewI4Scalar = 0;
      I8 NewI8Scalar = 0;
      R4 NewR4Scalar = 0;
      R8 NewR8Scalar = 0;

      HostArray1DI4 NewI4Vert("NewI4Vert", NVertLevels);
      HostArray1DI8 NewI8Vert("NewI8Vert", NVertLevels);
      HostArray1DR4 NewR4Vert("NewR4Vert", NVertLevels);
      HostArray1DR8 NewR8Vert("NewR8Vert", NVertLevels);

      HostArray2DI4 NewI4Cell("NewI4Cell", NCellsSize, NVertLevels);
      HostArray2DI8 NewI8Cell("NewI8Cell", NCellsSize, NVertLevels);
      HostArray2DR4 NewR4Cell("NewR4Cell", NCellsSize, NVertLevels);
      HostArray2DR8 NewR8Cell("NewR8Cell", NCellsSize, NVertLevels);

      HostArray2DI4 NewI4Edge("NewI4Edge", NEdgesSize, NVertLevels);
      HostArray2DI8 NewI8Edge("NewI8Edge", NEdgesSize, NVertLevels);
      HostArray2DR4 NewR4Edge("NewR4Edge", NEdgesSize, NVertLevels);
      HostArray2DR8 NewR8Edge("NewR8Edge", NEdgesSize, NVertLevels);

      HostArray2DI4 NewI4Vrtx("NewI4Vrtx", NVerticesSize, NVertLevels);
      HostArray2DI8 NewI8Vrtx("NewI8Vrtx", NVerticesSize, NVertLevels);
      HostArray2DR4 NewR4Vrtx("NewR4Vrtx", NVerticesSize, NVertLevels);
      HostArray2DR8 NewR8Vrtx("NewR8Vrtx", NVerticesSize, NVertLevels);

      // Read non-distributed variables
      Err = IO::readNDVar(&NewI4Scalar, "ScalarI4", InFileID, VarIDScalarI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading I4 scalar FAIL");
      }

      Err = IO::readNDVar(&NewI8Scalar, "ScalarI8", InFileID, VarIDScalarI8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading I8 scalar FAIL");
      }

      Err = IO::readNDVar(&NewR4Scalar, "ScalarR4", InFileID, VarIDScalarR4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading R4 scalar FAIL");
      }

      Err = IO::readNDVar(&NewR8Scalar, "ScalarR8", InFileID, VarIDScalarR8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading R8 scalar FAIL");
      }

      Err = IO::readNDVar(NewI4Vert.data(), "I4Vert", InFileID, VarIDI4Vert);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading I4Vert vector FAIL");
      }

      Err = IO::readNDVar(NewI8Vert.data(), "I8Vert", InFileID, VarIDI8Vert);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading I8Vert vector FAIL");
      }

      Err = IO::readNDVar(NewR4Vert.data(), "R4Vert", InFileID, VarIDR4Vert);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading R4Vert vector FAIL");
      }

      Err = IO::readNDVar(NewR8Vert.data(), "R8Vert", InFileID, VarIDR8Vert);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading R8Vert vector FAIL");
      }

      // Read distributed arrays
      Err = IO::readArray(NewI4Cell.data(), NCellsSize * NVertLevels, "CellI4",
                          InFileID, DecompCellI4, VarIDCellI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading I4 array on cells FAIL");
      }
      Err = IO::readArray(NewI8Cell.data(), NCellsSize * NVertLevels, "CellI8",
                          InFileID, DecompCellI8, VarIDCellI8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading I8 array on cells FAIL");
      }
      Err = IO::readArray(NewR4Cell.data(), NCellsSize * NVertLevels, "CellR4",
                          InFileID, DecompCellR4, VarIDCellR4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading R4 array on cells FAIL");
      }
      Err = IO::readArray(NewR8Cell.data(), NCellsSize * NVertLevels, "CellR8",
                          InFileID, DecompCellR8, VarIDCellR8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading R8 array on cells FAIL");
      }

      Err = IO::readArray(NewI4Edge.data(), NEdgesSize * NVertLevels, "EdgeI4",
                          InFileID, DecompEdgeI4, VarIDEdgeI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading I4 array on Edges FAIL");
      }
      Err = IO::readArray(NewI8Edge.data(), NEdgesSize * NVertLevels, "EdgeI8",
                          InFileID, DecompEdgeI8, VarIDEdgeI8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading I8 array on Edges FAIL");
      }
      Err = IO::readArray(NewR4Edge.data(), NEdgesSize * NVertLevels, "EdgeR4",
                          InFileID, DecompEdgeR4, VarIDEdgeR4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading R4 array on Edges FAIL");
      }
      Err = IO::readArray(NewR8Edge.data(), NEdgesSize * NVertLevels, "EdgeR8",
                          InFileID, DecompEdgeR8, VarIDEdgeR8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading R8 array on Edges FAIL");
      }

      Err = IO::readArray(NewI4Vrtx.data(), NVerticesSize * NVertLevels,
                          "VrtxI4", InFileID, DecompVrtxI4, VarIDVrtxI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading I4 array on vertices FAIL");
      }
      Err = IO::readArray(NewI8Vrtx.data(), NVerticesSize * NVertLevels,
                          "VrtxI8", InFileID, DecompVrtxI8, VarIDVrtxI8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading I8 array on vertices FAIL");
      }
      Err = IO::readArray(NewR4Vrtx.data(), NVerticesSize * NVertLevels,
                          "VrtxR4", InFileID, DecompVrtxR4, VarIDVrtxR4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading R4 array on vertices FAIL");
      }
      Err = IO::readArray(NewR8Vrtx.data(), NVerticesSize * NVertLevels,
                          "VrtxR8", InFileID, DecompVrtxR8, VarIDVrtxR8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading R8 array on vertices FAIL");
      }

      // Check that variables match the reference cases that were written
      // Only check the owned values for distributed arrays - these would need
      // to be followed by a halo update.

      int Err1 = 0;
      int Err2 = 0;
      int Err3 = 0;
      int Err4 = 0;

      if (NewI4Scalar != RefI4Scalar)
         Err1++;
      if (NewI8Scalar != RefI8Scalar)
         Err2++;
      if (NewR4Scalar != RefR4Scalar)
         Err3++;
      if (NewR8Scalar != RefR8Scalar)
         Err4++;
      if (Err1 == 0) {
         LOG_INFO("IOTest: read/write scalar I4 test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write scalar I4 test FAIL");
      }
      if (Err2 == 0) {
         LOG_INFO("IOTest: read/write scalar I8 test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write scalar I8 test FAIL");
      }
      if (Err3 == 0) {
         LOG_INFO("IOTest: read/write scalar R4 test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write scalar R4 test FAIL");
      }
      if (Err4 == 0) {
         LOG_INFO("IOTest: read/write scalar R8 test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write scalar R8 test FAIL");
      }

      Err1 = 0;
      Err2 = 0;
      Err3 = 0;
      Err4 = 0;
      for (int k = 0; k < NVertLevels; ++k) {
         if (NewI4Vert(k) != RefI4Vert(k))
            Err1++;
         if (NewI8Vert(k) != RefI8Vert(k))
            Err2++;
         if (NewR4Vert(k) != RefR4Vert(k))
            Err3++;
         if (NewR8Vert(k) != RefR8Vert(k))
            Err4++;
      }
      if (Err1 == 0) {
         LOG_INFO("IOTest: read/write vert I4 vector test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write vert I4 vector test FAIL");
      }
      if (Err2 == 0) {
         LOG_INFO("IOTest: read/write vert I8 vector test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write vert I8 vector test FAIL");
      }
      if (Err3 == 0) {
         LOG_INFO("IOTest: read/write vert R4 vector test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write vert R4 vector test FAIL");
      }
      if (Err4 == 0) {
         LOG_INFO("IOTest: read/write vert R8 vector test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write vert R8 vector test FAIL");
      }

      Err1 = 0;
      Err2 = 0;
      Err3 = 0;
      Err4 = 0;
      for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
         for (int k = 0; k < NVertLevels; ++k) {
            if (NewI4Cell(Cell, k) != RefI4Cell(Cell, k))
               Err1++;
            if (NewI8Cell(Cell, k) != RefI8Cell(Cell, k))
               Err2++;
            if (NewR4Cell(Cell, k) != RefR4Cell(Cell, k))
               Err3++;
            if (NewR8Cell(Cell, k) != RefR8Cell(Cell, k))
               Err4++;
         }
      }
      if (Err1 == 0) {
         LOG_INFO("IOTest: read/write array I4 on Cells test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write array I4 on Cells test FAIL");
      }
      if (Err2 == 0) {
         LOG_INFO("IOTest: read/write array I8 on Cells test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write array I8 on Cells test FAIL");
      }
      if (Err3 == 0) {
         LOG_INFO("IOTest: read/write array R4 on Cells test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write array R4 on Cells test FAIL");
      }
      if (Err4 == 0) {
         LOG_INFO("IOTest: read/write array R8 on Cells test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write array R8 on Cells test FAIL");
      }

      Err1 = 0;
      Err2 = 0;
      Err3 = 0;
      Err4 = 0;
      for (int Edge = 0; Edge < NEdgesOwned; ++Edge) {
         for (int k = 0; k < NVertLevels; ++k) {
            if (NewI4Edge(Edge, k) != RefI4Edge(Edge, k))
               Err1++;
            if (NewI8Edge(Edge, k) != RefI8Edge(Edge, k))
               Err2++;
            if (NewR4Edge(Edge, k) != RefR4Edge(Edge, k))
               Err3++;
            if (NewR8Edge(Edge, k) != RefR8Edge(Edge, k))
               Err4++;
         }
      }
      if (Err1 == 0) {
         LOG_INFO("IOTest: read/write array I4 on Edges test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write array I4 on Edges test FAIL");
      }
      if (Err2 == 0) {
         LOG_INFO("IOTest: read/write array I8 on Edges test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write array I8 on Edges test FAIL");
      }
      if (Err3 == 0) {
         LOG_INFO("IOTest: read/write array R4 on Edges test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write array R4 on Edges test FAIL");
      }
      if (Err4 == 0) {
         LOG_INFO("IOTest: read/write array R8 on Edges test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write array R8 on Edges test FAIL");
      }

      Err1 = 0;
      Err2 = 0;
      Err3 = 0;
      Err4 = 0;
      for (int Vrtx = 0; Vrtx < NVerticesOwned; ++Vrtx) {
         for (int k = 0; k < NVertLevels; ++k) {
            if (NewI4Vrtx(Vrtx, k) != RefI4Vrtx(Vrtx, k))
               Err1++;
            if (NewI8Vrtx(Vrtx, k) != RefI8Vrtx(Vrtx, k))
               Err2++;
            if (NewR4Vrtx(Vrtx, k) != RefR4Vrtx(Vrtx, k))
               Err3++;
            if (NewR8Vrtx(Vrtx, k) != RefR8Vrtx(Vrtx, k))
               Err4++;
         }
      }
      if (Err1 == 0) {
         LOG_INFO("IOTest: read/write array I4 on Vertices test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write array I4 on Vertices test FAIL");
      }
      if (Err2 == 0) {
         LOG_INFO("IOTest: read/write array I8 on Vertices test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write array I8 on Vertices test FAIL");
      }
      if (Err3 == 0) {
         LOG_INFO("IOTest: read/write array R4 on Vertices test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write array R4 on Vertices test FAIL");
      }
      if (Err4 == 0) {
         LOG_INFO("IOTest: read/write array R8 on Vertices test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write array R8 on Vertices test FAIL");
      }

      // Read array attributes
      I4 VarMetaI4New;
      I8 VarMetaI8New;
      R4 VarMetaR4New;
      R8 VarMetaR8New;
      std::string VarMetaDescrNew;

      Err = IO::readMeta("VarMetaI4", VarMetaI4New, InFileID, VarIDCellI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading var I4 metadata FAIL");
      }
      if (VarMetaI4New == VarMetaI4Ref) {
         LOG_INFO("IOTest: read/write var metadata I4 test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write var metadata I4 test FAIL");
      }
      Err = IO::readMeta("VarMetaI8", VarMetaI8New, InFileID, VarIDCellI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading var I8 metadata FAIL");
      }
      if (VarMetaI8New == VarMetaI8Ref) {
         LOG_INFO("IOTest: read/write var metadata I8 test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write var metadata I8 test FAIL");
      }
      Err = IO::readMeta("VarMetaR4", VarMetaR4New, InFileID, VarIDCellI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading var R4 metadata FAIL");
      }
      if (VarMetaR4New == VarMetaR4Ref) {
         LOG_INFO("IOTest: read/write var metadata R4 test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write var metadata R4 test FAIL");
      }
      Err = IO::readMeta("VarMetaR8", VarMetaR8New, InFileID, VarIDCellI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading var R8 metadata FAIL");
      }
      if (VarMetaR8New == VarMetaR8Ref) {
         LOG_INFO("IOTest: read/write var metadata R8 test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write var metadata R8 test FAIL");
      }
      Err =
          IO::readMeta("VarMetaDescr", VarMetaDescrNew, InFileID, VarIDCellI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error reading var string metadata FAIL");
      }
      if (VarMetaDescrNew == VarMetaDescrRef) {
         LOG_INFO("IOTest: read/write var metadata string test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("IOTest: read/write var metadata string test FAIL");
      }

      // Finished reading, close file
      Err = IO::closeFile(InFileID);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error closing input file FAIL");
      }

      // Test destruction of Decompositions
      Err = IO::destroyDecomp(DecompCellI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error destroying decomp cell I4 FAIL");
      }
      Err = IO::destroyDecomp(DecompCellI8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error destroying decomp cell I8 FAIL");
      }
      Err = IO::destroyDecomp(DecompCellR4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error destroying decomp cell R4 FAIL");
      }
      Err = IO::destroyDecomp(DecompCellR8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error destroying decomp cell R8 FAIL");
      }

      Err = IO::destroyDecomp(DecompEdgeI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error destroying decomp Edge I4 FAIL");
      }
      Err = IO::destroyDecomp(DecompEdgeI8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error destroying decomp Edge I8 FAIL");
      }
      Err = IO::destroyDecomp(DecompEdgeR4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error destroying decomp Edge R4 FAIL");
      }
      Err = IO::destroyDecomp(DecompEdgeR8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error destroying decomp Edge R8 FAIL");
      }

      Err = IO::destroyDecomp(DecompVrtxI4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error destroying decomp Vrtx I4 FAIL");
      }
      Err = IO::destroyDecomp(DecompVrtxI8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error destroying decomp Vrtx I8 FAIL");
      }
      Err = IO::destroyDecomp(DecompVrtxR4);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error destroying decomp Vrtx R4 FAIL");
      }
      Err = IO::destroyDecomp(DecompVrtxR8);
      if (Err != 0) {
         RetVal += 1;
         LOG_ERROR("IOTest: error destroying decomp Vrtx R8 FAIL");
      }

      // Exit environments
      Decomp::clear();
      MachEnv::removeAll();
      if (Err == 0)
         LOG_INFO("IOTest: Successful completion");
   }
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;
} // end of main
//===-----------------------------------------------------------------------===/
