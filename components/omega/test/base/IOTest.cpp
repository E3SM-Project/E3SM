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
#include "DataTypes.h"
#include "Decomp.h"
#include "Logging.h"
#include "MachEnv.h"
#include "mpi.h"

#include <iostream>

//------------------------------------------------------------------------------
// The initialization routine for IO testing. It calls various
// init routines, including the creation of the default decomposition.

int initIOTest() {

   int Err = 0;

   // Initialize the Machine Environment class - this also creates
   // the default MachEnv. Then retrieve the default environment and
   // some needed data members.
   OMEGA::MachEnv::init(MPI_COMM_WORLD);
   OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefaultEnv();
   MPI_Comm DefComm       = DefEnv->getComm();

   // Initialize the IO system
   Err = OMEGA::IO::init(DefComm);
   if (Err != 0)
      LOG_ERROR("IOTest: error initializing parallel IO");

   // Create the default decomposition (initializes the decomposition)
   Err = OMEGA::Decomp::init();
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
      OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefaultEnv();
      MPI_Comm Comm          = DefEnv->getComm();
      OMEGA::I4 MyTask       = DefEnv->getMyTask();
      OMEGA::I4 NumTasks     = DefEnv->getNumTasks();
      bool IsMaster          = DefEnv->isMasterTask();

      // Retrieve the default decomposition
      OMEGA::Decomp *DefDecomp = OMEGA::Decomp::getDefault();
      if (DefDecomp) { // true if non-null ptr
         LOG_INFO("IOTest: Default decomp retrieval PASS");
      } else {
         LOG_INFO("IOTest: Default decomp retrieval FAIL");
         return -1;
      }

      // Create Kokkos arrays of each type and at various mesh locations
      OMEGA::I4 NCellsSize      = DefDecomp->NCellsSize;
      OMEGA::I4 NEdgesSize      = DefDecomp->NEdgesSize;
      OMEGA::I4 NVerticesSize   = DefDecomp->NVerticesSize;
      OMEGA::I4 NCellsOwned     = DefDecomp->NCellsOwned;
      OMEGA::I4 NEdgesOwned     = DefDecomp->NEdgesOwned;
      OMEGA::I4 NVerticesOwned  = DefDecomp->NVerticesOwned;
      OMEGA::I4 NCellsGlobal    = DefDecomp->NCellsGlobal;
      OMEGA::I4 NEdgesGlobal    = DefDecomp->NEdgesGlobal;
      OMEGA::I4 NVerticesGlobal = DefDecomp->NVerticesGlobal;
      OMEGA::I4 NVertLevels     = 128;

      OMEGA::HostArray2DI4 RefI4Cell("RefI4Cell", NCellsSize, NVertLevels);
      OMEGA::HostArray2DI8 RefI8Cell("RefI8Cell", NCellsSize, NVertLevels);
      OMEGA::HostArray2DR4 RefR4Cell("RefR4Cell", NCellsSize, NVertLevels);
      OMEGA::HostArray2DR8 RefR8Cell("RefR8Cell", NCellsSize, NVertLevels);

      OMEGA::HostArray2DI4 RefI4Edge("RefI4Edge", NEdgesSize, NVertLevels);
      OMEGA::HostArray2DI8 RefI8Edge("RefI8Edge", NEdgesSize, NVertLevels);
      OMEGA::HostArray2DR4 RefR4Edge("RefR4Edge", NEdgesSize, NVertLevels);
      OMEGA::HostArray2DR8 RefR8Edge("RefR8Edge", NEdgesSize, NVertLevels);

      OMEGA::HostArray2DI4 RefI4Vrtx("RefI4Vrtx", NVerticesSize, NVertLevels);
      OMEGA::HostArray2DI8 RefI8Vrtx("RefI8Vrtx", NVerticesSize, NVertLevels);
      OMEGA::HostArray2DR4 RefR4Vrtx("RefR4Vrtx", NVerticesSize, NVertLevels);
      OMEGA::HostArray2DR8 RefR8Vrtx("RefR8Vrtx", NVerticesSize, NVertLevels);

      OMEGA::HostArray1DI4 CellIDH = DefDecomp->CellIDH;
      OMEGA::HostArray1DI4 EdgeIDH = DefDecomp->EdgeIDH;
      OMEGA::HostArray1DI4 VrtxIDH = DefDecomp->VertexIDH;

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

      Err = OMEGA::IO::createDecomp(DecompCellI4, OMEGA::IO::IOTypeI4, 2,
                                    CellDims, CellArraySize, OffsetCell,
                                    OMEGA::IO::DefaultRearr);
      if (Err != 0)
         LOG_ERROR("IOTest: error creating cell decomp I4 FAIL");
      Err = OMEGA::IO::createDecomp(DecompCellI8, OMEGA::IO::IOTypeI8, 2,
                                    CellDims, CellArraySize, OffsetCell,
                                    OMEGA::IO::DefaultRearr);
      if (Err != 0)
         LOG_ERROR("IOTest: error creating cell decomp I8 FAIL");
      Err = OMEGA::IO::createDecomp(DecompCellR4, OMEGA::IO::IOTypeR4, 2,
                                    CellDims, CellArraySize, OffsetCell,
                                    OMEGA::IO::DefaultRearr);
      if (Err != 0)
         LOG_ERROR("IOTest: error creating cell decomp R4 FAIL");
      Err = OMEGA::IO::createDecomp(DecompCellR8, OMEGA::IO::IOTypeR8, 2,
                                    CellDims, CellArraySize, OffsetCell,
                                    OMEGA::IO::DefaultRearr);
      if (Err != 0)
         LOG_ERROR("IOTest: error creating cell decomp R8 FAIL");

      Err = OMEGA::IO::createDecomp(DecompEdgeI4, OMEGA::IO::IOTypeI4, 2,
                                    EdgeDims, EdgeArraySize, OffsetEdge,
                                    OMEGA::IO::DefaultRearr);
      if (Err != 0)
         LOG_ERROR("IOTest: error creating edge decomp I4 FAIL");
      Err = OMEGA::IO::createDecomp(DecompEdgeI8, OMEGA::IO::IOTypeI8, 2,
                                    EdgeDims, EdgeArraySize, OffsetEdge,
                                    OMEGA::IO::DefaultRearr);
      if (Err != 0)
         LOG_ERROR("IOTest: error creating edge decomp I8 FAIL");
      Err = OMEGA::IO::createDecomp(DecompEdgeR4, OMEGA::IO::IOTypeR4, 2,
                                    EdgeDims, EdgeArraySize, OffsetEdge,
                                    OMEGA::IO::DefaultRearr);
      if (Err != 0)
         LOG_ERROR("IOTest: error creating edge decomp R4 FAIL");
      Err = OMEGA::IO::createDecomp(DecompEdgeR8, OMEGA::IO::IOTypeR8, 2,
                                    EdgeDims, EdgeArraySize, OffsetEdge,
                                    OMEGA::IO::DefaultRearr);
      if (Err != 0)
         LOG_ERROR("IOTest: error creating edge decomp R8 FAIL");

      Err = OMEGA::IO::createDecomp(DecompVrtxI4, OMEGA::IO::IOTypeI4, 2,
                                    VrtxDims, VrtxArraySize, OffsetVrtx,
                                    OMEGA::IO::DefaultRearr);
      if (Err != 0)
         LOG_ERROR("IOTest: error creating vertex decomp I4 FAIL");
      Err = OMEGA::IO::createDecomp(DecompVrtxI8, OMEGA::IO::IOTypeI8, 2,
                                    VrtxDims, VrtxArraySize, OffsetVrtx,
                                    OMEGA::IO::DefaultRearr);
      if (Err != 0)
         LOG_ERROR("IOTest: error creating vertex decomp I8 FAIL");
      Err = OMEGA::IO::createDecomp(DecompVrtxR4, OMEGA::IO::IOTypeR4, 2,
                                    VrtxDims, VrtxArraySize, OffsetVrtx,
                                    OMEGA::IO::DefaultRearr);
      if (Err != 0)
         LOG_ERROR("IOTest: error creating vertex decomp R4 FAIL");
      Err = OMEGA::IO::createDecomp(DecompVrtxR8, OMEGA::IO::IOTypeR8, 2,
                                    VrtxDims, VrtxArraySize, OffsetVrtx,
                                    OMEGA::IO::DefaultRearr);
      if (Err != 0)
         LOG_ERROR("IOTest: error creating vertex decomp R8 FAIL");

      // Open a file for output
      int OutFileID;
      Err = OMEGA::IO::openFile(
          OutFileID, "IOTest.nc", OMEGA::IO::ModeWrite, OMEGA::IO::FmtDefault,
          OMEGA::IO::IfExists::Replace, OMEGA::IO::Precision::Double);
      if (Err != 0)
         LOG_ERROR("IOTest: error opening file for output FAIL");

      // Define array dimensions
      int DimCellID;
      int DimEdgeID;
      int DimVrtxID;
      int DimVertID;
      Err = OMEGA::IO::defineDim(OutFileID, "NVertLevels", NVertLevels,
                                 DimVertID);
      if (Err != 0)
         LOG_ERROR("IOTest: error defining vertical dimension FAIL");
      Err = OMEGA::IO::defineDim(OutFileID, "NCells", NCellsGlobal, DimCellID);
      if (Err != 0)
         LOG_ERROR("IOTest: error defining Cell dimension FAIL");
      Err = OMEGA::IO::defineDim(OutFileID, "NEdges", NEdgesGlobal, DimEdgeID);
      if (Err != 0)
         LOG_ERROR("IOTest: error defining Edge dimension FAIL");
      Err = OMEGA::IO::defineDim(OutFileID, "NVertices", NVerticesGlobal,
                                 DimVrtxID);
      if (Err != 0)
         LOG_ERROR("IOTest: error defining Vertex dimension FAIL");

      // Write some global file metadata
      OMEGA::I4 FileMetaI4Ref   = 2;
      OMEGA::I8 FileMetaI8Ref   = 4;
      OMEGA::R4 FileMetaR4Ref   = 6.789;
      OMEGA::R8 FileMetaR8Ref   = 1.23456789;
      std::string FileMetaDescr = "OMEGA IO Unit test file";

      Err = OMEGA::IO::writeMeta("FileMetaI4", FileMetaI4Ref, OutFileID,
                                 OMEGA::IO::GlobalID);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing global I4 metadata FAIL");
      Err = OMEGA::IO::writeMeta("FileMetaI8", FileMetaI8Ref, OutFileID,
                                 OMEGA::IO::GlobalID);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing global I8 metadata FAIL");
      Err = OMEGA::IO::writeMeta("FileMetaR4", FileMetaR4Ref, OutFileID,
                                 OMEGA::IO::GlobalID);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing global R4 metadata FAIL");
      Err = OMEGA::IO::writeMeta("FileMetaR8", FileMetaR8Ref, OutFileID,
                                 OMEGA::IO::GlobalID);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing global R8 metadata FAIL");
      Err = OMEGA::IO::writeMeta("FileMetaDescr", FileMetaDescr, OutFileID,
                                 OMEGA::IO::GlobalID);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing global char metadata FAIL");

      // Define variables/arrays
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

      int CellDimIDs[2] = {DimCellID, DimVertID};
      int EdgeDimIDs[2] = {DimEdgeID, DimVertID};
      int VrtxDimIDs[2] = {DimVrtxID, DimVertID};

      Err = OMEGA::IO::defineVar(OutFileID, "CellI4", OMEGA::IO::IOTypeI4, 2,
                                 CellDimIDs, VarIDCellI4);
      if (Err != 0)
         LOG_ERROR("IOTest: Error defining CellI4 array FAIL");
      Err = OMEGA::IO::defineVar(OutFileID, "CellI8", OMEGA::IO::IOTypeI8, 2,
                                 CellDimIDs, VarIDCellI8);
      if (Err != 0)
         LOG_ERROR("IOTest: Error defining CellI8 array FAIL");
      Err = OMEGA::IO::defineVar(OutFileID, "CellR4", OMEGA::IO::IOTypeR4, 2,
                                 CellDimIDs, VarIDCellR4);
      if (Err != 0)
         LOG_ERROR("IOTest: Error defining CellR4 array FAIL");
      Err = OMEGA::IO::defineVar(OutFileID, "CellR8", OMEGA::IO::IOTypeR8, 2,
                                 CellDimIDs, VarIDCellR8);
      if (Err != 0)
         LOG_ERROR("IOTest: Error defining CellR8 array FAIL");

      Err = OMEGA::IO::defineVar(OutFileID, "EdgeI4", OMEGA::IO::IOTypeI4, 2,
                                 EdgeDimIDs, VarIDEdgeI4);
      if (Err != 0)
         LOG_ERROR("IOTest: Error defining EdgeI4 array FAIL");
      Err = OMEGA::IO::defineVar(OutFileID, "EdgeI8", OMEGA::IO::IOTypeI8, 2,
                                 EdgeDimIDs, VarIDEdgeI8);
      if (Err != 0)
         LOG_ERROR("IOTest: Error defining EdgeI8 array FAIL");
      Err = OMEGA::IO::defineVar(OutFileID, "EdgeR4", OMEGA::IO::IOTypeR4, 2,
                                 EdgeDimIDs, VarIDEdgeR4);
      if (Err != 0)
         LOG_ERROR("IOTest: Error defining EdgeR4 array FAIL");
      Err = OMEGA::IO::defineVar(OutFileID, "EdgeR8", OMEGA::IO::IOTypeR8, 2,
                                 EdgeDimIDs, VarIDEdgeR8);
      if (Err != 0)
         LOG_ERROR("IOTest: Error defining EdgeR8 array FAIL");

      Err = OMEGA::IO::defineVar(OutFileID, "VrtxI4", OMEGA::IO::IOTypeI4, 2,
                                 VrtxDimIDs, VarIDVrtxI4);
      if (Err != 0)
         LOG_ERROR("IOTest: Error defining VrtxI4 array FAIL");
      Err = OMEGA::IO::defineVar(OutFileID, "VrtxI8", OMEGA::IO::IOTypeI8, 2,
                                 VrtxDimIDs, VarIDVrtxI8);
      if (Err != 0)
         LOG_ERROR("IOTest: Error defining VrtxI8 array FAIL");
      Err = OMEGA::IO::defineVar(OutFileID, "VrtxR4", OMEGA::IO::IOTypeR4, 2,
                                 VrtxDimIDs, VarIDVrtxR4);
      if (Err != 0)
         LOG_ERROR("IOTest: Error defining VrtxR4 array FAIL");
      Err = OMEGA::IO::defineVar(OutFileID, "VrtxR8", OMEGA::IO::IOTypeR8, 2,
                                 VrtxDimIDs, VarIDVrtxR8);
      if (Err != 0)
         LOG_ERROR("IOTest: Error defining VrtxR8 array FAIL");

      // Add Variable metadata just for one array
      OMEGA::I4 VarMetaI4Ref      = 3;
      OMEGA::I8 VarMetaI8Ref      = 5;
      OMEGA::R4 VarMetaR4Ref      = 5.789;
      OMEGA::R8 VarMetaR8Ref      = 2.23456789;
      std::string VarMetaDescrRef = "Test array for I4 on Cells";

      Err = OMEGA::IO::writeMeta("VarMetaI4", VarMetaI4Ref, OutFileID,
                                 VarIDCellI4);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing var I4 metadata FAIL");
      Err = OMEGA::IO::writeMeta("VarMetaI8", VarMetaI8Ref, OutFileID,
                                 VarIDCellI4);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing var I8 metadata FAIL");
      Err = OMEGA::IO::writeMeta("VarMetaR4", VarMetaR4Ref, OutFileID,
                                 VarIDCellI4);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing var R4 metadata FAIL");
      Err = OMEGA::IO::writeMeta("VarMetaR8", VarMetaR8Ref, OutFileID,
                                 VarIDCellI4);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing var R8 metadata FAIL");
      Err = OMEGA::IO::writeMeta("VarMetaDescr", VarMetaDescrRef, OutFileID,
                                 VarIDCellI4);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing var char metadata FAIL");

      // Exit define mode
      Err = OMEGA::IO::endDefinePhase(OutFileID);
      if (Err != 0)
         LOG_ERROR("IOTest: error ending define mode FAIL");

      // Write variables
      OMEGA::I4 FillI4 = -999;
      OMEGA::I8 FillI8 = -999999;
      OMEGA::R4 FillR4 = -1.234e30;
      OMEGA::R8 FillR8 = -1.23456789e30;

      Err =
          OMEGA::IO::writeArray(RefI4Cell.data(), NCellsSize * NVertLevels,
                                &FillI4, OutFileID, DecompCellI4, VarIDCellI4);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing I4 array on cells FAIL");
      Err =
          OMEGA::IO::writeArray(RefI8Cell.data(), NCellsSize * NVertLevels,
                                &FillI8, OutFileID, DecompCellI8, VarIDCellI8);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing I8 array on cells FAIL");
      Err =
          OMEGA::IO::writeArray(RefR4Cell.data(), NCellsSize * NVertLevels,
                                &FillR4, OutFileID, DecompCellR4, VarIDCellR4);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing R4 array on cells FAIL");
      Err =
          OMEGA::IO::writeArray(RefR8Cell.data(), NCellsSize * NVertLevels,
                                &FillR8, OutFileID, DecompCellR8, VarIDCellR8);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing R8 array on cells FAIL");

      Err =
          OMEGA::IO::writeArray(RefI4Edge.data(), NEdgesSize * NVertLevels,
                                &FillI4, OutFileID, DecompEdgeI4, VarIDEdgeI4);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing I4 array on Edges FAIL");
      Err =
          OMEGA::IO::writeArray(RefI8Edge.data(), NEdgesSize * NVertLevels,
                                &FillI8, OutFileID, DecompEdgeI8, VarIDEdgeI8);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing I8 array on Edges FAIL");
      Err =
          OMEGA::IO::writeArray(RefR4Edge.data(), NEdgesSize * NVertLevels,
                                &FillR4, OutFileID, DecompEdgeR4, VarIDEdgeR4);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing R4 array on Edges FAIL");
      Err =
          OMEGA::IO::writeArray(RefR8Edge.data(), NEdgesSize * NVertLevels,
                                &FillR8, OutFileID, DecompEdgeR8, VarIDEdgeR8);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing R8 array on Edges FAIL");

      Err =
          OMEGA::IO::writeArray(RefI4Vrtx.data(), NVerticesSize * NVertLevels,
                                &FillI4, OutFileID, DecompVrtxI4, VarIDVrtxI4);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing I4 array on vertices FAIL");
      Err =
          OMEGA::IO::writeArray(RefI8Vrtx.data(), NVerticesSize * NVertLevels,
                                &FillI8, OutFileID, DecompVrtxI8, VarIDVrtxI8);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing I8 array on vertices FAIL");
      Err =
          OMEGA::IO::writeArray(RefR4Vrtx.data(), NVerticesSize * NVertLevels,
                                &FillR4, OutFileID, DecompVrtxR4, VarIDVrtxR4);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing R4 array on vertices FAIL");
      Err =
          OMEGA::IO::writeArray(RefR8Vrtx.data(), NVerticesSize * NVertLevels,
                                &FillR8, OutFileID, DecompVrtxR8, VarIDVrtxR8);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing R8 array on vertices FAIL");

      // Finished writing, close file
      Err = OMEGA::IO::closeFile(OutFileID);
      if (Err != 0)
         LOG_ERROR("IOTest: error closing output file FAIL");

      // Open a file for reading to verify read/write
      int InFileID;
      Err = OMEGA::IO::openFile(InFileID, "IOTest.nc", OMEGA::IO::ModeRead);
      if (Err != 0)
         LOG_ERROR("IOTest: error opening file for reading FAIL");

      // Get dimension lengths to verify read/write of dimension info
      OMEGA::I4 NVertLevelsNew =
          OMEGA::IO::getDimLength(InFileID, "NVertLevels");
      if (NVertLevelsNew == NVertLevels) {
         LOG_INFO("IOTest: read/write vert dimension test PASS");
      } else {
         LOG_INFO("IOTest: read/write vert dimension test FAIL");
      }

      OMEGA::I4 NCellsNew = OMEGA::IO::getDimLength(InFileID, "NCells");
      if (NCellsNew == NCellsGlobal) {
         LOG_INFO("IOTest: read/write cell dimension test PASS");
      } else {
         LOG_INFO("IOTest: read/write cell dimension test FAIL");
      }

      OMEGA::I4 NEdgesNew = OMEGA::IO::getDimLength(InFileID, "NEdges");
      if (NEdgesNew == NEdgesGlobal) {
         LOG_INFO("IOTest: read/write edge dimension test PASS");
      } else {
         LOG_INFO("IOTest: read/write edge dimension test FAIL");
      }

      OMEGA::I4 NVerticesNew = OMEGA::IO::getDimLength(InFileID, "NVertices");
      if (NVerticesNew == NVerticesGlobal) {
         LOG_INFO("IOTest: read/write vertex dimension test PASS");
      } else {
         LOG_INFO("IOTest: read/write vertex dimension test FAIL");
      }

      // Read global attributes
      OMEGA::I4 FileMetaI4New;
      OMEGA::I8 FileMetaI8New;
      OMEGA::R4 FileMetaR4New;
      OMEGA::R8 FileMetaR8New;
      std::string FileMetaDescrNew;

      Err = OMEGA::IO::readMeta("FileMetaI4", FileMetaI4New, InFileID,
                                OMEGA::IO::GlobalID);
      if (Err != 0)
         LOG_ERROR("IOTest: error reading file I4 metadata FAIL");
      if (FileMetaI4New == FileMetaI4Ref) {
         LOG_INFO("IOTest: read/write file metadata I4 test PASS");
      } else {
         LOG_INFO("IOTest: read/write file metadata I4 test FAIL");
      }

      Err = OMEGA::IO::readMeta("FileMetaI8", FileMetaI8New, InFileID,
                                OMEGA::IO::GlobalID);
      if (Err != 0)
         LOG_ERROR("IOTest: error reading file I8 metadata FAIL");
      if (FileMetaI8New == FileMetaI8Ref) {
         LOG_INFO("IOTest: read/write file metadata I8 test PASS");
      } else {
         LOG_INFO("IOTest: read/write file metadata I8 test FAIL");
      }

      Err = OMEGA::IO::readMeta("FileMetaR4", FileMetaR4New, InFileID,
                                OMEGA::IO::GlobalID);
      if (Err != 0)
         LOG_ERROR("IOTest: error reading file R4 metadata FAIL");
      if (FileMetaR4New == FileMetaR4Ref) {
         LOG_INFO("IOTest: read/write file metadata R4 test PASS");
      } else {
         LOG_INFO("IOTest: read/write file metadata R4 test FAIL");
      }

      Err = OMEGA::IO::readMeta("FileMetaR8", FileMetaR8New, InFileID,
                                OMEGA::IO::GlobalID);
      if (Err != 0)
         LOG_ERROR("IOTest: error reading file R8 metadata FAIL");
      if (FileMetaR8New == FileMetaR8Ref) {
         LOG_INFO("IOTest: read/write file metadata R8 test PASS");
      } else {
         LOG_INFO("IOTest: read/write file metadata R8 test FAIL");
      }

      Err = OMEGA::IO::readMeta("FileMetaDescr", FileMetaDescrNew, InFileID,
                                OMEGA::IO::GlobalID);
      if (Err != 0)
         LOG_ERROR("IOTest: error reading file string metadata FAIL");
      if (FileMetaDescrNew == FileMetaDescr) {
         LOG_INFO("IOTest: read/write file metadata string test PASS");
      } else {
         LOG_INFO("IOTest: read/write file metadata string test FAIL");
      }

      // Read arrays
      OMEGA::HostArray2DI4 NewI4Cell("NewI4Cell", NCellsSize, NVertLevels);
      OMEGA::HostArray2DI8 NewI8Cell("NewI8Cell", NCellsSize, NVertLevels);
      OMEGA::HostArray2DR4 NewR4Cell("NewR4Cell", NCellsSize, NVertLevels);
      OMEGA::HostArray2DR8 NewR8Cell("NewR8Cell", NCellsSize, NVertLevels);

      OMEGA::HostArray2DI4 NewI4Edge("NewI4Edge", NEdgesSize, NVertLevels);
      OMEGA::HostArray2DI8 NewI8Edge("NewI8Edge", NEdgesSize, NVertLevels);
      OMEGA::HostArray2DR4 NewR4Edge("NewR4Edge", NEdgesSize, NVertLevels);
      OMEGA::HostArray2DR8 NewR8Edge("NewR8Edge", NEdgesSize, NVertLevels);

      OMEGA::HostArray2DI4 NewI4Vrtx("NewI4Vrtx", NVerticesSize, NVertLevels);
      OMEGA::HostArray2DI8 NewI8Vrtx("NewI8Vrtx", NVerticesSize, NVertLevels);
      OMEGA::HostArray2DR4 NewR4Vrtx("NewR4Vrtx", NVerticesSize, NVertLevels);
      OMEGA::HostArray2DR8 NewR8Vrtx("NewR8Vrtx", NVerticesSize, NVertLevels);

      Err = OMEGA::IO::readArray(NewI4Cell.data(), NCellsSize * NVertLevels,
                                 "CellI4", InFileID, DecompCellI4, VarIDCellI4);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing I4 array on cells FAIL");
      Err = OMEGA::IO::readArray(NewI8Cell.data(), NCellsSize * NVertLevels,
                                 "CellI8", InFileID, DecompCellI8, VarIDCellI8);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing I8 array on cells FAIL");
      Err = OMEGA::IO::readArray(NewR4Cell.data(), NCellsSize * NVertLevels,
                                 "CellR4", InFileID, DecompCellR4, VarIDCellR4);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing R4 array on cells FAIL");
      Err = OMEGA::IO::readArray(NewR8Cell.data(), NCellsSize * NVertLevels,
                                 "CellR8", InFileID, DecompCellR8, VarIDCellR8);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing R8 array on cells FAIL");

      Err = OMEGA::IO::readArray(NewI4Edge.data(), NEdgesSize * NVertLevels,
                                 "EdgeI4", InFileID, DecompEdgeI4, VarIDEdgeI4);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing I4 array on Edges FAIL");
      Err = OMEGA::IO::readArray(NewI8Edge.data(), NEdgesSize * NVertLevels,
                                 "EdgeI8", InFileID, DecompEdgeI8, VarIDEdgeI8);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing I8 array on Edges FAIL");
      Err = OMEGA::IO::readArray(NewR4Edge.data(), NEdgesSize * NVertLevels,
                                 "EdgeR4", InFileID, DecompEdgeR4, VarIDEdgeR4);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing R4 array on Edges FAIL");
      Err = OMEGA::IO::readArray(NewR8Edge.data(), NEdgesSize * NVertLevels,
                                 "EdgeR8", InFileID, DecompEdgeR8, VarIDEdgeR8);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing R8 array on Edges FAIL");

      Err = OMEGA::IO::readArray(NewI4Vrtx.data(), NVerticesSize * NVertLevels,
                                 "VrtxI4", InFileID, DecompVrtxI4, VarIDVrtxI4);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing I4 array on vertices FAIL");
      Err = OMEGA::IO::readArray(NewI8Vrtx.data(), NVerticesSize * NVertLevels,
                                 "VrtxI8", InFileID, DecompVrtxI8, VarIDVrtxI8);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing I8 array on vertices FAIL");
      Err = OMEGA::IO::readArray(NewR4Vrtx.data(), NVerticesSize * NVertLevels,
                                 "VrtxR4", InFileID, DecompVrtxR4, VarIDVrtxR4);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing R4 array on vertices FAIL");
      Err = OMEGA::IO::readArray(NewR8Vrtx.data(), NVerticesSize * NVertLevels,
                                 "VrtxR8", InFileID, DecompVrtxR8, VarIDVrtxR8);
      if (Err != 0)
         LOG_ERROR("IOTest: error writing R8 array on vertices FAIL");

      // Check that arrays match the reference cases that were written
      // Only check the owned values - these would need to be followed by
      // a halo update.

      int Err1 = 0;
      int Err2 = 0;
      int Err3 = 0;
      int Err4 = 0;
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
         LOG_INFO("IOTest: read/write array I4 on Cells test FAIL");
      }
      if (Err2 == 0) {
         LOG_INFO("IOTest: read/write array I8 on Cells test PASS");
      } else {
         LOG_INFO("IOTest: read/write array I8 on Cells test FAIL");
      }
      if (Err3 == 0) {
         LOG_INFO("IOTest: read/write array R4 on Cells test PASS");
      } else {
         LOG_INFO("IOTest: read/write array R4 on Cells test FAIL");
      }
      if (Err4 == 0) {
         LOG_INFO("IOTest: read/write array R8 on Cells test PASS");
      } else {
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
         LOG_INFO("IOTest: read/write array I4 on Edges test FAIL");
      }
      if (Err2 == 0) {
         LOG_INFO("IOTest: read/write array I8 on Edges test PASS");
      } else {
         LOG_INFO("IOTest: read/write array I8 on Edges test FAIL");
      }
      if (Err3 == 0) {
         LOG_INFO("IOTest: read/write array R4 on Edges test PASS");
      } else {
         LOG_INFO("IOTest: read/write array R4 on Edges test FAIL");
      }
      if (Err4 == 0) {
         LOG_INFO("IOTest: read/write array R8 on Edges test PASS");
      } else {
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
         LOG_INFO("IOTest: read/write array I4 on Vertices test FAIL");
      }
      if (Err2 == 0) {
         LOG_INFO("IOTest: read/write array I8 on Vertices test PASS");
      } else {
         LOG_INFO("IOTest: read/write array I8 on Vertices test FAIL");
      }
      if (Err3 == 0) {
         LOG_INFO("IOTest: read/write array R4 on Vertices test PASS");
      } else {
         LOG_INFO("IOTest: read/write array R4 on Vertices test FAIL");
      }
      if (Err4 == 0) {
         LOG_INFO("IOTest: read/write array R8 on Vertices test PASS");
      } else {
         LOG_INFO("IOTest: read/write array R8 on Vertices test FAIL");
      }

      // Read array attributes
      OMEGA::I4 VarMetaI4New;
      OMEGA::I8 VarMetaI8New;
      OMEGA::R4 VarMetaR4New;
      OMEGA::R8 VarMetaR8New;
      std::string VarMetaDescrNew;

      Err =
          OMEGA::IO::readMeta("VarMetaI4", VarMetaI4New, InFileID, VarIDCellI4);
      if (Err != 0)
         LOG_ERROR("IOTest: error reading var I4 metadata FAIL");
      if (VarMetaI4New == VarMetaI4Ref) {
         LOG_INFO("IOTest: read/write var metadata I4 test PASS");
      } else {
         LOG_INFO("IOTest: read/write var metadata I4 test FAIL");
      }
      Err =
          OMEGA::IO::readMeta("VarMetaI8", VarMetaI8New, InFileID, VarIDCellI4);
      if (Err != 0)
         LOG_ERROR("IOTest: error reading var I8 metadata FAIL");
      if (VarMetaI8New == VarMetaI8Ref) {
         LOG_INFO("IOTest: read/write var metadata I8 test PASS");
      } else {
         LOG_INFO("IOTest: read/write var metadata I8 test FAIL");
      }
      Err =
          OMEGA::IO::readMeta("VarMetaR4", VarMetaR4New, InFileID, VarIDCellI4);
      if (Err != 0)
         LOG_ERROR("IOTest: error reading var R4 metadata FAIL");
      if (VarMetaR4New == VarMetaR4Ref) {
         LOG_INFO("IOTest: read/write var metadata R4 test PASS");
      } else {
         LOG_INFO("IOTest: read/write var metadata R4 test FAIL");
      }
      Err =
          OMEGA::IO::readMeta("VarMetaR8", VarMetaR8New, InFileID, VarIDCellI4);
      if (Err != 0)
         LOG_ERROR("IOTest: error reading var R8 metadata FAIL");
      if (VarMetaR8New == VarMetaR8Ref) {
         LOG_INFO("IOTest: read/write var metadata R8 test PASS");
      } else {
         LOG_INFO("IOTest: read/write var metadata R8 test FAIL");
      }
      Err = OMEGA::IO::readMeta("VarMetaDescr", VarMetaDescrNew, InFileID,
                                VarIDCellI4);
      if (Err != 0)
         LOG_ERROR("IOTest: error reading var string metadata FAIL");
      if (VarMetaDescrNew == VarMetaDescrRef) {
         LOG_INFO("IOTest: read/write var metadata string test PASS");
      } else {
         LOG_INFO("IOTest: read/write var metadata string test FAIL");
      }

      // Finished reading, close file
      Err = OMEGA::IO::closeFile(InFileID);
      if (Err != 0)
         LOG_ERROR("IOTest: error closing input file FAIL");

      // Test destruction of Decompositions
      Err = OMEGA::IO::destroyDecomp(DecompCellI4);
      if (Err != 0)
         LOG_ERROR("IOTest: error destroying decomp cell I4 FAIL");
      Err = OMEGA::IO::destroyDecomp(DecompCellI8);
      if (Err != 0)
         LOG_ERROR("IOTest: error destroying decomp cell I8 FAIL");
      Err = OMEGA::IO::destroyDecomp(DecompCellR4);
      if (Err != 0)
         LOG_ERROR("IOTest: error destroying decomp cell R4 FAIL");
      Err = OMEGA::IO::destroyDecomp(DecompCellR8);
      if (Err != 0)
         LOG_ERROR("IOTest: error destroying decomp cell R8 FAIL");

      Err = OMEGA::IO::destroyDecomp(DecompEdgeI4);
      if (Err != 0)
         LOG_ERROR("IOTest: error destroying decomp Edge I4 FAIL");
      Err = OMEGA::IO::destroyDecomp(DecompEdgeI8);
      if (Err != 0)
         LOG_ERROR("IOTest: error destroying decomp Edge I8 FAIL");
      Err = OMEGA::IO::destroyDecomp(DecompEdgeR4);
      if (Err != 0)
         LOG_ERROR("IOTest: error destroying decomp Edge R4 FAIL");
      Err = OMEGA::IO::destroyDecomp(DecompEdgeR8);
      if (Err != 0)
         LOG_ERROR("IOTest: error destroying decomp Edge R8 FAIL");

      Err = OMEGA::IO::destroyDecomp(DecompVrtxI4);
      if (Err != 0)
         LOG_ERROR("IOTest: error destroying decomp Vrtx I4 FAIL");
      Err = OMEGA::IO::destroyDecomp(DecompVrtxI8);
      if (Err != 0)
         LOG_ERROR("IOTest: error destroying decomp Vrtx I8 FAIL");
      Err = OMEGA::IO::destroyDecomp(DecompVrtxR4);
      if (Err != 0)
         LOG_ERROR("IOTest: error destroying decomp Vrtx R4 FAIL");
      Err = OMEGA::IO::destroyDecomp(DecompVrtxR8);
      if (Err != 0)
         LOG_ERROR("IOTest: error destroying decomp Vrtx R8 FAIL");

      // Exit environments
      OMEGA::Decomp::clear();
      OMEGA::MachEnv::removeAll();
      if (Err == 0)
         LOG_INFO("IOTest: Successful completion");
   }
   Kokkos::finalize();
   MPI_Finalize();

} // end of main
//===-----------------------------------------------------------------------===/
