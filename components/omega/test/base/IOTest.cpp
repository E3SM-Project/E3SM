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
#include "Pacer.h"
#include "mpi.h"

#include <iostream>

using namespace OMEGA;

//------------------------------------------------------------------------------
// The initialization routine for IO testing. It calls various
// init routines, including the creation of the default decomposition.

void initIOTest() {

   // Initialize the Machine Environment class - this also creates
   // the default MachEnv. Then retrieve the default environment and
   // some needed data members.
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // Initialize the Logging system
   initLogging(DefEnv);
   LOG_INFO("----- Base IO Unit Testing -----");

   // Open config file
   Config("Omega");
   Config::readAll("omega.yml");

   // Initialize the IO system
   IO::init(DefComm);

   // Create the default decomposition (initializes the decomposition)
   Decomp::init();

   return;
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
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");
   {
      // Call initialization routine to create the default decomposition
      // and initialize the parallel IO library
      Error Err;
      initIOTest();

      // Retrieve the default decomposition
      Decomp *DefDecomp = Decomp::getDefault();

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
      I4 NVertLayers     = 128;

      HostArray1DI4 RefI4Vert("RefI4Vert", NVertLayers);
      HostArray1DI8 RefI8Vert("RefI8Vert", NVertLayers);
      HostArray1DR4 RefR4Vert("RefR4Vert", NVertLayers);
      HostArray1DR8 RefR8Vert("RefR8Vert", NVertLayers);
      HostArray1DR8 RefR8Time("RefR8Time", NVertLayers);

      HostArray2DI4 RefI4Cell("RefI4Cell", NCellsSize, NVertLayers);
      HostArray2DI8 RefI8Cell("RefI8Cell", NCellsSize, NVertLayers);
      HostArray2DR4 RefR4Cell("RefR4Cell", NCellsSize, NVertLayers);
      HostArray2DR8 RefR8Cell("RefR8Cell", NCellsSize, NVertLayers);
      HostArray2DR8 RefR8Tim2("RefR8Tim2", NCellsSize, NVertLayers);

      HostArray2DI4 RefI4Edge("RefI4Edge", NEdgesSize, NVertLayers);
      HostArray2DI8 RefI8Edge("RefI8Edge", NEdgesSize, NVertLayers);
      HostArray2DR4 RefR4Edge("RefR4Edge", NEdgesSize, NVertLayers);
      HostArray2DR8 RefR8Edge("RefR8Edge", NEdgesSize, NVertLayers);

      HostArray2DI4 RefI4Vrtx("RefI4Vrtx", NVerticesSize, NVertLayers);
      HostArray2DI8 RefI8Vrtx("RefI8Vrtx", NVerticesSize, NVertLayers);
      HostArray2DR4 RefR4Vrtx("RefR4Vrtx", NVerticesSize, NVertLayers);
      HostArray2DR8 RefR8Vrtx("RefR8Vrtx", NVerticesSize, NVertLayers);

      HostArray1DI4 CellIDH = DefDecomp->CellIDH;
      HostArray1DI4 EdgeIDH = DefDecomp->EdgeIDH;
      HostArray1DI4 VrtxIDH = DefDecomp->VertexIDH;

      // Create local non-distributed arrays and scalars
      I4 RefI4Scalar = -1;
      I8 RefI8Scalar = -2;
      R4 RefR4Scalar = -3.1;
      R8 RefR8Scalar = -4.56789;
      for (int K = 0; K < NVertLayers; ++K) {
         RefI4Vert(K) = K;
         RefI8Vert(K) = K * 2;
         RefR4Vert(K) = K * 3.1;
         RefR8Vert(K) = K * 4.1234567;
         RefR8Time(K) = K * 5.1234567;
      }

      // Offset arrays - initialize to -1, corresponding to entries
      // that should not be written;
      std::vector<int> OffsetCell(NCellsSize * NVertLayers, -1);
      std::vector<int> OffsetEdge(NEdgesSize * NVertLayers, -1);
      std::vector<int> OffsetVrtx(NVerticesSize * NVertLayers, -1);
      for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
         int GlobalCellAdd = CellIDH(Cell) - 1; // 0-based offset
         for (int k = 0; k < NVertLayers; ++k) {
            RefI4Cell(Cell, k)    = GlobalCellAdd * k;
            RefI8Cell(Cell, k)    = GlobalCellAdd * k * 1000000000;
            RefR4Cell(Cell, k)    = GlobalCellAdd * k * 123.45;
            RefR8Cell(Cell, k)    = GlobalCellAdd * k * 1.23456789;
            RefR8Tim2(Cell, k)    = GlobalCellAdd * k * 2.23456789;
            int VectorAdd         = Cell * NVertLayers + k;
            OffsetCell[VectorAdd] = GlobalCellAdd * NVertLayers + k;
         }
      }

      for (int Edge = 0; Edge < NEdgesOwned; ++Edge) {
         int GlobalEdgeAdd = EdgeIDH(Edge) - 1;
         for (int k = 0; k < NVertLayers; ++k) {
            RefI4Edge(Edge, k)    = GlobalEdgeAdd * k;
            RefI8Edge(Edge, k)    = GlobalEdgeAdd * k * 1000000000;
            RefR4Edge(Edge, k)    = GlobalEdgeAdd * k * 123.45;
            RefR8Edge(Edge, k)    = GlobalEdgeAdd * k * 1.23456789;
            int VectorAdd         = Edge * NVertLayers + k;
            OffsetEdge[VectorAdd] = GlobalEdgeAdd * NVertLayers + k;
         }
      }

      for (int Vrtx = 0; Vrtx < NVerticesOwned; ++Vrtx) {
         int GlobalVrtxAdd = VrtxIDH(Vrtx) - 1;
         for (int k = 0; k < NVertLayers; ++k) {
            RefI4Vrtx(Vrtx, k)    = GlobalVrtxAdd * k;
            RefI8Vrtx(Vrtx, k)    = GlobalVrtxAdd * k * 1000000000;
            RefR4Vrtx(Vrtx, k)    = GlobalVrtxAdd * k * 123.45;
            RefR8Vrtx(Vrtx, k)    = GlobalVrtxAdd * k * 1.23456789;
            int VectorAdd         = Vrtx * NVertLayers + k;
            OffsetVrtx[VectorAdd] = GlobalVrtxAdd * NVertLayers + k;
         }
      }

      // Create the needed decomposition offsets

      std::vector<int> CellDims{NCellsGlobal, NVertLayers};
      std::vector<int> EdgeDims{NEdgesGlobal, NVertLayers};
      std::vector<int> VrtxDims{NVerticesGlobal, NVertLayers};
      int CellArraySize = NCellsSize * NVertLayers;
      int EdgeArraySize = NEdgesSize * NVertLayers;
      int VrtxArraySize = NVerticesSize * NVertLayers;

      int DecompCellI4 =
          IO::createDecomp(IO::IOTypeI4, 2, CellDims, CellArraySize, OffsetCell,
                           IO::DefaultRearr);
      int DecompCellI8 =
          IO::createDecomp(IO::IOTypeI8, 2, CellDims, CellArraySize, OffsetCell,
                           IO::DefaultRearr);
      int DecompCellR4 =
          IO::createDecomp(IO::IOTypeR4, 2, CellDims, CellArraySize, OffsetCell,
                           IO::DefaultRearr);
      int DecompCellR8 =
          IO::createDecomp(IO::IOTypeR8, 2, CellDims, CellArraySize, OffsetCell,
                           IO::DefaultRearr);
      int DecompEdgeI4 =
          IO::createDecomp(IO::IOTypeI4, 2, EdgeDims, EdgeArraySize, OffsetEdge,
                           IO::DefaultRearr);
      int DecompEdgeI8 =
          IO::createDecomp(IO::IOTypeI8, 2, EdgeDims, EdgeArraySize, OffsetEdge,
                           IO::DefaultRearr);
      int DecompEdgeR4 =
          IO::createDecomp(IO::IOTypeR4, 2, EdgeDims, EdgeArraySize, OffsetEdge,
                           IO::DefaultRearr);
      int DecompEdgeR8 =
          IO::createDecomp(IO::IOTypeR8, 2, EdgeDims, EdgeArraySize, OffsetEdge,
                           IO::DefaultRearr);
      int DecompVrtxI4 =
          IO::createDecomp(IO::IOTypeI4, 2, VrtxDims, VrtxArraySize, OffsetVrtx,
                           IO::DefaultRearr);
      int DecompVrtxI8 =
          IO::createDecomp(IO::IOTypeI8, 2, VrtxDims, VrtxArraySize, OffsetVrtx,
                           IO::DefaultRearr);
      int DecompVrtxR4 =
          IO::createDecomp(IO::IOTypeR4, 2, VrtxDims, VrtxArraySize, OffsetVrtx,
                           IO::DefaultRearr);
      int DecompVrtxR8 =
          IO::createDecomp(IO::IOTypeR8, 2, VrtxDims, VrtxArraySize, OffsetVrtx,
                           IO::DefaultRearr);

      // Open a file for output
      int OutFileID;
      bool NewFile;
      IO::openFileWrite(OutFileID, "IOTest.nc", NewFile, IO::IfExists::Replace,
                        IO::FmtDefault);

      // Define array dimensions
      int DimCellID;
      int DimEdgeID;
      int DimVrtxID;
      int DimVertID;
      int DimTimeID; // unlimited time dim
      DimVertID = IO::defineDim(OutFileID, "NVertLayers", NVertLayers);
      DimCellID = IO::defineDim(OutFileID, "NCells", NCellsGlobal);
      DimEdgeID = IO::defineDim(OutFileID, "NEdges", NEdgesGlobal);
      DimVrtxID = IO::defineDim(OutFileID, "NVertices", NVerticesGlobal);
      DimTimeID = IO::defineDim(OutFileID, "Time", IO::Unlimited);

      // Write some global file metadata
      I4 FileMetaI4Ref          = 2;
      I8 FileMetaI8Ref          = 4;
      R4 FileMetaR4Ref          = 6.789;
      R8 FileMetaR8Ref          = 1.23456789;
      std::string FileMetaDescr = "OMEGA IO Unit test file";

      IO::writeMeta("FileMetaI4", FileMetaI4Ref, OutFileID, IO::GlobalID);
      IO::writeMeta("FileMetaI8", FileMetaI8Ref, OutFileID, IO::GlobalID);
      IO::writeMeta("FileMetaR4", FileMetaR4Ref, OutFileID, IO::GlobalID);
      IO::writeMeta("FileMetaR8", FileMetaR8Ref, OutFileID, IO::GlobalID);
      IO::writeMeta("FileMetaDescr", FileMetaDescr, OutFileID, IO::GlobalID);
      IO::writeMeta("StringLiteral", "MyString", OutFileID, IO::GlobalID);

      // Define variables/arrays
      // Test every data type up to 2-d and use two of the arrays to
      // test the unlimited time dimension

      int VertDimIDs[1] = {DimVertID};
      int CellDimIDs[2] = {DimCellID, DimVertID};
      int EdgeDimIDs[2] = {DimEdgeID, DimVertID};
      int VrtxDimIDs[2] = {DimVrtxID, DimVertID};
      int TimeDimIDs[2] = {DimTimeID, DimVertID};
      int Tim2DimIDs[3] = {DimTimeID, DimCellID, DimVertID};

      int VarIDScalarI4 =
          IO::defineVar(OutFileID, "ScalarI4", IO::IOTypeI4, 0, nullptr);
      int VarIDScalarI8 =
          IO::defineVar(OutFileID, "ScalarI8", IO::IOTypeI8, 0, nullptr);
      int VarIDScalarR4 =
          IO::defineVar(OutFileID, "ScalarR4", IO::IOTypeR4, 0, nullptr);
      int VarIDScalarR8 =
          IO::defineVar(OutFileID, "ScalarR8", IO::IOTypeR8, 0, nullptr);
      int VarIDI4Vert =
          IO::defineVar(OutFileID, "I4Vert", IO::IOTypeI4, 1, VertDimIDs);
      int VarIDI8Vert =
          IO::defineVar(OutFileID, "I8Vert", IO::IOTypeI8, 1, VertDimIDs);
      int VarIDR4Vert =
          IO::defineVar(OutFileID, "R4Vert", IO::IOTypeR4, 1, VertDimIDs);
      int VarIDR8Time =
          IO::defineVar(OutFileID, "R8Time", IO::IOTypeR8, 2, TimeDimIDs);
      int VarIDCellI4 =
          IO::defineVar(OutFileID, "CellI4", IO::IOTypeI4, 2, CellDimIDs);
      int VarIDCellI8 =
          IO::defineVar(OutFileID, "CellI8", IO::IOTypeI8, 2, CellDimIDs);
      int VarIDCellR4 =
          IO::defineVar(OutFileID, "CellR4", IO::IOTypeR4, 2, CellDimIDs);
      int VarIDTimeR8 =
          IO::defineVar(OutFileID, "TimeR8", IO::IOTypeR8, 3, Tim2DimIDs);
      int VarIDEdgeI4 =
          IO::defineVar(OutFileID, "EdgeI4", IO::IOTypeI4, 2, EdgeDimIDs);
      int VarIDEdgeI8 =
          IO::defineVar(OutFileID, "EdgeI8", IO::IOTypeI8, 2, EdgeDimIDs);
      int VarIDEdgeR4 =
          IO::defineVar(OutFileID, "EdgeR4", IO::IOTypeR4, 2, EdgeDimIDs);
      int VarIDEdgeR8 =
          IO::defineVar(OutFileID, "EdgeR8", IO::IOTypeR8, 2, EdgeDimIDs);
      int VarIDVrtxI4 =
          IO::defineVar(OutFileID, "VrtxI4", IO::IOTypeI4, 2, VrtxDimIDs);
      int VarIDVrtxI8 =
          IO::defineVar(OutFileID, "VrtxI8", IO::IOTypeI8, 2, VrtxDimIDs);
      int VarIDVrtxR4 =
          IO::defineVar(OutFileID, "VrtxR4", IO::IOTypeR4, 2, VrtxDimIDs);
      int VarIDVrtxR8 =
          IO::defineVar(OutFileID, "VrtxR8", IO::IOTypeR8, 2, VrtxDimIDs);

      // Add Variable metadata just for one array
      I4 VarMetaI4Ref             = 3;
      I8 VarMetaI8Ref             = 5;
      R4 VarMetaR4Ref             = 5.789;
      R8 VarMetaR8Ref             = 2.23456789;
      std::string VarMetaDescrRef = "Test array for I4 on Cells";

      IO::writeMeta("VarMetaI4", VarMetaI4Ref, OutFileID, VarIDCellI4);
      IO::writeMeta("VarMetaI8", VarMetaI8Ref, OutFileID, VarIDCellI4);
      IO::writeMeta("VarMetaR4", VarMetaR4Ref, OutFileID, VarIDCellI4);
      IO::writeMeta("VarMetaR8", VarMetaR8Ref, OutFileID, VarIDCellI4);
      IO::writeMeta("VarMetaDescr", VarMetaDescrRef, OutFileID, VarIDCellI4);

      // Exit define mode
      IO::endDefinePhase(OutFileID);

      // Write variables
      I4 FillI4 = -999;
      I8 FillI8 = -999999;
      R4 FillR4 = -1.234e30;
      R8 FillR8 = -1.23456789e30;

      // Write non-distributed variables
      IO::writeNDVar(&RefI4Scalar, OutFileID, VarIDScalarI4);
      IO::writeNDVar(&RefI8Scalar, OutFileID, VarIDScalarI8);
      IO::writeNDVar(&RefR4Scalar, OutFileID, VarIDScalarR4);
      IO::writeNDVar(&RefR8Scalar, OutFileID, VarIDScalarR8);
      IO::writeNDVar(RefI4Vert.data(), OutFileID, VarIDI4Vert);
      IO::writeNDVar(RefI8Vert.data(), OutFileID, VarIDI8Vert);
      IO::writeNDVar(RefR4Vert.data(), OutFileID, VarIDR4Vert);
      // Write R8 arrays as two time slices - the first frame here with
      // the second frame written after re-open
      std::vector<int> DimLengths(1);
      DimLengths[0] = NVertLayers;
      IO::writeNDVar(RefR8Vert.data(), OutFileID, VarIDR8Time, 0, &DimLengths);

      // Write distributed arrays
      IO::writeArray(RefI4Cell.data(), NCellsSize * NVertLayers, &FillI4,
                     OutFileID, DecompCellI4, VarIDCellI4);
      IO::writeArray(RefI8Cell.data(), NCellsSize * NVertLayers, &FillI8,
                     OutFileID, DecompCellI8, VarIDCellI8);
      IO::writeArray(RefR4Cell.data(), NCellsSize * NVertLayers, &FillR4,
                     OutFileID, DecompCellR4, VarIDCellR4);
      // Write R8 cell data as two time slices - this is first frame
      // Second frame written below
      IO::writeArray(RefR8Cell.data(), NCellsSize * NVertLayers, &FillR8,
                     OutFileID, DecompCellR8, VarIDTimeR8, 0);
      IO::writeArray(RefI4Edge.data(), NEdgesSize * NVertLayers, &FillI4,
                     OutFileID, DecompEdgeI4, VarIDEdgeI4);
      IO::writeArray(RefI8Edge.data(), NEdgesSize * NVertLayers, &FillI8,
                     OutFileID, DecompEdgeI8, VarIDEdgeI8);
      IO::writeArray(RefR4Edge.data(), NEdgesSize * NVertLayers, &FillR4,
                     OutFileID, DecompEdgeR4, VarIDEdgeR4);
      IO::writeArray(RefR8Edge.data(), NEdgesSize * NVertLayers, &FillR8,
                     OutFileID, DecompEdgeR8, VarIDEdgeR8);
      IO::writeArray(RefI4Vrtx.data(), NVerticesSize * NVertLayers, &FillI4,
                     OutFileID, DecompVrtxI4, VarIDVrtxI4);
      IO::writeArray(RefI8Vrtx.data(), NVerticesSize * NVertLayers, &FillI8,
                     OutFileID, DecompVrtxI8, VarIDVrtxI8);
      IO::writeArray(RefR4Vrtx.data(), NVerticesSize * NVertLayers, &FillR4,
                     OutFileID, DecompVrtxR4, VarIDVrtxR4);
      IO::writeArray(RefR8Vrtx.data(), NVerticesSize * NVertLayers, &FillR8,
                     OutFileID, DecompVrtxR8, VarIDVrtxR8);

      // Finished writing, close file
      IO::closeFile(OutFileID);

      // Re-open to write additional frames for multi-frame fields
      // Open a file for output
      IO::openFileWrite(OutFileID, "IOTest.nc", NewFile, IO::IfExists::Append,
                        IO::FmtDefault);

      // write second frame data
      IO::writeNDVar(RefR8Time.data(), OutFileID, VarIDR8Time, 1, &DimLengths);
      IO::writeArray(RefR8Tim2.data(), NCellsSize * NVertLayers, &FillR8,
                     OutFileID, DecompCellR8, VarIDTimeR8, 1);

      // Finished writing again, close file
      IO::closeFile(OutFileID);

      // Open a file for reading to verify read/write
      int InFileID;
      IO::openFileRead(InFileID, "IOTest.nc");

      // Get dimension lengths to verify read/write of dimension info
      I4 NVertLayersID;
      I4 NVertLayersNew;
      Err = IO::getDimFromFile(InFileID, "NVertLayers", NVertLayersID,
                               NVertLayersNew);
      CHECK_ERROR_ABORT(Err, "IOTest: failed to get NVertLayers from file");

      I4 NCellsNewID;
      I4 NCellsNew;
      Err = IO::getDimFromFile(InFileID, "NCells", NCellsNewID, NCellsNew);
      CHECK_ERROR_ABORT(Err, "IOTest: failed to get NCells from file");

      I4 NEdgesNewID;
      I4 NEdgesNew;
      Err = IO::getDimFromFile(InFileID, "NEdges", NEdgesNewID, NEdgesNew);
      CHECK_ERROR_ABORT(Err, "IOTest: failed to get NEdges from file");

      I4 NVerticesNewID;
      I4 NVerticesNew;
      Err = IO::getDimFromFile(InFileID, "NVertices", NVerticesNewID,
                               NVerticesNew);
      CHECK_ERROR_ABORT(Err, "IOTest: failed to get NVertices from file");

      I4 NTimeNewID;
      I4 NTimeNew;
      Err = IO::getDimFromFile(InFileID, "Time", NTimeNewID, NTimeNew);
      CHECK_ERROR_ABORT(Err, "IOTest: failed to get NTime dimension from file");

      // Read global attributes
      I4 FileMetaI4New;
      I8 FileMetaI8New;
      R4 FileMetaR4New;
      R8 FileMetaR8New;
      std::string FileMetaDescrNew;

      Err = IO::readMeta("FileMetaI4", FileMetaI4New, InFileID, IO::GlobalID);
      CHECK_ERROR_ABORT(Err, "IOTest FAIL: error reading file I4 metadata");
      if (FileMetaI4New != FileMetaI4Ref)
         ABORT_ERROR("IOTest: read I4 file metadata FAIL with bad data");

      Err = IO::readMeta("FileMetaI8", FileMetaI8New, InFileID, IO::GlobalID);
      CHECK_ERROR_ABORT(Err, "IOTest FAIL: error reading file I8 metadata");
      if (FileMetaI8New != FileMetaI8Ref)
         ABORT_ERROR("IOTest: read I8 file metadata FAIL with bad data");

      Err = IO::readMeta("FileMetaR4", FileMetaR4New, InFileID, IO::GlobalID);
      CHECK_ERROR_ABORT(Err, "IOTest FAIL: error reading file R4 metadata");
      if (FileMetaR4New != FileMetaR4Ref)
         ABORT_ERROR("IOTest: read R4 file metadata FAIL with bad data");

      Err = IO::readMeta("FileMetaR8", FileMetaR8New, InFileID, IO::GlobalID);
      CHECK_ERROR_ABORT(Err, "IOTest FAIL: error reading file R8 metadata");
      if (FileMetaR8New != FileMetaR8Ref)
         ABORT_ERROR("IOTest: read R8 file metadata FAIL with bad data");

      Err = IO::readMeta("FileMetaDescr", FileMetaDescrNew, InFileID,
                         IO::GlobalID);
      CHECK_ERROR_ABORT(Err, "IOTest FAIL: error reading file string metadata");
      if (FileMetaDescrNew != FileMetaDescr)
         ABORT_ERROR("IOTest: read string file metadata FAIL with bad data");

      std::string MyStringNew;
      Err = IO::readMeta("StringLiteral", MyStringNew, InFileID, IO::GlobalID);
      CHECK_ERROR_ABORT(Err, "IOTest FAIL: reading string literal metadata");
      if (MyStringNew != "MyString")
         ABORT_ERROR("IOTest: read string literal metadata FAIL with bad data");

      // Read new variables

      I4 NewI4Scalar = 0;
      I8 NewI8Scalar = 0;
      R4 NewR4Scalar = 0;
      R8 NewR8Scalar = 0;

      HostArray1DI4 NewI4Vert("NewI4Vert", NVertLayers);
      HostArray1DI8 NewI8Vert("NewI8Vert", NVertLayers);
      HostArray1DR4 NewR4Vert("NewR4Vert", NVertLayers);
      HostArray1DR8 NewR8Vert("NewR8Vert", NVertLayers);
      HostArray1DR8 NewR8Time("NewR8Time", NVertLayers);

      HostArray2DI4 NewI4Cell("NewI4Cell", NCellsSize, NVertLayers);
      HostArray2DI8 NewI8Cell("NewI8Cell", NCellsSize, NVertLayers);
      HostArray2DR4 NewR4Cell("NewR4Cell", NCellsSize, NVertLayers);
      HostArray2DR8 NewR8Cell("NewR8Cell", NCellsSize, NVertLayers);
      HostArray2DR8 NewR8Tim2("NewR8Tim2", NCellsSize, NVertLayers);

      HostArray2DI4 NewI4Edge("NewI4Edge", NEdgesSize, NVertLayers);
      HostArray2DI8 NewI8Edge("NewI8Edge", NEdgesSize, NVertLayers);
      HostArray2DR4 NewR4Edge("NewR4Edge", NEdgesSize, NVertLayers);
      HostArray2DR8 NewR8Edge("NewR8Edge", NEdgesSize, NVertLayers);

      HostArray2DI4 NewI4Vrtx("NewI4Vrtx", NVerticesSize, NVertLayers);
      HostArray2DI8 NewI8Vrtx("NewI8Vrtx", NVerticesSize, NVertLayers);
      HostArray2DR4 NewR4Vrtx("NewR4Vrtx", NVerticesSize, NVertLayers);
      HostArray2DR8 NewR8Vrtx("NewR8Vrtx", NVerticesSize, NVertLayers);

      // Read non-distributed variables
      Err = IO::readNDVar(&NewI4Scalar, "ScalarI4", InFileID, VarIDScalarI4);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read scalar I4 NDVar")

      Err = IO::readNDVar(&NewI8Scalar, "ScalarI8", InFileID, VarIDScalarI8);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read scalar I8 NDVar")

      Err = IO::readNDVar(&NewR4Scalar, "ScalarR4", InFileID, VarIDScalarR4);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read scalar R4 NDVar")

      Err = IO::readNDVar(&NewR8Scalar, "ScalarR8", InFileID, VarIDScalarR8);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read scalar R8 NDVar")

      Err = IO::readNDVar(NewI4Vert.data(), "I4Vert", InFileID, VarIDI4Vert);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read I4 vertical NDVar")

      Err = IO::readNDVar(NewI8Vert.data(), "I8Vert", InFileID, VarIDI8Vert);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read I8 vertical NDVar")

      Err = IO::readNDVar(NewR4Vert.data(), "R4Vert", InFileID, VarIDR4Vert);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read R4 vertical NDVar")

      // Read R8 data as two time slices
      Err = IO::readNDVar(NewR8Vert.data(), "R8Time", InFileID, VarIDR8Time, 0,
                          &DimLengths);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read R8 NDVar time slice 0")
      Err = IO::readNDVar(NewR8Time.data(), "R8Time", InFileID, VarIDR8Time, 1,
                          &DimLengths);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read R8 NDVar time slice 1")

      // Read distributed arrays
      Err = IO::readArray(NewI4Cell.data(), NCellsSize * NVertLayers, "CellI4",
                          InFileID, DecompCellI4, VarIDCellI4);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read I4 cell array")
      Err = IO::readArray(NewI8Cell.data(), NCellsSize * NVertLayers, "CellI8",
                          InFileID, DecompCellI8, VarIDCellI8);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read I8 cell array")
      Err = IO::readArray(NewR4Cell.data(), NCellsSize * NVertLayers, "CellR4",
                          InFileID, DecompCellR4, VarIDCellR4);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read R4 cell array")
      // Read R8 Cell data as two time slices
      Err = IO::readArray(NewR8Cell.data(), NCellsSize * NVertLayers, "TimeR8",
                          InFileID, DecompCellR8, VarIDTimeR8, 0);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read R8 cell array, time 0")
      Err = IO::readArray(NewR8Tim2.data(), NCellsSize * NVertLayers, "TimeR8",
                          InFileID, DecompCellR8, VarIDTimeR8, 1);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read R8 cell array, time 1")

      Err = IO::readArray(NewI4Edge.data(), NEdgesSize * NVertLayers, "EdgeI4",
                          InFileID, DecompEdgeI4, VarIDEdgeI4);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read I4 edge array")
      Err = IO::readArray(NewI8Edge.data(), NEdgesSize * NVertLayers, "EdgeI8",
                          InFileID, DecompEdgeI8, VarIDEdgeI8);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read I8 edge array")
      Err = IO::readArray(NewR4Edge.data(), NEdgesSize * NVertLayers, "EdgeR4",
                          InFileID, DecompEdgeR4, VarIDEdgeR4);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read R4 edge array")
      Err = IO::readArray(NewR8Edge.data(), NEdgesSize * NVertLayers, "EdgeR8",
                          InFileID, DecompEdgeR8, VarIDEdgeR8);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read R8 edge array")

      Err = IO::readArray(NewI4Vrtx.data(), NVerticesSize * NVertLayers,
                          "VrtxI4", InFileID, DecompVrtxI4, VarIDVrtxI4);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read I4 vertex array")
      Err = IO::readArray(NewI8Vrtx.data(), NVerticesSize * NVertLayers,
                          "VrtxI8", InFileID, DecompVrtxI8, VarIDVrtxI8);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read I8 vertex array")
      Err = IO::readArray(NewR4Vrtx.data(), NVerticesSize * NVertLayers,
                          "VrtxR4", InFileID, DecompVrtxR4, VarIDVrtxR4);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read R4 vertex array")
      Err = IO::readArray(NewR8Vrtx.data(), NVerticesSize * NVertLayers,
                          "VrtxR8", InFileID, DecompVrtxR8, VarIDVrtxR8);
      CHECK_ERROR_ABORT(Err, "IOTest: Failed to read R8 vertex array")

      // Check that variables match the reference cases that were written
      // Only check the owned values for distributed arrays - these would need
      // to be followed by a halo update.

      if (NewI4Scalar != RefI4Scalar)
         ABORT_ERROR("IOTest: read/write scalar I4 FAIL");
      if (NewI8Scalar != RefI8Scalar)
         ABORT_ERROR("IOTest: read/write scalar I8 test FAIL");
      if (NewR4Scalar != RefR4Scalar)
         ABORT_ERROR("IOTest: read/write scalar R4 test FAIL");
      if (NewR8Scalar != RefR8Scalar)
         ABORT_ERROR("IOTest: read/write scalar R8 test FAIL");

      int Err1 = 0;
      int Err2 = 0;
      int Err3 = 0;
      int Err4 = 0;
      int Err5 = 0;

      for (int k = 0; k < NVertLayers; ++k) {
         if (NewI4Vert(k) != RefI4Vert(k))
            Err1++;
         if (NewI8Vert(k) != RefI8Vert(k))
            Err2++;
         if (NewR4Vert(k) != RefR4Vert(k))
            Err3++;
         if (NewR8Vert(k) != RefR8Vert(k))
            Err4++;
         if (NewR8Time(k) != RefR8Time(k))
            Err5++;
      }
      if (Err1 > 0)
         ABORT_ERROR("IOTest: read/write vert I4 vector test FAIL");
      if (Err2 > 0)
         ABORT_ERROR("IOTest: read/write vert I8 vector test FAIL");
      if (Err3 > 0)
         ABORT_ERROR("IOTest: read/write vert R4 vector test FAIL");
      if (Err4 > 0)
         ABORT_ERROR("IOTest: read/write vert R8 vector test frame 0 FAIL");
      if (Err5 > 0)
         ABORT_ERROR("IOTest: read/write vert R8 vector test frame 1 FAIL");

      Err1 = 0;
      Err2 = 0;
      Err3 = 0;
      Err4 = 0;
      Err5 = 0;
      for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
         for (int k = 0; k < NVertLayers; ++k) {
            if (NewI4Cell(Cell, k) != RefI4Cell(Cell, k))
               Err1++;
            if (NewI8Cell(Cell, k) != RefI8Cell(Cell, k))
               Err2++;
            if (NewR4Cell(Cell, k) != RefR4Cell(Cell, k))
               Err3++;
            if (NewR8Cell(Cell, k) != RefR8Cell(Cell, k))
               Err4++;
            if (NewR8Tim2(Cell, k) != RefR8Tim2(Cell, k))
               Err5++;
         }
      }
      if (Err1 > 0)
         ABORT_ERROR("IOTest: read/write array I4 on Cells test FAIL");
      if (Err2 > 0)
         ABORT_ERROR("IOTest: read/write array I8 on Cells test FAIL");
      if (Err3 > 0)
         ABORT_ERROR("IOTest: read/write array R4 on Cells test FAIL");
      if (Err4 > 0)
         ABORT_ERROR("IOTest: read/write array R8 on Cells test frame 0 FAIL");
      if (Err5 > 0)
         ABORT_ERROR("IOTest: read/write array R8 on Cells test frame 1 FAIL");

      Err1 = 0;
      Err2 = 0;
      Err3 = 0;
      Err4 = 0;
      for (int Edge = 0; Edge < NEdgesOwned; ++Edge) {
         for (int k = 0; k < NVertLayers; ++k) {
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
      if (Err1 > 0)
         ABORT_ERROR("IOTest: read/write array I4 on Edges test FAIL");
      if (Err2 > 0)
         ABORT_ERROR("IOTest: read/write array I8 on Edges test FAIL");
      if (Err3 > 0)
         ABORT_ERROR("IOTest: read/write array R4 on Edges test FAIL");
      if (Err4 > 0)
         ABORT_ERROR("IOTest: read/write array R8 on Edges test FAIL");

      Err1 = 0;
      Err2 = 0;
      Err3 = 0;
      Err4 = 0;
      for (int Vrtx = 0; Vrtx < NVerticesOwned; ++Vrtx) {
         for (int k = 0; k < NVertLayers; ++k) {
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
      if (Err1 > 0)
         ABORT_ERROR("IOTest: read/write array I4 on Vertices test FAIL");
      if (Err2 > 0)
         ABORT_ERROR("IOTest: read/write array I8 on Vertices test FAIL");
      if (Err3 > 0)
         ABORT_ERROR("IOTest: read/write array R4 on Vertices test FAIL");
      if (Err4 > 0)
         ABORT_ERROR("IOTest: read/write array R8 on Vertices test FAIL");

      // Read array attributes
      I4 VarMetaI4New;
      I8 VarMetaI8New;
      R4 VarMetaR4New;
      R8 VarMetaR8New;
      std::string VarMetaDescrNew;

      Err = IO::readMeta("VarMetaI4", VarMetaI4New, InFileID, VarIDCellI4);
      CHECK_ERROR_ABORT(Err, "IOTest: error reading var I4 metadata FAIL");
      if (VarMetaI4New != VarMetaI4Ref)
         ABORT_ERROR("IOTest: read var metadata I4 test FAIL with bad data");

      Err = IO::readMeta("VarMetaI8", VarMetaI8New, InFileID, VarIDCellI4);
      CHECK_ERROR_ABORT(Err, "IOTest: error reading var I8 metadata FAIL");
      if (VarMetaI8New != VarMetaI8Ref)
         ABORT_ERROR("IOTest: read var metadata I8 test FAIL with bad data");

      Err = IO::readMeta("VarMetaR4", VarMetaR4New, InFileID, VarIDCellI4);
      CHECK_ERROR_ABORT(Err, "IOTest: error reading var R4 metadata FAIL");
      if (VarMetaR4New != VarMetaR4Ref)
         ABORT_ERROR("IOTest: read var metadata R4 test FAIL with bad data");

      Err = IO::readMeta("VarMetaR8", VarMetaR8New, InFileID, VarIDCellI4);
      CHECK_ERROR_ABORT(Err, "IOTest: error reading var R8 metadata FAIL");
      if (VarMetaR8New != VarMetaR8Ref)
         ABORT_ERROR("IOTest: read var metadata R8 test FAIL with bad data");

      Err =
          IO::readMeta("VarMetaDescr", VarMetaDescrNew, InFileID, VarIDCellI4);
      CHECK_ERROR_ABORT(Err, "IOTest: error reading var string metadata FAIL");
      if (VarMetaDescrNew != VarMetaDescrRef)
         ABORT_ERROR("IOTest: read var metadata string test FAIL bad data");

      // Finished reading, close file
      IO::closeFile(InFileID);

      // Test destruction of Decompositions
      IO::destroyDecomp(DecompCellI4);
      IO::destroyDecomp(DecompCellI8);
      IO::destroyDecomp(DecompCellR4);
      IO::destroyDecomp(DecompCellR8);
      IO::destroyDecomp(DecompEdgeI4);
      IO::destroyDecomp(DecompEdgeI8);
      IO::destroyDecomp(DecompEdgeR4);
      IO::destroyDecomp(DecompEdgeR8);
      IO::destroyDecomp(DecompVrtxI4);
      IO::destroyDecomp(DecompVrtxI8);
      IO::destroyDecomp(DecompVrtxR4);
      IO::destroyDecomp(DecompVrtxR8);

      // Exit environments
      Decomp::clear();
      MachEnv::removeAll();

      LOG_INFO("IOTest: Successful completion");
   }
   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
