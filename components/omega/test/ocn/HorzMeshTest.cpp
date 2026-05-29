//===-- Test driver for OMEGA Decomp -----------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA mesh class
///
/// This driver tests that the OMEGA mesh class member variables are read in
/// correctly from a sample shperical mesh file.
//
//===-----------------------------------------------------------------------===/

#include "HorzMesh.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Error.h"
#include "Field.h"
#include "Halo.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "TimeMgr.h"
#include "mpi.h"

using namespace OMEGA;

//------------------------------------------------------------------------------
// The initialization routine for Mesh testing. It calls various
// init routines, including the creation of the default decomposition.

void initHorzMeshTest() {

   // Initialize the Machine Environment class - this also creates
   // the default MachEnv. Then retrieve the default environment and
   // some needed data members.
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // Initialize the Logging system
   initLogging(DefEnv);
   LOG_INFO("------ Horz Mesh Unit Tests ------");

   // Open config file
   Config("Omega");
   Config::readAll("omega.yml");

   // Initialize the IO system
   IO::init(DefComm);

   // Create the default decomposition (initializes the decomposition)
   Decomp::init();

   // Initialize the default halo
   int Err = Halo::init();
   if (Err != 0)
      ABORT_ERROR("HorzMeshTest: error initializing default halo");

   // Create a dummy model clock
   Calendar::init("No Leap");
   TimeInstant StartTime(0, 1, 1, 0, 0, 0.0);
   TimeInterval TimeStep(1, TimeUnits::Hours);
   Clock ModelClockTmp(StartTime, TimeStep);
   Clock *ModelClock = &ModelClockTmp;

   // Initialize the default mesh
   Field::init(ModelClock);
   IOStream::init(ModelClock);
   HorzMesh::init(ModelClock);

   return;
}

//------------------------------------------------------------------------------
// Computes the distance of a x,y,z coordinate from the origin
R8 distance(R8 X, R8 Y, R8 Z) {

   R8 Dist = sqrt(X * X + Y * Y + Z * Z);
   return Dist;
}

//------------------------------------------------------------------------------
// Computes the distance between two lon/lat points on the sphere
R8 sphereDistance(R8 Lon1, R8 Lat1, R8 Lon2, R8 Lat2) {

   R8 Arg = sqrt(pow(sin(0.5 * (Lat1 - Lat2)), 2) +
                 cos(Lat2) * cos(Lat1) * pow(sin(0.5 * (Lon1 - Lon2)), 2));
   return 2.0 * asin(Arg);
}

//------------------------------------------------------------------------------
// Computes the longitude of a point given its Cartesian coordinates
R8 computeLon(R8 X, R8 Y, R8 Z) {

   R8 Lon = atan2(Y, X);
   R8 Pi  = 4.0 * atan(1.0);
   if (Lon < 0.0) {
      Lon = 2.0 * Pi + Lon;
   }
   return Lon;
}

//------------------------------------------------------------------------------
// Computes the latitude of a point given its Cartesian coordinates
R8 computeLat(R8 X, R8 Y, R8 Z) {

   R8 Dist = distance(X, Y, Z);
   R8 Lat  = asin(Z / Dist);

   return Lat;
}

//------------------------------------------------------------------------------
// Computes coriolis parameter for a given latitude
R8 coriolis(R8 Lat) {

   R8 Omega = 7.29212e-5;
   R8 F     = 2.0 * Omega * sin(Lat);

   return F;
}

//------------------------------------------------------------------------------
// The test driver for Mesh-> This tests the decomposition of a sample
// horizontal domain and verifies the mesh is read in correctly.
//
int main(int argc, char *argv[]) {

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");
   {

      int Err = 0;
      R8 Tol  = 1e-6;
      R8 Pi   = 4.0 * atan(1.0);

      // Call initialization routine to create the default decomposition
      initHorzMeshTest();

      // Get MPI vars if needed
      MachEnv *DefEnv = MachEnv::getDefault();
      MPI_Comm Comm   = DefEnv->getComm();
      I4 MyTask       = DefEnv->getMyTask();
      I4 NumTasks     = DefEnv->getNumTasks();
      bool IsMaster   = DefEnv->isMasterTask();

      // Test retrieval of the default decomposition
      Decomp *DefDecomp = Decomp::getDefault();
      if (DefDecomp == nullptr)
         ABORT_ERROR("HorzMeshTest: Default decomp retrieval FAIL");

      // Retrieve default mesh
      HorzMesh *Mesh = HorzMesh::getDefault();

      // Test sum of local mesh cells
      // Get the global sum of all local cell counts
      // Tests that the correct cell counts have been retrieved from the Decomp
      // object
      I4 SumCells;
      I4 LocCells;
      LocCells = Mesh->NCellsOwned;
      Err = MPI_Allreduce(&LocCells, &SumCells, 1, MPI_INT32_T, MPI_SUM, Comm);
      if (Err != MPI_SUCCESS)
         ABORT_ERROR("HorzMeshTest: MPI error summing cells");

      if (SumCells != Mesh->NCellsGlobal)
         ABORT_ERROR("HorzMeshTest: Sum cell ID test FAIL {} {}", SumCells,
                     Mesh->NCellsGlobal);

      // Test that cell centers are on sphere
      // Check that all cell centers are a uniform distance from the origin
      // Tests that the Cartesian coordinates for cell centers have been read in
      // corectly
      R8 SphereRadius =
          distance(Mesh->XCellH(0), Mesh->YCellH(0), Mesh->ZCellH(0));
      R8 Dist;
      I4 Count = 0;
      for (int Cell = 0; Cell < LocCells; Cell++) {
         Dist = distance(Mesh->XCellH(Cell), Mesh->YCellH(Cell),
                         Mesh->ZCellH(Cell));
         if (abs(SphereRadius - Dist) > Tol)
            Count++;
      }

      if (Count > 0)
         ABORT_ERROR("HorzMeshTest: Cell sphere radius test FAIL");

      // Test lon/lat coordinates of cell centers
      // Convert Cartesian coordinates to lon/lat and check these agree with the
      // values that have been read in
      // Tests that the lon/lat coordinates for cell
      // centers have been read in correctly
      R8 Lon;
      R8 Lat;
      Count = 0;
      for (int Cell = 0; Cell < LocCells; Cell++) {

         Lon = computeLon(Mesh->XCellH(Cell), Mesh->YCellH(Cell),
                          Mesh->ZCellH(Cell));
         Lat = computeLat(Mesh->XCellH(Cell), Mesh->YCellH(Cell),
                          Mesh->ZCellH(Cell));

         if (abs(Lon - Mesh->LonCellH(Cell)) > Tol)
            Count++;
         if (abs(Lat - Mesh->LatCellH(Cell)) > Tol)
            Count++;
      }

      if (Count > 0)
         ABORT_ERROR("HorzMeshTest: Cell lon/lat test FAIL");

      // Test sum of local mesh edges
      // Get the global sum of all local edge counts
      // Tests that the correct edge counts have been retrieved from the Decomp
      // object
      I4 SumEdges;
      I4 LocEdges;
      LocEdges = Mesh->NEdgesOwned;
      Err = MPI_Allreduce(&LocEdges, &SumEdges, 1, MPI_INT32_T, MPI_SUM, Comm);

      if (SumEdges != Mesh->NEdgesGlobal)
         ABORT_ERROR("HorzMeshTest: Sum edge ID test FAIL {} {}", SumEdges,
                     Mesh->NEdgesGlobal);

      // Test that edge coordinates are on sphere
      // Check that all edge centers are a uniform distance from the origin
      // Tests that the Cartesian coordinates for edge centers have been read in
      // correctly
      SphereRadius =
          distance(Mesh->XEdgeH(0), Mesh->YEdgeH(0), Mesh->ZEdgeH(0));
      Count = 0;
      for (int Edge = 0; Edge < LocEdges; Edge++) {
         Dist = distance(Mesh->XEdgeH(Edge), Mesh->YEdgeH(Edge),
                         Mesh->ZEdgeH(Edge));
         if (abs(SphereRadius - Dist) > Tol)
            Count++;
      }

      if (Count > 0)
         ABORT_ERROR("HorzMeshTest: Edge sphere radius test FAIL");

      // Test lon/lat coordinates of edge centers
      // Convert Cartesian coordinates to lon/lat and check these agree with the
      // values that have been read in
      // Tests that the lon/lat coordinates for edge centers have been read in
      // correctly
      Count = 0;
      for (int Edge = 0; Edge < LocEdges; Edge++) {

         Lon = computeLon(Mesh->XEdgeH(Edge), Mesh->YEdgeH(Edge),
                          Mesh->ZEdgeH(Edge));
         Lat = computeLat(Mesh->XEdgeH(Edge), Mesh->YEdgeH(Edge),
                          Mesh->ZEdgeH(Edge));

         if (abs(Lon - Mesh->LonEdgeH(Edge)) > Tol)
            Count++;
         if (abs(Lat - Mesh->LatEdgeH(Edge)) > Tol)
            Count++;
      }

      if (Count > 0)
         ABORT_ERROR("HorzMeshTest: Edge lon/lat test FAIL");

      // Test sum of local mesh vertices
      // Get the global sum of all local vertex counts
      // Tests that the correct vertex counts have been retrieved from the
      // Decomp object
      I4 SumVertices;
      I4 LocVertices;
      LocVertices = Mesh->NVerticesOwned;
      Err = MPI_Allreduce(&LocVertices, &SumVertices, 1, MPI_INT32_T, MPI_SUM,
                          Comm);
      if (Err != MPI_SUCCESS)
         ABORT_ERROR("HorzMeshTest: MPI error summing vertices");

      if (SumVertices != DefDecomp->NVerticesGlobal)
         ABORT_ERROR("HorzMeshTest: Sum vertex ID test FAIL {} {}", SumVertices,
                     DefDecomp->NVerticesGlobal);

      // Test that vertex coordinates are on sphere
      // Check that all vertices are a uniform distance from the origin
      // Tests that the Cartesian coordinates for vertices have been read in
      // correctly
      SphereRadius =
          distance(Mesh->XVertexH(0), Mesh->YVertexH(0), Mesh->ZVertexH(0));
      Count = 0;
      for (int Vertex = 0; Vertex < LocVertices; Vertex++) {
         Dist = distance(Mesh->XVertexH(Vertex), Mesh->YVertexH(Vertex),
                         Mesh->ZVertexH(Vertex));
         if (abs(SphereRadius - Dist) > Tol)
            Count++;
      }

      if (Count > 0)
         ABORT_ERROR("HorzMeshTest: Vertex sphere radius test FAIL");

      // Test lon/lat coordinates of vertices
      // Convert Cartesian coordinates to lon/lat and check these agree with the
      // values that have been read in
      // Tests that the lon/lat coordinates for vertices have been read in
      // correctly
      Count = 0;
      for (int Vertex = 0; Vertex < LocVertices; Vertex++) {

         Lon = computeLon(Mesh->XVertexH(Vertex), Mesh->YVertexH(Vertex),
                          Mesh->ZVertexH(Vertex));
         Lat = computeLat(Mesh->XVertexH(Vertex), Mesh->YVertexH(Vertex),
                          Mesh->ZVertexH(Vertex));

         if (abs(Lon - Mesh->LonVertexH(Vertex)) > Tol)
            Count++;

         if (abs(Lat - Mesh->LatVertexH(Vertex)) > Tol)
            Count++;
      }

      if (Count > 0)
         ABORT_ERROR("HorzMeshTest: Vertex lon/lat test FAIL");

      // Test cell areas
      // Find the global sum of all the local cell areas
      // and compares to reasonable value for Earth's ocean area
      // Tests that cell areas have been read in correctly
      R8 LocSumArea = 0;
      R8 SumCellArea;
      for (int Cell = 0; Cell < LocCells; Cell++) {
         LocSumArea += Mesh->AreaCellH(Cell);
      }
      Err = MPI_Allreduce(&LocSumArea, &SumCellArea, 1, MPI_DOUBLE, MPI_SUM,
                          Comm);
      if (Err != MPI_SUCCESS)
         ABORT_ERROR("HorzMeshTest: MPI error summing cell area");

      R8 OceanArea = 3.61e14;
      if (abs(SumCellArea - OceanArea) / OceanArea > 0.05)
         ABORT_ERROR("HorzMeshTest: Cell area test FAIL");

      // Test triangle areas
      // Find the global sum of all the local triangle areas
      // and compare to resonable value for the Earth's ocean area
      // Tests that the triangle areas have been read in correctly
      LocSumArea = 0;
      R8 SumTriangleArea;
      for (int Vertex = 0; Vertex < LocVertices; Vertex++) {
         LocSumArea += Mesh->AreaTriangleH(Vertex);
      }
      Err = MPI_Allreduce(&LocSumArea, &SumTriangleArea, 1, MPI_DOUBLE, MPI_SUM,
                          Comm);
      if (Err != MPI_SUCCESS)
         ABORT_ERROR("HorzMeshTest: MPI error summing triangle area");

      if (abs(SumTriangleArea - OceanArea) / OceanArea > 0.05)
         ABORT_ERROR("HorzMeshTest: Triangle area test FAIL");

      // Test kite areas
      // Find the local sum of all the local kite areas
      // and compare to reasonable value for the Earth's ocean area
      // Tests that the kite areas have been read in correctly
      LocSumArea = 0;
      R8 SumKiteArea;
      for (int Vertex = 0; Vertex < LocVertices; Vertex++) {
         for (int i = 0; i < Mesh->VertexDegree; i++) {
            LocSumArea += Mesh->KiteAreasOnVertexH(Vertex, i);
         }
      }
      Err = MPI_Allreduce(&LocSumArea, &SumKiteArea, 1, MPI_DOUBLE, MPI_SUM,
                          Comm);
      if (Err != MPI_SUCCESS)
         ABORT_ERROR("HorzMeshTest: MPI error summing kite area");

      if (abs(SumKiteArea - OceanArea) / OceanArea > 0.05)
         ABORT_ERROR("HorzMeshTest: Kite area test FAIL");

      // Test DcEdge
      // Compute spherical distance between cell centers and compare to value
      // that was read in Tests that the distances between cell centers have
      // been read in correctly
      Count = 0;
      for (int Edge = 0; Edge < LocEdges; Edge++) {
         int Cell1 = Mesh->CellsOnEdgeH(Edge, 0);
         int Cell2 = Mesh->CellsOnEdgeH(Edge, 1);

         if ((Cell1 < DefDecomp->NCellsAll) && (Cell2 < DefDecomp->NCellsAll)) {

            R8 Dc =
                sphereDistance(Mesh->LonCellH(Cell1), Mesh->LatCellH(Cell1),
                               Mesh->LonCellH(Cell2), Mesh->LatCellH(Cell2));
            Dc = SphereRadius * Dc;

            if (abs((Dc - Mesh->DcEdgeH(Edge)) / Mesh->DcEdgeH(Edge)) > Tol) {
               Count++;
            }
         }
      }

      if (Count > 0)
         ABORT_ERROR("HorzMeshTest: DcEdge test FAIL");

      // Test DvEdge
      // Compute spherical distance between vertices on edges and compare to
      // value that was read in Tests that the distances between vertices have
      // been read in correctly
      Count = 0;
      for (int Edge = 0; Edge < LocEdges; Edge++) {
         int Vertex1 = Mesh->VerticesOnEdgeH(Edge, 0);
         int Vertex2 = Mesh->VerticesOnEdgeH(Edge, 1);

         if ((Vertex1 < DefDecomp->NVerticesAll) &&
             (Vertex2 < DefDecomp->NVerticesAll)) {

            R8 Dv = sphereDistance(
                Mesh->LonVertexH(Vertex1), Mesh->LatVertexH(Vertex1),
                Mesh->LonVertexH(Vertex2), Mesh->LatVertexH(Vertex2));

            Dv = SphereRadius * Dv;

            if (abs((Dv - Mesh->DvEdgeH(Edge)) / Mesh->DvEdgeH(Edge)) > Tol) {
               Count++;
            }
         }
      }

      if (Count > 0)
         ABORT_ERROR("HorzMeshTest: DvEdge test FAIL");

      // Test AngleEdge
      // Check that the range of edge angles is between (-Pi, Pi)
      // Tests that the edge angles have been read in correctly
      Count = 0;
      for (int Edge = 0; Edge < LocEdges; Edge++) {
         if (abs(Mesh->AngleEdgeH(Edge)) > Pi) {
            Count++;
         }
      }

      if (Count > 0)
         ABORT_ERROR("HorzMeshTest: angleEdge test FAIL");

      // Test FCell
      // Compute the Coriolis parameter for cell centers and compare with values
      // that were read in
      // Tests that the cell Coriolis values were read in correctly
      Count = 0;
      for (int Cell = 0; Cell < LocCells; Cell++) {
         R8 F = coriolis(Mesh->LatCellH(Cell));

         if (abs(F - Mesh->FCellH(Cell)) > Tol) {
            Count++;
         }
      }

      if (Count > 0)
         ABORT_ERROR("HorzMeshTest: fCell test FAIL");

      // Test FVertex
      // Compute the Coriolis parameter for vertices and compare with values
      // that were read in Tests that the vertex Coriolis values were read in
      // correctly
      Count = 0;
      for (int Vertex = 0; Vertex < LocVertices; Vertex++) {

         R8 F = coriolis(Mesh->LatVertexH(Vertex));

         if (abs(F - Mesh->FVertexH(Vertex)) > Tol) {
            Count++;
         }
      }

      if (Count > 0)
         ABORT_ERROR("HorzMeshTest: fVertex test FAIL");

      // Test FEdge
      // Compute the Coriolis parameter for edges and compare with values that
      // were read in
      // Tests that the edge Coriolis values were read in correctly
      Count = 0;
      for (int Edge = 0; Edge < LocEdges; Edge++) {
         R8 F = coriolis(Mesh->LatEdgeH(Edge));

         if (abs(F - Mesh->FEdgeH(Edge)) > Tol) {
            Count++;
         }
      }

      if (Count > 0)
         ABORT_ERROR("HorzMeshTest: fEdge test FAIL");

      // Test weightsOnEdge
      // Check the range of the edge weights
      // Tests that the edge weights were read in correctly
      Count = 0;
      for (int Edge = 0; Edge < LocEdges; Edge++) {
         for (int i = 0; i < Mesh->MaxEdges2; i++) {
            if (abs(Mesh->WeightsOnEdgeH(Edge, i)) > 1.0) {
               Count++;
            }
         }
      }

      if (Count > 0)
         ABORT_ERROR("HorzMeshTest: weightsOnEdge test FAIL");

      // Test edgeSignOnCell
      // Check that the sign corresponds with convention
      // Tests that the edge sign values were calculated correctly
      Count = 0;
      for (int Edge = 0; Edge < LocEdges; Edge++) {
         int Cell0 = Mesh->CellsOnEdgeH(Edge, 0);
         int IEdge0;
         for (int i = 0; i < Mesh->NEdgesOnCellH(Cell0); i++) {
            if (Mesh->EdgesOnCellH(Cell0, i) == Edge) {
               IEdge0 = i;
               break;
            }
         }
         if (abs(Mesh->EdgeSignOnCellH(Cell0, IEdge0) + 1.0) > Tol) {
            Count++;
         }

         int Cell1 = Mesh->CellsOnEdgeH(Edge, 1);
         if (Cell1 < DefDecomp->NCellsAll) {
            int IEdge1;
            for (int i = 0; i < Mesh->NEdgesOnCellH(Cell1); i++) {
               if (Mesh->EdgesOnCellH(Cell1, i) == Edge) {
                  IEdge1 = i;
                  break;
               }
            }
            if (abs(Mesh->EdgeSignOnCellH(Cell1, IEdge1) - 1.0) > Tol) {
               Count++;
            }
         }
      }

      if (Count > 0)
         ABORT_ERROR("HorzMeshTest: edgeSignOnCell test FAIL");

      // Test edgeSignOnVertex
      // Check that the sign corresponds with convention
      // Tests that the edge sign vlues were calculated correctly
      Count = 0;
      for (int Edge = 0; Edge < LocEdges; Edge++) {
         int Vertex0 = Mesh->VerticesOnEdgeH(Edge, 0);
         int IEdge0;
         for (int i = 0; i < Mesh->VertexDegree; i++) {
            if (Mesh->EdgesOnVertexH(Vertex0, i) == Edge) {
               IEdge0 = i;
               break;
            }
         }
         if (abs(Mesh->EdgeSignOnVertexH(Vertex0, IEdge0) + 1.0) > Tol) {
            Count++;
         }

         int Vertex1 = Mesh->VerticesOnEdgeH(Edge, 1);
         int IEdge1;
         for (int i = 0; i < Mesh->VertexDegree; i++) {
            if (Mesh->EdgesOnVertexH(Vertex1, i) == Edge) {
               IEdge1 = i;
               break;
            }
         }
         if (abs(Mesh->EdgeSignOnVertexH(Vertex1, IEdge1) - 1.0) > Tol) {
            Count++;
         }
      }

      if (Count > 0)
         ABORT_ERROR("HorzMeshTest: edgeSignOnVertex test FAIL");

      // Test cell halo values
      // Perform halo exhange on owned cell only array and compare
      // read values
      // Tests that halo values are read in correctly
      Halo *DefHalo = Halo::getDefault();
      HostArray1DR8 XCellTest("XCellTest", Mesh->NCellsSize);
      deepCopy(XCellTest, Mesh->XCellH);

      for (int Cell = Mesh->NCellsOwned; Cell < Mesh->NCellsAll; Cell++) {
         XCellTest(Cell) = 0.0;
      }
      DefHalo->exchangeFullArrayHalo(XCellTest, OnCell);

      Count = 0;
      for (int Cell = 0; Cell < Mesh->NCellsAll; Cell++) {
         if (Mesh->XCellH(Cell) != XCellTest(Cell)) {
            Count++;
            break;
         }
      }

      if (Count > 0)
         ABORT_ERROR("HorzMeshTest: cell halo exhange FAIL");

      // Test edge halo values
      // Perform halo exhange on owned edge only array and compare
      // read values
      // Tests that halo values are read in correctly
      HostArray1DR8 XEdgeTest("XEdgeTest", Mesh->NEdgesSize);
      deepCopy(XEdgeTest, Mesh->XEdgeH);

      for (int Edge = Mesh->NEdgesOwned; Edge < Mesh->NEdgesAll; Edge++) {
         XEdgeTest(Edge) = 0.0;
      }
      DefHalo->exchangeFullArrayHalo(XEdgeTest, OnEdge);

      Count = 0;
      for (int Edge = 0; Edge < Mesh->NEdgesAll; Edge++) {
         if (Mesh->XEdgeH(Edge) != XEdgeTest(Edge)) {
            Count++;
            break;
         }
      }

      if (Count > 0)
         ABORT_ERROR("HorzMeshTest: edge halo exhange FAIL");

      // Test vertex halo values
      // Perform halo exhange on owned vertex only array and compare
      // read values
      // Tests that halo values are read in correctly
      HostArray1DR8 XVertexTest("XVertexTest", Mesh->NVerticesSize);
      deepCopy(XVertexTest, Mesh->XVertexH);

      for (int Vertex = Mesh->NVerticesOwned; Vertex < Mesh->NVerticesAll;
           Vertex++) {
         XVertexTest(Vertex) = 0.0;
      }
      DefHalo->exchangeFullArrayHalo(XVertexTest, OnVertex);

      Count = 0;
      for (int Vertex = 0; Vertex < Mesh->NVerticesAll; Vertex++) {
         if (Mesh->XVertexH(Vertex) != XVertexTest(Vertex)) {
            Count++;
            break;
         }
      }

      if (Count > 0)
         ABORT_ERROR("HorzMeshTest: vertex halo exhange FAIL");

      // Test that all Cell IDs are accounted for by
      // summing the list of owned values by all tasks. The result should
      // be the sum of the integers from 1 to NCellsGlobal
      OMEGA::I4 RefSumIDs = 0, LocSumIDs = 0, SumIDs = 0;
      OMEGA::HostArray1DI4 CellIDH = OMEGA::createHostMirrorCopy(Mesh->CellID);
      for (int n = 0; n < Mesh->NCellsGlobal; ++n)
         RefSumIDs += n + 1;
      for (int n = 0; n < Mesh->NCellsOwned; ++n)
         LocSumIDs += CellIDH(n);
      Err = MPI_Allreduce(&LocSumIDs, &SumIDs, 1, MPI_INT32_T, MPI_SUM, Comm);
      if (Err != MPI_SUCCESS)
         ABORT_ERROR("HorzMeshTest: MPI error summing cell IDs");

      if (SumIDs != RefSumIDs)
         ABORT_ERROR("DecompTest: Sum cell ID test FAIL {} {}", SumIDs,
                     RefSumIDs);

      // Finalize Omega objects
      IOStream::finalize();
      HorzMesh::clear();
      Field::clear();
      Dimension::clear();
      Halo::clear();
      Decomp::clear();
      MachEnv::removeAll();
   }
   Pacer::finalize();
   Kokkos::finalize();

   // Success if we made it here
   LOG_INFO("------ Horz Mesh unit tests successful ------");
   MPI_Finalize();

   // Return success if we made it here
   return 0;

} // end of main
//===-----------------------------------------------------------------------===/
