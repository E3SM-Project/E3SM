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
#include "Halo.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "mpi.h"

#include <iostream>

//------------------------------------------------------------------------------
// The initialization routine for Mesh testing. It calls various
// init routines, including the creation of the default decomposition.

int initHorzMeshTest() {

   int Err = 0;

   // Initialize the Machine Environment class - this also creates
   // the default MachEnv. Then retrieve the default environment and
   // some needed data members.
   OMEGA::MachEnv::init(MPI_COMM_WORLD);
   OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefault();
   MPI_Comm DefComm       = DefEnv->getComm();

   // Initialize the Logging system
   OMEGA::initLogging(DefEnv);

   // Open config file
   OMEGA::Config("Omega");
   Err = OMEGA::Config::readAll("omega.yml");
   if (Err != 0) {
      LOG_CRITICAL("HorzMeshTest: Error reading config file");
      return Err;
   }

   // Initialize the IO system
   Err = OMEGA::IO::init(DefComm);
   if (Err != 0)
      LOG_ERROR("HorzMeshTest: error initializing parallel IO");

   // Create the default decomposition (initializes the decomposition)
   Err = OMEGA::Decomp::init();
   if (Err != 0)
      LOG_ERROR("HorzMeshTest: error initializing default decomposition");

   // Initialize the default halo
   Err = OMEGA::Halo::init();
   if (Err != 0)
      LOG_ERROR("HorzMeshTest: error initializing default halo");

   // Initialize the default mesh
   Err = OMEGA::HorzMesh::init();
   if (Err != 0)
      LOG_ERROR("HorzMeshTest: error initializing default mesh");

   return Err;
}

//------------------------------------------------------------------------------
// Computes the distance of a x,y,z coordinate from the origin
OMEGA::R8 distance(OMEGA::R8 x, OMEGA::R8 y, OMEGA::R8 z) {

   OMEGA::R8 dist;

   dist = sqrt(x * x + y * y + z * z);

   return dist;
}

//------------------------------------------------------------------------------
// Computes the distance between two lon/lat points on the sphere
OMEGA::R8 sphereDistance(OMEGA::R8 lon1, OMEGA::R8 lat1, OMEGA::R8 lon2,
                         OMEGA::R8 lat2) {

   OMEGA::R8 arg;

   arg = sqrt(pow(sin(0.5 * (lat1 - lat2)), 2) +
              cos(lat2) * cos(lat1) * pow(sin(0.5 * (lon1 - lon2)), 2));
   return 2.0 * asin(arg);
}

//------------------------------------------------------------------------------
// Computes the longitude of a point given its Cartesian coordinates
OMEGA::R8 computeLon(OMEGA::R8 x, OMEGA::R8 y, OMEGA::R8 z) {

   OMEGA::R8 lon;
   lon = atan2(y, x);

   OMEGA::R8 pi;
   pi = 4.0 * atan(1.0);

   if (lon < 0.0) {
      lon = 2.0 * pi + lon;
   }

   return lon;
}

//------------------------------------------------------------------------------
// Computes the latitude of a point given its Cartesian coordinates
OMEGA::R8 computeLat(OMEGA::R8 x, OMEGA::R8 y, OMEGA::R8 z) {

   OMEGA::R8 dist;
   dist = distance(x, y, z);

   OMEGA::R8 lat;
   lat = asin(z / dist);

   return lat;
}

//------------------------------------------------------------------------------
// Computes coriolis parameter for a given latitude
OMEGA::R8 coriolis(OMEGA::R8 lat) {

   OMEGA::R8 f;
   OMEGA::R8 omega = 7.29212e-5;

   f = 2.0 * omega * sin(lat);

   return f;
}

//------------------------------------------------------------------------------
// The test driver for Mesh-> This tests the decomposition of a sample
// horizontal domain and verifies the mesh is read in correctly.
//
int main(int argc, char *argv[]) {

   int RetVal = 0;

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   {

      OMEGA::R8 tol = 1e-6;
      OMEGA::R8 pi  = 4.0 * atan(1.0);

      // Call initialization routine to create the default decomposition
      int Err = initHorzMeshTest();
      if (Err != 0)
         LOG_CRITICAL("HorzMeshTest: Error initializing");

      // Get MPI vars if needed
      OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefault();
      MPI_Comm Comm          = DefEnv->getComm();
      OMEGA::I4 MyTask       = DefEnv->getMyTask();
      OMEGA::I4 NumTasks     = DefEnv->getNumTasks();
      bool IsMaster          = DefEnv->isMasterTask();

      // Test retrieval of the default decomposition
      OMEGA::Decomp *DefDecomp = OMEGA::Decomp::getDefault();
      if (DefDecomp) { // true if non-null ptr
         LOG_INFO("HorzMeshTest: Default decomp retrieval PASS");
      } else {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: Default decomp retrieval FAIL");
      }

      // Retrieve default mesh
      OMEGA::HorzMesh *Mesh = OMEGA::HorzMesh::getDefault();

      // Test sum of local mesh cells
      // Get the global sum of all local cell counts
      // Tests that the correct cell counts have been retrieved from the Decomp
      // object
      OMEGA::I4 SumCells;
      OMEGA::I4 LocCells;
      LocCells = Mesh->NCellsOwned;
      Err = MPI_Allreduce(&LocCells, &SumCells, 1, MPI_INT32_T, MPI_SUM, Comm);

      if (SumCells == DefDecomp->NCellsGlobal) {
         LOG_INFO("HorzMeshTest: Sum cell ID test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: Sum cell ID test FAIL {} {}", SumCells,
                  DefDecomp->NCellsGlobal);
      }

      // Test that cell centers are on sphere
      // Check that all cell centers are a uniform distance from the origin
      // Tests that the Cartesian coordinates for cell centers have been read in
      // corectly
      OMEGA::R8 sphere_radius =
          distance(Mesh->XCellH(0), Mesh->YCellH(0), Mesh->ZCellH(0));
      OMEGA::R8 dist;
      OMEGA::I4 count = 0;
      for (int Cell = 0; Cell < LocCells; Cell++) {
         dist = distance(Mesh->XCellH(Cell), Mesh->YCellH(Cell),
                         Mesh->ZCellH(Cell));
         if (abs(sphere_radius - dist) > tol)
            count++;
      }

      if (count > 0) {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: Cell sphere radius test FAIL");
      } else {
         LOG_INFO("HorzMeshTest: Cell sphere radius test PASS");
      }

      // Test lon/lat coordinates of cell centers
      // Convert Cartesian coordinates to lon/lat and check these agree with the
      // values that have been read in
      // Tests that the lon/lat coordinates for cell
      // centers have been read in correctly
      OMEGA::R8 lon;
      OMEGA::R8 lat;
      count = 0;
      for (int Cell = 0; Cell < LocCells; Cell++) {

         lon = computeLon(Mesh->XCellH(Cell), Mesh->YCellH(Cell),
                          Mesh->ZCellH(Cell));
         lat = computeLat(Mesh->XCellH(Cell), Mesh->YCellH(Cell),
                          Mesh->ZCellH(Cell));

         if (abs(lon - Mesh->LonCellH(Cell)) > tol)
            count++;
         if (abs(lat - Mesh->LatCellH(Cell)) > tol)
            count++;
      }

      if (count > 0) {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: Cell lon/lat test FAIL");
      } else {
         LOG_INFO("HorzMeshTest: Cell lon/lat test PASS");
      }

      // Test sum of local mesh edges
      // Get the global sum of all local edge counts
      // Tests that the correct edge counts have been retrieved from the Decomp
      // object
      OMEGA::I4 SumEdges;
      OMEGA::I4 LocEdges;
      LocEdges = Mesh->NEdgesOwned;
      Err = MPI_Allreduce(&LocEdges, &SumEdges, 1, MPI_INT32_T, MPI_SUM, Comm);

      if (SumEdges == DefDecomp->NEdgesGlobal) {
         LOG_INFO("HorzMeshTest: Sum edge ID test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: Sum edge ID test FAIL {} {}", SumEdges,
                  DefDecomp->NEdgesGlobal);
      }

      // Test that edge coordinates are on sphere
      // Check that all edge centers are a uniform distance from the origin
      // Tests that the Cartesian coordinates for edge centers have been read in
      // correctly
      sphere_radius =
          distance(Mesh->XEdgeH(0), Mesh->YEdgeH(0), Mesh->ZEdgeH(0));
      count = 0;
      for (int Edge = 0; Edge < LocEdges; Edge++) {
         dist = distance(Mesh->XEdgeH(Edge), Mesh->YEdgeH(Edge),
                         Mesh->ZEdgeH(Edge));
         if (abs(sphere_radius - dist) > tol)
            count++;
      }

      if (count > 0) {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: Edge sphere radius test FAIL");
      } else {
         LOG_INFO("HorzMeshTest: Edge sphere radius test PASS");
      }

      // Test lon/lat coordinates of edge centers
      // Convert Cartesian coordinates to lon/lat and check these agree with the
      // values that have been read in
      // Tests that the lon/lat coordinates for edge centers have been read in
      // correctly
      count = 0;
      for (int Edge = 0; Edge < LocEdges; Edge++) {

         lon = computeLon(Mesh->XEdgeH(Edge), Mesh->YEdgeH(Edge),
                          Mesh->ZEdgeH(Edge));
         lat = computeLat(Mesh->XEdgeH(Edge), Mesh->YEdgeH(Edge),
                          Mesh->ZEdgeH(Edge));

         if (abs(lon - Mesh->LonEdgeH(Edge)) > tol)
            count++;
         if (abs(lat - Mesh->LatEdgeH(Edge)) > tol)
            count++;
      }

      if (count > 0) {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: Edge lon/lat test FAIL");
      } else {
         LOG_INFO("HorzMeshTest: Edge lon/lat test PASS");
      }

      // Test sum of local mesh vertices
      // Get the global sum of all local vertex counts
      // Tests that the correct vertex counts have been retrieved from the
      // Decomp object
      OMEGA::I4 SumVertices;
      OMEGA::I4 LocVertices;
      LocVertices = Mesh->NVerticesOwned;
      Err = MPI_Allreduce(&LocVertices, &SumVertices, 1, MPI_INT32_T, MPI_SUM,
                          Comm);

      if (SumVertices == DefDecomp->NVerticesGlobal) {
         LOG_INFO("HorzMeshTest: Sum vertex ID test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: Sum vertex ID test FAIL {} {}", SumVertices,
                  DefDecomp->NVerticesGlobal);
      }

      // Test that vertex coordinates are on sphere
      // Check that all vertices are a uniform distance from the origin
      // Tests that the Cartesian coordinates for vertices have been read in
      // correctly
      sphere_radius =
          distance(Mesh->XVertexH(0), Mesh->YVertexH(0), Mesh->ZVertexH(0));
      count = 0;
      for (int Vertex = 0; Vertex < LocVertices; Vertex++) {
         dist = distance(Mesh->XVertexH(Vertex), Mesh->YVertexH(Vertex),
                         Mesh->ZVertexH(Vertex));
         if (abs(sphere_radius - dist) > tol)
            count++;
      }

      if (count > 0) {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: Vertex sphere radius test FAIL");
      } else {
         LOG_INFO("HorzMeshTest: Vertex sphere radius test PASS");
      }

      // Test lon/lat coordinates of vertices
      // Convert Cartesian coordinates to lon/lat and check these agree with the
      // values that have been read in
      // Tests that the lon/lat coordinates for vertices have been read in
      // correctly
      count = 0;
      for (int Vertex = 0; Vertex < LocVertices; Vertex++) {

         lon = computeLon(Mesh->XVertexH(Vertex), Mesh->YVertexH(Vertex),
                          Mesh->ZVertexH(Vertex));
         lat = computeLat(Mesh->XVertexH(Vertex), Mesh->YVertexH(Vertex),
                          Mesh->ZVertexH(Vertex));

         if (abs(lon - Mesh->LonVertexH(Vertex)) > tol)
            count++;

         if (abs(lat - Mesh->LatVertexH(Vertex)) > tol)
            count++;
      }

      if (count > 0) {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: Vertex lon/lat test FAIL");
      } else {
         LOG_INFO("HorzMeshTest: Vertex lon/lat test PASS");
      }

      // Test bounds of bathymetry
      // Find minimum and maximum values of the bottom depth
      // and compares to reasonable values
      // Tests that the bottom depth has been read in correctly
      OMEGA::R8 MaxBathy = -1e10;
      OMEGA::R8 MinBathy = 1e10;
      for (int Cell = 0; Cell < LocCells; Cell++) {
         if (Mesh->BottomDepthH(Cell) < MinBathy) {
            MinBathy = Mesh->BottomDepthH(Cell);
         }
         if (Mesh->BottomDepthH(Cell) > MaxBathy) {
            MaxBathy = Mesh->BottomDepthH(Cell);
         }
      }

      if ((MinBathy > 0) && (MaxBathy < 11000.0)) {
         LOG_INFO("HorzMeshTest: Bathy min/max test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: Bathy min/max test FAIL");
      }

      // Test cell areas
      // Find the global sum of all the local cell areas
      // and compares to reasonable value for Earth's ocean area
      // Tests that cell areas have been read in correctly
      OMEGA::R8 LocSumArea = 0;
      OMEGA::R8 SumCellArea;
      for (int Cell = 0; Cell < LocCells; Cell++) {
         LocSumArea += Mesh->AreaCellH(Cell);
      }
      Err = MPI_Allreduce(&LocSumArea, &SumCellArea, 1, MPI_DOUBLE, MPI_SUM,
                          Comm);

      OMEGA::R8 OceanArea = 3.61e14;
      if (abs(SumCellArea - OceanArea) / OceanArea < 0.05) {
         LOG_INFO("HorzMeshTest: Cell area test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: Cell area test FAIL");
      }

      // Test triangle areas
      // Find the global sum of all the local triangle areas
      // and compare to resonable value for the Earth's ocean area
      // Tests that the triangle areas have been read in correctly
      LocSumArea = 0;
      OMEGA::R8 SumTriangleArea;
      for (int Vertex = 0; Vertex < LocVertices; Vertex++) {
         LocSumArea += Mesh->AreaTriangleH(Vertex);
      }
      Err = MPI_Allreduce(&LocSumArea, &SumTriangleArea, 1, MPI_DOUBLE, MPI_SUM,
                          Comm);

      if (abs(SumTriangleArea - OceanArea) / OceanArea < 0.05) {
         LOG_INFO("HorzMeshTest: Triangle area test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: Triangle area test FAIL");
      }

      // Test kite areas
      // Find the local sum of all the local kite areas
      // and compare to reasonable value for the Earth's ocean area
      // Tests that the kite areas have been read in correctly
      LocSumArea = 0;
      OMEGA::R8 SumKiteArea;
      for (int Vertex = 0; Vertex < LocVertices; Vertex++) {
         for (int i = 0; i < Mesh->VertexDegree; i++) {
            LocSumArea += Mesh->KiteAreasOnVertexH(Vertex, i);
         }
      }
      Err = MPI_Allreduce(&LocSumArea, &SumKiteArea, 1, MPI_DOUBLE, MPI_SUM,
                          Comm);

      if (abs(SumKiteArea - OceanArea) / OceanArea < 0.05) {
         LOG_INFO("HorzMeshTest: Kite area test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: Kite area test FAIL");
      }

      // Test dcEdge
      // Compute spherical distance between cell centers and compare to value
      // that was read in Tests that the distances between cell centers have
      // been read in correctly
      count = 0;
      for (int Edge = 0; Edge < LocEdges; Edge++) {
         int Cell1 = Mesh->CellsOnEdgeH(Edge, 0);
         int Cell2 = Mesh->CellsOnEdgeH(Edge, 1);

         if ((Cell1 < DefDecomp->NCellsAll) && (Cell2 < DefDecomp->NCellsAll)) {

            OMEGA::R8 dc =
                sphereDistance(Mesh->LonCellH(Cell1), Mesh->LatCellH(Cell1),
                               Mesh->LonCellH(Cell2), Mesh->LatCellH(Cell2));
            dc = sphere_radius * dc;

            if (abs((dc - Mesh->DcEdgeH(Edge)) / Mesh->DcEdgeH(Edge)) > tol) {
               count++;
            }
         }
      }

      if (count == 0) {
         LOG_INFO("HorzMeshTest: dcEdge test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: dcEdge test FAIL");
      }

      // Test dvEdge
      // Compute spherical distance between vertices on edges and compare to
      // value that was read in Tests that the distances between vertices have
      // been read in correctly
      count = 0;
      for (int Edge = 0; Edge < LocEdges; Edge++) {
         int Vertex1 = Mesh->VerticesOnEdgeH(Edge, 0);
         int Vertex2 = Mesh->VerticesOnEdgeH(Edge, 1);

         if ((Vertex1 < DefDecomp->NVerticesAll) &&
             (Vertex2 < DefDecomp->NVerticesAll)) {

            OMEGA::R8 dv = sphereDistance(
                Mesh->LonVertexH(Vertex1), Mesh->LatVertexH(Vertex1),
                Mesh->LonVertexH(Vertex2), Mesh->LatVertexH(Vertex2));

            dv = sphere_radius * dv;

            if (abs((dv - Mesh->DvEdgeH(Edge)) / Mesh->DvEdgeH(Edge)) > tol) {
               count++;
            }
         }
      }

      if (count == 0) {
         LOG_INFO("HorzMeshTest: dvEdge test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: dvEdge test FAIL");
      }

      // Test angleEdge
      // Check that the range of edge angles is between (-pi, pi)
      // Tests that the edge angles have been read in correctly
      count = 0;
      for (int Edge = 0; Edge < LocEdges; Edge++) {
         if (abs(Mesh->AngleEdgeH(Edge)) > pi) {
            count++;
         }
      }

      if (count == 0) {
         LOG_INFO("HorzMeshTest: angleEdge test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: angleEdge test FAIL");
      }

      // Test fCell
      // Compute the Coriolis parameter for cell centers and compare with values
      // that were read in
      // Tests that the cell Coriolis values were read in correctly
      count = 0;
      for (int Cell = 0; Cell < LocCells; Cell++) {
         OMEGA::R8 f = coriolis(Mesh->LatCellH(Cell));

         if (abs(f - Mesh->FCellH(Cell)) > tol) {
            count++;
         }
      }

      if (count == 0) {
         LOG_INFO("HorzMeshTest: fCell test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: fCell test FAIL");
      }

      // Test fVertex
      // Compute the Coriolis parameter for vertices and compare with values
      // that were read in Tests that the vertex Coriolis values were read in
      // correctly
      count = 0;
      for (int Vertex = 0; Vertex < LocVertices; Vertex++) {

         OMEGA::R8 f = coriolis(Mesh->LatVertexH(Vertex));

         if (abs(f - Mesh->FVertexH(Vertex)) > tol) {
            count++;
         }
      }

      if (count == 0) {
         LOG_INFO("HorzMeshTest: fVertex test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: fVertex test FAIL");
      }

      // Test fEdge
      // Compute the Coriolis parameter for edges and compare with values that
      // were read in
      // Tests that the edge Coriolis values were read in correctly
      count = 0;
      for (int Edge = 0; Edge < LocEdges; Edge++) {
         OMEGA::R8 f = coriolis(Mesh->LatEdgeH(Edge));

         if (abs(f - Mesh->FEdgeH(Edge)) > tol) {
            count++;
         }
      }

      if (count == 0) {
         LOG_INFO("HorzMeshTest: fEdge test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: fEdge test FAIL");
      }

      // Test weightsOnEdge
      // Check the range of the edge weights
      // Tests that the edge weights were read in correctly
      count = 0;
      for (int Edge = 0; Edge < LocEdges; Edge++) {
         for (int i = 0; i < Mesh->MaxEdges2; i++) {
            if (abs(Mesh->WeightsOnEdgeH(Edge, i)) > 1.0) {
               count++;
            }
         }
      }

      if (count == 0) {
         LOG_INFO("HorzMeshTest: weightsOnEdge test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: weightsOnEdge test FAIL");
      }
      // Test edgeSignOnCell
      // Check that the sign corresponds with convention
      // Tests that the edge sign values were calculated correctly
      count = 0;
      for (int Edge = 0; Edge < LocEdges; Edge++) {
         int Cell0 = Mesh->CellsOnEdgeH(Edge, 0);
         int iEdge0;
         for (int i = 0; i < Mesh->NEdgesOnCellH(Cell0); i++) {
            if (Mesh->EdgesOnCellH(Cell0, i) == Edge) {
               iEdge0 = i;
               break;
            }
         }
         if (abs(Mesh->EdgeSignOnCellH(Cell0, iEdge0) + 1.0) > tol) {
            count++;
         }

         int Cell1 = Mesh->CellsOnEdgeH(Edge, 1);
         if (Cell1 < DefDecomp->NCellsAll) {
            int iEdge1;
            for (int i = 0; i < Mesh->NEdgesOnCellH(Cell1); i++) {
               if (Mesh->EdgesOnCellH(Cell1, i) == Edge) {
                  iEdge1 = i;
                  break;
               }
            }
            if (abs(Mesh->EdgeSignOnCellH(Cell1, iEdge1) - 1.0) > tol) {
               count++;
            }
         }
      }

      if (count == 0) {
         LOG_INFO("HorzMeshTest: edgeSignOnCell test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: edgeSignOnCell test FAIL");
      }

      // Test edgeSignOnVertex
      // Check that the sign corresponds with convention
      // Tests that the edge sign vlues were calculated correctly
      count = 0;
      for (int Edge = 0; Edge < LocEdges; Edge++) {
         int Vertex0 = Mesh->VerticesOnEdgeH(Edge, 0);
         int iEdge0;
         for (int i = 0; i < Mesh->VertexDegree; i++) {
            if (Mesh->EdgesOnVertexH(Vertex0, i) == Edge) {
               iEdge0 = i;
               break;
            }
         }
         if (abs(Mesh->EdgeSignOnVertexH(Vertex0, iEdge0) + 1.0) > tol) {
            count++;
         }

         int Vertex1 = Mesh->VerticesOnEdgeH(Edge, 1);
         int iEdge1;
         for (int i = 0; i < Mesh->VertexDegree; i++) {
            if (Mesh->EdgesOnVertexH(Vertex1, i) == Edge) {
               iEdge1 = i;
               break;
            }
         }
         if (abs(Mesh->EdgeSignOnVertexH(Vertex1, iEdge1) - 1.0) > tol) {
            count++;
         }
      }

      if (count == 0) {
         LOG_INFO("HorzMeshTest: edgeSignOnVertex test PASS");
      } else {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: edgeSignOnVertex test FAIL");
      }

      // Test cell halo values
      // Perform halo exhange on owned cell only array and compare
      // read values
      // Tests that halo values are read in correctly
      OMEGA::Halo *DefHalo = OMEGA::Halo::getDefault();
      OMEGA::HostArray1DR8 XCellTest("XCellTest", Mesh->NCellsSize);
      // Mesh->XCellH.deep_copy_to(XCellTest);
      OMEGA::deepCopy(XCellTest, Mesh->XCellH);

      for (int Cell = Mesh->NCellsOwned; Cell < Mesh->NCellsAll; Cell++) {
         XCellTest(Cell) = 0.0;
      }
      DefHalo->exchangeFullArrayHalo(XCellTest, OMEGA::OnCell);

      count = 0;
      for (int Cell = 0; Cell < Mesh->NCellsAll; Cell++) {
         if (Mesh->XCellH(Cell) != XCellTest(Cell)) {
            count++;
            break;
         }
      }

      if (count == 0) {
         LOG_INFO("HorzMeshTest: cell halo exhange PASS");
      } else {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: cell halo exhange FAIL");
      }

      // Test edge halo values
      // Perform halo exhange on owned edge only array and compare
      // read values
      // Tests that halo values are read in correctly
      OMEGA::HostArray1DR8 XEdgeTest("XEdgeTest", Mesh->NEdgesSize);
      // Mesh->XEdgeH.deep_copy_to(XEdgeTest);
      OMEGA::deepCopy(XEdgeTest, Mesh->XEdgeH);

      for (int Edge = Mesh->NEdgesOwned; Edge < Mesh->NEdgesAll; Edge++) {
         XEdgeTest(Edge) = 0.0;
      }
      DefHalo->exchangeFullArrayHalo(XEdgeTest, OMEGA::OnEdge);

      count = 0;
      for (int Edge = 0; Edge < Mesh->NEdgesAll; Edge++) {
         if (Mesh->XEdgeH(Edge) != XEdgeTest(Edge)) {
            count++;
            break;
         }
      }

      if (count == 0) {
         LOG_INFO("HorzMeshTest: edge halo exhange PASS");
      } else {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: edge halo exhange FAIL");
      }

      // Test vertex halo values
      // Perform halo exhange on owned vertex only array and compare
      // read values
      // Tests that halo values are read in correctly
      OMEGA::HostArray1DR8 XVertexTest("XVertexTest", Mesh->NVerticesSize);
      // Mesh->XVertexH.deep_copy_to(XVertexTest);
      OMEGA::deepCopy(XVertexTest, Mesh->XVertexH);

      for (int Vertex = Mesh->NVerticesOwned; Vertex < Mesh->NVerticesAll;
           Vertex++) {
         XVertexTest(Vertex) = 0.0;
      }
      DefHalo->exchangeFullArrayHalo(XVertexTest, OMEGA::OnVertex);

      count = 0;
      for (int Vertex = 0; Vertex < Mesh->NVerticesAll; Vertex++) {
         if (Mesh->XVertexH(Vertex) != XVertexTest(Vertex)) {
            count++;
            break;
         }
      }

      if (count == 0) {
         LOG_INFO("HorzMeshTest: vertex halo exhange PASS");
      } else {
         RetVal += 1;
         LOG_INFO("HorzMeshTest: vertex halo exhange FAIL");
      }
      // Finalize Omega objects
      OMEGA::HorzMesh::clear();
      OMEGA::Dimension::clear();
      OMEGA::Halo::clear();
      OMEGA::Decomp::clear();
      OMEGA::MachEnv::removeAll();

      // MPI_Status status;
      if (Err == 0)
         LOG_INFO("HorzMeshTest: Successful completion");
   }
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
