//===-- Test driver for OMEGA Decomp -----------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA domain decomposition (Decomp)
///
/// This driver tests the OMEGA domain decomposition, decomposing the
/// horizontal domain and creating a number of index-space arrays for
/// locating and describing mesh locations within a parallel distributed
/// memory.
///
//
//===-----------------------------------------------------------------------===/

#include "Decomp.h"
#include "HorzMesh.h"
#include "DataTypes.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "mpi.h"


#include <iostream>

//------------------------------------------------------------------------------
// The initialization routine for Decomp testing. It calls various
// init routines, including the creation of the default decomposition.

int initHorzMeshTest() {

   int Err = 0;

   // Initialize the Machine Environment class - this also creates
   // the default MachEnv. Then retrieve the default environment and
   // some needed data members.
   OMEGA::MachEnv::init(MPI_COMM_WORLD);
   OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefaultEnv();
   MPI_Comm DefComm       = DefEnv->getComm();

   // Initialize the IO system
   Err = OMEGA::IOInit(DefComm);
   if (Err != 0)
      LOG_ERROR("HorzMeshTest: error initializing parallel IO");

   // Create the default decomposition (initializes the decomposition)
   Err = OMEGA::Decomp::init();
   if (Err != 0)
      LOG_ERROR("HorzMeshTest: error initializing default decomposition");

   return Err;
}

OMEGA::R8 distance(OMEGA::R8 x, OMEGA::R8 y, OMEGA::R8 z) {
   
   OMEGA::R8 dist;

   dist = sqrt(x*x + y*y + z*z);

   return dist;

}

OMEGA::R8 sphereDistance(OMEGA::R8 lon1, OMEGA::R8 lat1, OMEGA::R8 lon2, OMEGA::R8 lat2) {

   OMEGA::R8 arg;

   arg = sqrt( pow(sin(0.5*(lat1-lat2)),2) +
                    cos(lat2)*cos(lat1)*pow(sin(0.5*(lon1-lon2)),2));
   return 2.0*asin(arg);

}

OMEGA::R8 computeLon(OMEGA::R8 x, OMEGA::R8 y, OMEGA::R8 z) {

   OMEGA::R8 lon;
   lon = atan2(y,x);

   OMEGA::R8 pi;
   pi = 4.0*atan(1.0);

   if (lon < 0.0) {
     lon = 2.0*pi + lon;
   }

   return lon;
}

OMEGA::R8 computeLat(OMEGA::R8 x, OMEGA::R8 y, OMEGA::R8 z) {

   OMEGA::R8 dist;
   dist = distance(x, y, z);
   
   OMEGA::R8 lat;
   lat = asin(z/dist);

   return lat;
}

OMEGA::R8 coriolis(OMEGA::R8 lat) {

   OMEGA::R8 f;
   OMEGA::R8 omega = 7.29212e-5;

   f = 2.0*omega*sin(lat);

   return f;

}

//------------------------------------------------------------------------------
// The test driver for Decomp. This tests the decomposition of a sample
// horizontal domain and verifies the mesh is decomposed correctly.
//
int main(int argc, char *argv[]) {

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);
   yakl::init();

   OMEGA::R8 tol = 1e-6;
   OMEGA::R8 pi = 4.0*atan(1.0);

   // Call initialization routine to create the default decomposition
   int Err = initHorzMeshTest();
   if (Err != 0)
      LOG_CRITICAL("HorzMeshTest: Error initializing");

   // Get MPI vars if needed
   OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefaultEnv();
   MPI_Comm Comm          = DefEnv->getComm();
   OMEGA::I4 MyTask       = DefEnv->getMyTask();
   OMEGA::I4 NumTasks     = DefEnv->getNumTasks();
   bool IsMaster          = DefEnv->isMasterTask();

   // Test retrieval of the default decomposition
   OMEGA::Decomp *DefDecomp = OMEGA::Decomp::getDefault();
   if (DefDecomp) { // true if non-null ptr
      LOG_INFO("HorzMeshTest: Default decomp retrieval PASS");
   } else {
      LOG_INFO("HorzMeshTest: Default decomp retrieval FAIL");
      return -1;
   }

   // Initialize mesh
   OMEGA::HorzMesh Mesh(DefDecomp);

   // Test sum of local mesh cells
   OMEGA::I4 SumCells;
   OMEGA::I4 LocCells;
   LocCells = Mesh.NCellsOwned;
   Err = MPI_Allreduce(&LocCells, &SumCells, 1, MPI_INT32_T, MPI_SUM, Comm);


   if (SumCells == DefDecomp->NCellsGlobal){
     LOG_INFO("HorzMeshTest: Sum cell ID test PASS");
   } else {
      LOG_INFO("HorzMeshTest: Sum cell ID test FAIL {} {}", SumCells,
              DefDecomp->NCellsGlobal);
   }

   // Test that cell centers are on sphere
   OMEGA::R8 sphere_radius = distance(Mesh.XCellH(0), Mesh.YCellH(0), Mesh.ZCellH(0));
   OMEGA::R8 dist;
   OMEGA::I4 count = 0;
   for (int Cell = 0; Cell < LocCells; Cell++) {
       dist = distance(Mesh.XCellH(Cell), Mesh.YCellH(Cell), Mesh.ZCellH(Cell));
       if ( abs(sphere_radius - dist) > tol)
          count++;
   }

   if (count > 0) {
     LOG_INFO("HorzMeshTest: Cell sphere radius test FAIL");
   } else {
     LOG_INFO("HorzMeshTest: Cell sphere radius test PASS");
   }

   // Test lon/lat coordinates
   OMEGA::R8 lon;
   OMEGA::R8 lat;
   count = 0;
   for (int Cell = 0; Cell < LocCells; Cell++) {

      lon = computeLon(Mesh.XCellH(Cell), Mesh.YCellH(Cell), Mesh.ZCellH(Cell));
      lat = computeLat(Mesh.XCellH(Cell), Mesh.YCellH(Cell), Mesh.ZCellH(Cell));

      if (abs(lon - Mesh.LonCellH(Cell)) > tol)
         count++;

      if (abs(lat - Mesh.LatCellH(Cell)) > tol)
         count++ ;

   }

   if (count > 0) {
     LOG_INFO("HorzMeshTest: Cell lon/lat test FAIL");
     return -1;
   } else {
     LOG_INFO("HorzMeshTest: Cell lon/lat test PASS");
   }

   // Test sum of local mesh edges 
   OMEGA::I4 SumEdges;
   OMEGA::I4 LocEdges;
   LocEdges = Mesh.NEdgesOwned;
   Err = MPI_Allreduce(&LocEdges, &SumEdges, 1, MPI_INT32_T, MPI_SUM, Comm);


   if (SumEdges == DefDecomp->NEdgesGlobal){
     LOG_INFO("HorzMeshTest: Sum edge ID test PASS");
   } else {
      LOG_INFO("HorzMeshTest: Sum edge ID test FAIL {} {}", SumEdges,
              DefDecomp->NEdgesGlobal);
      return -1;
   }

   // Test that edge coordinates are on sphere
   sphere_radius = distance(Mesh.XEdgeH(0), Mesh.YEdgeH(0), Mesh.ZEdgeH(0));
   count = 0;
   for (int Edge = 0; Edge < LocEdges; Edge++) {
       dist = distance(Mesh.XEdgeH(Edge), Mesh.YEdgeH(Edge), Mesh.ZEdgeH(Edge));
       if ( abs(sphere_radius - dist) > tol)
          count++;
   }

   if (count > 0) {
     LOG_INFO("HorzMeshTest: Edge sphere radius test FAIL");
     return -1;
   } else {
     LOG_INFO("HorzMeshTest: Edge sphere radius test PASS");
   }

   // Test sum of local mesh vertices 
   OMEGA::I4 SumVertices;
   OMEGA::I4 LocVertices;
   LocVertices = Mesh.NVerticesOwned;
   Err = MPI_Allreduce(&LocVertices, &SumVertices, 1, MPI_INT32_T, MPI_SUM, Comm);


   if (SumVertices == DefDecomp->NVerticesGlobal){
     LOG_INFO("HorzMeshTest: Sum vertex ID test PASS");
   } else {
      LOG_INFO("HorzMeshTest: Sum vertex ID test FAIL {} {}", SumVertices,
              DefDecomp->NVerticesGlobal);
      return -1;
   }

   // Test that cell centers are on sphere
   sphere_radius = distance(Mesh.XVertexH(0), Mesh.YVertexH(0), Mesh.ZVertexH(0));
   count = 0;
   for (int Vertex = 0; Vertex < LocVertices; Vertex++) {
       dist = distance(Mesh.XVertexH(Vertex), Mesh.YVertexH(Vertex), Mesh.ZVertexH(Vertex));
       if ( abs(sphere_radius - dist) > tol)
          count++;
   }

   if (count > 0) {
     LOG_INFO("HorzMeshTest: Vertex sphere radius test FAIL");
     return -1;
   } else {
     LOG_INFO("HorzMeshTest: Vertex sphere radius test PASS");
   }

   // Test bounds of bathymetry
   OMEGA::R8 MaxBathy = -1e10;
   OMEGA::R8 MinBathy = 1e10;
   for (int Cell = 0; Cell < LocCells; Cell++) {
      if (Mesh.BottomDepthH(Cell) < MinBathy) {
         MinBathy = Mesh.BottomDepthH(Cell);
      }
      if (Mesh.BottomDepthH(Cell) > MaxBathy) {
         MaxBathy = Mesh.BottomDepthH(Cell);
      }
   }

   if ((MinBathy > 0) && (MaxBathy < 11000.0)) {
      LOG_INFO("HorzMeshTest: Bathy min/max test PASS");
   } else {
      LOG_INFO("HorzMeshTest: Bathy min/max test FAIL");
      return -1;
   }

   // Test cell areas
   OMEGA::R8 LocSumArea = 0;
   OMEGA::R8 SumCellArea;
   for (int Cell = 0; Cell < LocCells; Cell++) {
      LocSumArea += Mesh.AreaCellH(Cell);
   }
   Err = MPI_Allreduce(&LocSumArea, &SumCellArea, 1, MPI_DOUBLE, MPI_SUM, Comm);

   OMEGA::R8 OceanArea = 3.61e14;
   if (abs(SumCellArea - OceanArea)/OceanArea < 0.05) { 
      LOG_INFO("HorzMeshTest: Cell area test PASS");
   } else {
      LOG_INFO("HorzMeshTest: Cell area test FAIL");
      return -1;
   }

   // Test triangle areas
   LocSumArea = 0;
   OMEGA::R8 SumTriangleArea;
   for (int Vertex = 0; Vertex < LocVertices; Vertex++) {
      LocSumArea += Mesh.AreaTriangleH(Vertex);
   }
   Err = MPI_Allreduce(&LocSumArea, &SumTriangleArea, 1, MPI_DOUBLE, MPI_SUM, Comm);

   if (abs(SumTriangleArea - OceanArea)/OceanArea < 0.05) { 
      LOG_INFO("HorzMeshTest: Triangle area test PASS");
   } else {
      LOG_INFO("HorzMeshTest: Triangle area test FAIL");
      return -1;
   }

   // Test dcEdge
   count = 0;
   for (int Edge = 0; Edge < LocEdges; Edge++) {
      int Cell1 = Mesh.CellsOnEdgeH(Edge,0);
      int Cell2 = Mesh.CellsOnEdgeH(Edge,1);

      if ((Cell1<DefDecomp->NCellsAll) && (Cell2<DefDecomp->NCellsAll)) {

         OMEGA::R8 dc = sphere_radius*sphereDistance(Mesh.LonCellH(Cell1), Mesh.LatCellH(Cell1),
                                                     Mesh.LonCellH(Cell2), Mesh.LatCellH(Cell2)); 

         if (abs((dc - Mesh.DcEdgeH(Edge))/Mesh.DcEdgeH(Edge)) > tol) {
            count++;
         }
      }

   }

   if ( count == 0 ) {
      LOG_INFO("HorzMeshTest: dcEdge test PASS");
   } else {
      LOG_INFO("HorzMeshTest: dcEdge test FAIL");
      return -1;
   }
   
   // Test dvEdge
   count = 0;
   for (int Edge = 0; Edge < LocEdges; Edge++) {
      int Vertex1 = Mesh.VerticesOnEdgeH(Edge,0);
      int Vertex2 = Mesh.VerticesOnEdgeH(Edge,1);

      if ((Vertex1<DefDecomp->NVerticesAll) && (Vertex2<DefDecomp->NVerticesAll)) {

         OMEGA::R8 dv = sphere_radius*sphereDistance(Mesh.LonVertexH(Vertex1), Mesh.LatVertexH(Vertex1),
                                                     Mesh.LonVertexH(Vertex2), Mesh.LatVertexH(Vertex2)); 

         if (abs((dv - Mesh.DvEdgeH(Edge))/Mesh.DvEdgeH(Edge)) > tol) {
            count++;
         }
      }

   }

   if ( count == 0 ) {
      LOG_INFO("HorzMeshTest: dvEdge test PASS");
   } else {
      LOG_INFO("HorzMeshTest: dvEdge test FAIL");
      return -1;
   }

   // Test fCell
   count = 0;
   for (int Cell = 0; Cell < LocCells; Cell++) {
      OMEGA::R8 f = coriolis(Mesh.LatCellH(Cell));

      if (abs(f-Mesh.FCellH(Cell)) > tol) {
         count++;
      }
   }

   if ( count == 0 ) {
      LOG_INFO("HorzMeshTest: fCell test PASS");
   } else {
      LOG_INFO("HorzMeshTest: fCell test FAIL");
      return -1;
   }

   // Test fVertex
   count = 0;
   for (int Vertex = 0; Vertex < LocVertices; Vertex++) {
      
      OMEGA::R8 f = coriolis(Mesh.LatVertexH(Vertex));

      if (abs(f-Mesh.FVertexH(Vertex)) > tol) {
         count++;
      }

   }

   if ( count == 0 ) {
      LOG_INFO("HorzMeshTest: fVertex test PASS");
   } else {
      LOG_INFO("HorzMeshTest: fVertex test FAIL");
      return -1;
   }

   // Test fEdge
   count = 0;
   for (int Edge = 0; Edge < LocEdges; Edge++) {
      OMEGA::R8 f = coriolis(Mesh.LatEdgeH(Edge));

      if (abs(f-Mesh.FEdgeH(Edge)) > tol) {
         count++;
      }
   }

   if ( count == 0 ) {
      LOG_INFO("HorzMeshTest: fEdge test PASS");
   } else {
      LOG_INFO("HorzMeshTest: fEdge test FAIL");
      return -1;
   }

   // Test that device arrays are identical

   // MPI_Status status;
   if (Err == 0)
      LOG_INFO("HorzMeshTest: Successful completion");
   yakl::finalize();
   MPI_Finalize();

} // end of main
//===-----------------------------------------------------------------------===/
