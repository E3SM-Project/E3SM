/*
Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
and the University Corporation for Atmospheric Research (UCAR).

Unless noted otherwise source code is licensed under the BSD license.
Additional copyright and license information can be found in the LICENSE file
distributed with this code, or at http://mpas-dev.github.com/license.html
*/

// ===================================================
//! Includes
// ===================================================

#include <algorithm>
#include <set>
#include <fstream>
#include <sstream>
#include <limits>
#include "Interface_velocity_solver.hpp"

// ===================================================
//! Namespaces
// ===================================================

#define changeTrianglesOwnership


// ice_problem pointer

int Ordering = 0; //ordering ==0 means that the mesh is extruded layerwise, whereas ordering==1 means that the mesh is extruded columnwise.
MPI_Comm comm, reducedComm;
bool isDomainEmpty = true;
bool first_time_step = true;
int nCells_F, nEdges_F, nVertices_F;
int nCellsSolve_F, nEdgesSolve_F, nVerticesSolve_F;
int nVertices, nEdges, nTriangles, globalVertexStride, globalEdgeStride,globalTriangleStride;
int maxNEdgesOnCell_F;
int const *cellsOnEdge_F, *cellsOnVertex_F, *verticesOnCell_F,
    *verticesOnEdge_F, *edgesOnCell_F, *indexToCellID_F, *indexToEdgeID_F, *indexToVertexID_F, *nEdgesOnCells_F,
    *verticesMask_F, *cellsMask_F, *dirichletCellsMask_F;
std::vector<double> layersRatio, levelsNormalizedThickness;
int nLayers;
double const *xCell_F, *yCell_F, *zCell_F, *xVertex_F,  *yVertex_F, *zVertex_F, *areaTriangle_F;
std::vector<double> xCellProjected, yCellProjected, zCellProjected;
const double unit_length = 1000;
const double T0 = 273.15;
const double secondsInAYear = 31536000.0;  // This may vary slightly in MPAS, but this should be close enough for how this is used.
double minThickness = 1e-3; //[km]
double thermal_thickness_limit; //[km]
const double minBeta = 1e-5;
double rho_ice;
double rho_ocean;
//unsigned char dynamic_ice_bit_value;
//unsigned char ice_present_bit_value;
int dynamic_ice_bit_value;
int ice_present_bit_value;

// global variables used for handling logging
char albany_log_filename[128];
int original_stdout; // the location of stdout before we captured it
int original_stderr; // the location of stderr before we captured it
int Interface_stdout; // the location of stdout as we use it here

//void *phgGrid = 0;
std::vector<int> indexToTriangleID,
    verticesOnTria, trianglesOnEdge, verticesOnEdge,
    trianglesProcIds, reduced_ranks;
std::vector<int> indexToVertexID, vertexToFCell, vertexProcIDs, triangleToFVertex, indexToEdgeID, edgeToFEdge,
    fVertexToTriangleID, fCellToVertex, iceMarginEdgesLIds, dirichletNodesIDs;
std::vector<double> dissipationHeatOnPrisms, velocityOnVertices, velocityOnCells,
    elevationData, thicknessData, betaData, bedTopographyData, stiffnessFactorData, effecPressData, muFrictionData, temperatureDataOnPrisms, smbData, thicknessOnCells, bodyForceOnBasalCell;
std::vector<bool> isVertexBoundary, isBoundaryEdge;

// only needed for creating ASCII mesh
std::vector<double> thicknessUncertaintyData;
std::vector<double> smbUncertaintyData;
std::vector<double> bmbData, bmbUncertaintyData;
std::vector<double> observedVeloXData, observedVeloYData, observedVeloUncertaintyData;
std::vector<double> observedDHDtData, observedDHDtUncertaintyData;
std::vector<double> surfaceAirTemperatureData, basalHeatFluxData;
std::vector<int> indexToCellIDData;

int numBoundaryEdges;
double radius;

exchangeList_Type const *sendCellsList_F = 0, *recvCellsList_F = 0;
exchangeList_Type const *sendEdgesList_F = 0, *recvEdgesList_F = 0;
exchangeList_Type const *sendVerticesList_F = 0, *recvVerticesList_F = 0;
exchangeList_Type sendVerticesListReversed, recvVerticesListReversed,
    sendCellsListReversed, recvCellsListReversed;

exchange::exchange(int _procID, int const* vec_first, int const* vec_last,
    int fieldDim) :
    procID(_procID), vec(vec_first, vec_last), buffer(
        fieldDim * (vec_last - vec_first)), doubleBuffer(
        fieldDim * (vec_last - vec_first)) {
}

extern "C" {

// ===================================================
//! Interface functions
// ===================================================

int velocity_solver_init_mpi(int* fComm) {
  // get MPI_Comm from Fortran
  comm = MPI_Comm_f2c(*fComm);
  reducedComm = MPI_COMM_NULL;  // initialize to null so we can check if set
  return velocity_solver_init_mpi__(comm);
}


void velocity_solver_set_parameters(double const* gravity_F, double const* ice_density_F, double const* ocean_density_F,
                         double const* sea_level_F, double const* flowParamA_F,
                         double const* flowLawExponent_F, double const* dynamic_thickness_F,
                         double const* clausius_clapeyron_coeff,
                         double const* thermal_thickness_limit_F,
                         int const* li_mask_ValueDynamicIce, int const* li_mask_ValueIce,
                         bool const* use_GLP_F) {
  // This function sets parameter values used by MPAS on the C/C++ side
  rho_ice = *ice_density_F;
  rho_ocean = *ocean_density_F;
  thermal_thickness_limit = *thermal_thickness_limit_F / unit_length; // Import with Albany scaling
  dynamic_ice_bit_value = *li_mask_ValueDynamicIce;
  ice_present_bit_value = *li_mask_ValueIce;
  velocity_solver_set_physical_parameters__(*gravity_F, rho_ice, *ocean_density_F, *sea_level_F/unit_length, *flowParamA_F*std::pow(unit_length,4)*secondsInAYear, 
                                            *flowLawExponent_F, *dynamic_thickness_F/unit_length, *use_GLP_F, *clausius_clapeyron_coeff);
}



void velocity_solver_export_2d_data(double const* lowerSurface_F,
    double const* thickness_F, double const* beta_F) {
  if (isDomainEmpty)
    return;
}

void velocity_solver_set_grid_data(int const* _nCells_F, int const* _nEdges_F,
    int const* _nVertices_F, int const* _nLevels, int const* _nCellsSolve_F,
    int const* _nEdgesSolve_F, int const* _nVerticesSolve_F,
    int const* _maxNEdgesOnCell_F, double const* radius_F,
    int const* _cellsOnEdge_F, int const* _cellsOnVertex_F,
    int const* _verticesOnCell_F, int const* _verticesOnEdge_F,
    int const* _edgesOnCell_F, int const* _nEdgesOnCells_F,
    int const* _indexToCellID_F,
    int const* _indexToEdgeID_F,
    int const* _indexToVertexID_F,
    double const* _xCell_F, double const* _yCell_F, double const* _zCell_F,
    double const* _xVertex_F, double const* _yVertex_F, double const* _zVertex_F,
    double const* _areaTriangle_F,
    int const* sendCellsArray_F, int const* recvCellsArray_F,
    int const* sendEdgesArray_F, int const* recvEdgesArray_F,
    int const* sendVerticesArray_F, int const* recvVerticesArray_F) {

  nCells_F = *_nCells_F;
  nEdges_F = *_nEdges_F;
  nVertices_F = *_nVertices_F;
  nLayers = *_nLevels-1;
  nCellsSolve_F = *_nCellsSolve_F;
  nEdgesSolve_F = *_nEdgesSolve_F;
  nVerticesSolve_F = *_nVerticesSolve_F;
  maxNEdgesOnCell_F = *_maxNEdgesOnCell_F;
  radius = *radius_F;
  cellsOnEdge_F = _cellsOnEdge_F;
  cellsOnVertex_F = _cellsOnVertex_F;
  verticesOnCell_F = _verticesOnCell_F;
  verticesOnEdge_F = _verticesOnEdge_F;
  edgesOnCell_F = _edgesOnCell_F;
  nEdgesOnCells_F = _nEdgesOnCells_F;
  indexToEdgeID_F = _indexToEdgeID_F;
  indexToCellID_F = _indexToCellID_F;
  indexToVertexID_F = _indexToVertexID_F;
  xCell_F = _xCell_F;
  yCell_F = _yCell_F;
  zCell_F = _zCell_F;
  xVertex_F = _xVertex_F;
  yVertex_F = _yVertex_F;
  zVertex_F = _zVertex_F;
  areaTriangle_F = _areaTriangle_F;

  thicknessOnCells.resize(nCellsSolve_F);

  sendCellsList_F = new exchangeList_Type(unpackMpiArray(sendCellsArray_F));
  recvCellsList_F = new exchangeList_Type(unpackMpiArray(recvCellsArray_F));
  sendEdgesList_F = new exchangeList_Type(unpackMpiArray(sendEdgesArray_F));
  recvEdgesList_F = new exchangeList_Type(unpackMpiArray(recvEdgesArray_F));
  sendVerticesList_F = new exchangeList_Type(
      unpackMpiArray(sendVerticesArray_F));
  recvVerticesList_F = new exchangeList_Type(
      unpackMpiArray(recvVerticesArray_F));

  if (radius > 10) {
    xCellProjected.resize(nCells_F);
    yCellProjected.resize(nCells_F);
    zCellProjected.assign(nCells_F, 0.);
    for (int i = 0; i < nCells_F; i++) {
      double r = std::sqrt(
          xCell_F[i] * xCell_F[i] + yCell_F[i] * yCell_F[i]
              + zCell_F[i] * zCell_F[i]);
      xCellProjected[i] = radius * std::asin(xCell_F[i] / r);
      yCellProjected[i] = radius * std::asin(yCell_F[i] / r);
    }
    xCell_F = &xCellProjected[0];
    yCell_F = &yCellProjected[0];
    zCell_F = &zCellProjected[0];
  }
}

void velocity_solver_init_l1l2(double const* levelsRatio_F) {
}




void velocity_solver_solve_l1l2(double const* lowerSurface_F,
    double const* thickness_F, double const* beta_F, double const* temperature_F,
    double* const dirichletVelocityXValue, double* const dirichletVelocitYValue,
    double* u_normal_F, double* xVelocityOnCell, double* yVelocityOnCell) {
}

void velocity_solver_export_l1l2_velocity() {}

void velocity_solver_init_fo(double const *levelsRatio_F) {

  velocityOnVertices.resize(2 * nVertices * (nLayers + 1), 0.);
  velocityOnCells.resize(2 * nCells_F * (nLayers + 1), 0.);

  if (isDomainEmpty)
    return;

  layersRatio.resize(nLayers);
  // !!Indexing of layers is reversed
  for (int i = 0; i < nLayers; i++)
    layersRatio[i] = levelsRatio_F[nLayers - 1 - i];

  mapCellsToVertices(velocityOnCells, velocityOnVertices, 2, nLayers, Ordering);
}

void velocity_solver_solve_fo(double const* bedTopography_F, double const* lowerSurface_F,
    double const* thickness_F, double * beta_F,
    double const* smb_F, double const* temperature_F, double const* stiffnessFactor_F,
    double const* effecPress_F, double const* muFriction_F,
    double* const dirichletVelocityXValue, double* const dirichletVelocitYValue,
    double* u_normal_F, double* bodyForce_F, double* dissipation_heat_F,
    double* xVelocityOnCell, double* yVelocityOnCell, double const* deltat,
    int *error) {

  std::fill(u_normal_F, u_normal_F + nEdges_F * (nLayers+1), 0.);
  //import velocity from initial guess and from dirichlet values.
  int sizeVelOnCell = nCells_F * (nLayers + 1);
  for(int iCell=0; iCell<nCells_F; ++iCell) {
    for(int il=0; il<nLayers + 1; ++il) {
      int ilReversed = nLayers - il;
      int indexReversed = iCell * (nLayers+1) + ilReversed;
      int index = iCell * (nLayers + 1) +il;
      if(dirichletCellsMask_F[indexReversed]!=0) {
        velocityOnCells[index] = dirichletVelocityXValue[indexReversed];
        velocityOnCells[index+sizeVelOnCell] = dirichletVelocitYValue[indexReversed];
      }
      else {
        velocityOnCells[index] = xVelocityOnCell[indexReversed];
        velocityOnCells[index+sizeVelOnCell] = yVelocityOnCell[indexReversed];
      }
    }
  }
  mapCellsToVertices(velocityOnCells, velocityOnVertices, 2, nLayers, Ordering);

  if (!isDomainEmpty) {

    std::vector<std::pair<int, int> > marineBdyExtensionMap;
    importFields(marineBdyExtensionMap, bedTopography_F, lowerSurface_F, thickness_F, beta_F, stiffnessFactor_F, effecPress_F, muFriction_F, temperature_F, smb_F,  minThickness);

    std::vector<double> regulThk(thicknessData);
    for (int index = 0; index < nVertices; index++)
      regulThk[index] = std::max(1e-4, thicknessData[index]);

    dissipationHeatOnPrisms.resize(nLayers * indexToTriangleID.size());
    bodyForceOnBasalCell.resize(indexToTriangleID.size());

    std::cout << "\n\nTimeStep: "<< *deltat << "\n\n"<< std::endl;

    double dt = (*deltat)/secondsInAYear;
    int albany_error;
    velocity_solver_solve_fo__(nLayers, globalVertexStride, globalTriangleStride,
        Ordering, first_time_step, indexToVertexID, indexToTriangleID, minBeta,
        regulThk, levelsNormalizedThickness, elevationData, thicknessData,
        betaData, bedTopographyData, smbData,
        stiffnessFactorData, effecPressData, muFrictionData,
        temperatureDataOnPrisms, bodyForceOnBasalCell, dissipationHeatOnPrisms, velocityOnVertices,
        albany_error, dt);
    *error=albany_error;
  }

  exportDissipationHeat(dissipation_heat_F);

  if (bodyForce_F!=nullptr) {
    exportBodyForce(bodyForce_F);
  }
  exportBeta(beta_F);

  mapVerticesToCells(velocityOnVertices, &velocityOnCells[0], 2, nLayers,
      Ordering);

  //computing x, yVelocityOnCell
  for(int iCell=0; iCell<nCells_F; ++iCell)
    for(int il=0; il<nLayers + 1; ++il) {
      int ilReversed = nLayers - il;
      int indexReversed = iCell * (nLayers+1) + ilReversed;
      int index = iCell * (nLayers + 1) +il;
    xVelocityOnCell[indexReversed] = velocityOnCells[index];
    yVelocityOnCell[indexReversed] = velocityOnCells[index+sizeVelOnCell];
  }
    get_prism_velocity_on_FEdges(u_normal_F, velocityOnCells, edgeToFEdge);
    allToAll(u_normal_F, sendEdgesList_F, recvEdgesList_F, nLayers+1);
  std::vector<double> velOnEdges(nEdges * (nLayers+1));
  for (int i = 0; i < nEdges; i++) {
    for (int il = 0; il < nLayers+1; il++) {
      velOnEdges[i * (nLayers+1) + il] = u_normal_F[edgeToFEdge[i] * (nLayers+1) + il];
    }
  }
  first_time_step = false;
}


void velocity_solver_export_fo_velocity() {

  if (isDomainEmpty)
    return;

  velocity_solver_export_fo_velocity__(reducedComm);
}

void velocity_solver_finalize() {
  velocity_solver_finalize__();
  delete sendCellsList_F;
  delete recvCellsList_F;
  delete sendEdgesList_F;
  delete recvEdgesList_F;
  delete sendVerticesList_F;
  delete recvVerticesList_F;
}

/*duality:
 *
 *   mpas(F) |  Albany LandIce (C++)
 *  ---------|---------
 *   cell    |  vertex
 *   vertex  |  triangle
 *   edge    |  edge
 *
 */

void velocity_solver_compute_2d_grid(int const* _verticesMask_F, int const* _cellsMask_F, int const* _dirichletCellsMask_F) {
  int numProcs, me;
  verticesMask_F = _verticesMask_F;
  cellsMask_F = _cellsMask_F;
  verticesMask_F = _verticesMask_F;
  dirichletCellsMask_F = _dirichletCellsMask_F;

  MPI_Comm_size(comm, &numProcs);
  MPI_Comm_rank(comm, &me);
  std::vector<int> partialOffset(numProcs + 1), globalOffsetTriangles(
      numProcs + 1), globalOffsetVertices(numProcs + 1), globalOffsetEdge(
      numProcs + 1);


  // First, we compute the FE triangles belonging to this processor.
  // If changeTrianglesOwnership is not define, the triangles belonging to this
  // processor will be the subset of the triangles (MPAS vertices) owned by this proc
  // that contain dynamic ice.
  // If changeTrianglesOwnership is define, we rearrange the ownership of the triangles
  // to improve the quality of the FE mesh and avoid corner cases (see below).

  triangleToFVertex.clear();
  triangleToFVertex.reserve(nVertices_F);
  std::vector<int> fVertexToTriangle(nVertices_F, NotAnId);

  //vector containing proc ranks for owned and shared FE triangles
  trianglesProcIds.assign(nVertices_F,NotAnId);

  //vector containing proc ranks for owned and shared MPAS cells
  std::vector<int> fCellsProcIds(nCells_F);
  getProcIds(fCellsProcIds, recvCellsList_F);
#ifdef changeTrianglesOwnership
  std::vector<int> fVerticesProcIds(nVertices_F);
  getProcIds(fVerticesProcIds, recvVerticesList_F);
  for (int i(0); i < nVertices_F; i++) {
    int cellWithMinID=nCellsSolve_F;
    if ((verticesMask_F[i] & dynamic_ice_bit_value)) {
      int minCellId = std::numeric_limits<int>::max();
      int minCellIdProc(0);

      int cellProc[3];
      bool invalidCell=false;
      for (int j = 0; j < 3; j++) {
        int iCell = cellsOnVertex_F[3 * i + j] - 1;
        if(iCell >= nCells_F) {
          invalidCell = true;
          break;
        }
        int cellID = indexToCellID_F[iCell];
        cellProc[j] = fCellsProcIds[iCell];
        if(cellID < minCellId) {
          minCellId = cellID;
          cellWithMinID = iCell;
          minCellIdProc = cellProc[j];
        }
      }

      if(invalidCell) continue;

      // the proc that owns at least 2 nodes of the triangle i. If all nodes belong to different procs, procOwns2Nodes is set to -1
      int procOwns2Nodes = ((cellProc[0] ==  cellProc[1]) || (cellProc[0] ==  cellProc[2])) ? cellProc[0] :
                           (cellProc[1] == cellProc[2]) ? cellProc[1] : -1;

      int vertexProc = fVerticesProcIds[i];
      bool triangleOwnsANode = (cellProc[0] == vertexProc) || (cellProc[1] == vertexProc) || (cellProc[2] == vertexProc);

      //A triangle will be owned by a proc if:
      // 1. the proc owns at least 2 nodes of the triangle associated to that vertex, OR
      // 2. all the nodes of the triangle belong to three different procs, and the proc owns the fortran vertex  and a node OR
      // 3. the three nodes of the triangle and the fortran vertex belong to four different procs, and the proc owns the node with the minimum ID

      trianglesProcIds[i] = (procOwns2Nodes != -1) ? procOwns2Nodes :
                       triangleOwnsANode ? vertexProc :
                       minCellIdProc;

      if (trianglesProcIds[i] == me) {
        fVertexToTriangle[i] = triangleToFVertex.size();
        triangleToFVertex.push_back(i);
      } 
    }
  }
#else
  //in this case we just set the proc ranks for owned and shared FE triangles to the be the same as MPAS owned and shared vertices
  getProcIds(trianglesProcIds, recvVerticesList_F);
  for (int i(0); i < nVerticesSolve_F; i++) {
    if (verticesMask_F[i] & dynamic_ice_bit_value) {
      fVertexToTriangle[i] = triangleToFVertex.size();
      triangleToFVertex.push_back(i);
    }
  }
#endif

  nTriangles = triangleToFVertex.size();

  //Initialize the ice sheet problem with the number of FE triangles on this prov
  initialize_iceProblem(nTriangles);

  //Create a list of global IDs for FE triangles, using MPAS vertices IDs
  fVertexToTriangleID.assign(nVertices_F, NotAnId);
  for (int index(0); index < nTriangles; index++) {
    fVertexToTriangleID[triangleToFVertex[index]] = indexToVertexID_F[triangleToFVertex[index]];
  }

#ifdef changeTrianglesOwnership
  // because we change the ownership of some triangles, we need to first communicate back to the processors that used to own those triangles
  // the data of the newly owned triangles. We do this by defining "reversed" send and receive lists, communicate back using those lists, and
  // then communicate "forward" using the usual send and receive lists.
  // We could join these two step in one communication, but for the moment we do that separately
  createReverseExchangeLists(sendVerticesListReversed, recvVerticesListReversed, trianglesProcIds, indexToVertexID_F, recvVerticesList_F);
  allToAll(fVertexToTriangleID, &sendVerticesListReversed, &recvVerticesListReversed);
  allToAll(fVertexToTriangle, &sendVerticesListReversed, &recvVerticesListReversed);
  allToAll(trianglesProcIds, sendVerticesList_F, recvVerticesList_F);
#endif
  allToAll(fVertexToTriangleID, sendVerticesList_F, recvVerticesList_F);
  allToAll(fVertexToTriangle, sendVerticesList_F, recvVerticesList_F);

  //we define the vector of global triangles Ids and compute the stride between the largest and the smallest Id globally
  //This will be needed by the velocity solver to create the 3D FE mesh.
  indexToTriangleID.resize(nTriangles);
  int maxTriangleID=std::numeric_limits<int>::min(), minTriangleID=std::numeric_limits<int>::max(), maxGlobalTriangleID, minGlobalTriangleID;
  for (int index(0); index < nTriangles; index++) {
    indexToTriangleID[index] = fVertexToTriangleID[triangleToFVertex[index]];
    maxTriangleID = (indexToTriangleID[index] > maxTriangleID) ? indexToTriangleID[index] : maxTriangleID;
    minTriangleID = (indexToTriangleID[index] < minTriangleID) ? indexToTriangleID[index] : minTriangleID;
  }

  MPI_Allreduce(&maxTriangleID, &maxGlobalTriangleID, 1, MPI_INT, MPI_MAX, comm);
  MPI_Allreduce(&minTriangleID, &minGlobalTriangleID, 1, MPI_INT, MPI_MIN, comm);
  globalTriangleStride = maxGlobalTriangleID - minGlobalTriangleID +1;

  // Second, we compute the FE edges belonging to the FE triangles owned by this processor.
  // We first compute boundary edges, and then all the other edges.
  std::vector<int> fEdgeToEdge(nEdges_F);

  std::vector<int> fEdgeToEdgeID(nEdges_F, NotAnId);
  edgeToFEdge.clear();
  isBoundaryEdge.clear();
  trianglesOnEdge.clear();

  edgeToFEdge.reserve(nEdges_F);
  trianglesOnEdge.reserve(nEdges_F * 2);
  isBoundaryEdge.reserve(nEdges_F);

  //we compute boundary edges (boundary edges must be the first edges)
  for (int i = 0; i < nEdges_F; i++) {
    ID fVertex1(verticesOnEdge_F[2 * i] - 1), fVertex2(
        verticesOnEdge_F[2 * i + 1] - 1);

    // skip the (shared) edge when the associated MPAS vertices are not  valid
    if((fVertex1>=nVertices_F) || (fVertex2>=nVertices_F))
      continue;

    ID triaId_1 = fVertexToTriangleID[fVertex1];
    ID triaId_2 = fVertexToTriangleID[fVertex2];
    bool isboundary = (triaId_1 == NotAnId) || (triaId_2 == NotAnId);

    ID iTria1 = fVertexToTriangle[fVertex1];
    ID iTria2 = fVertexToTriangle[fVertex2];
    if (trianglesProcIds[fVertex1] != me) {
      std::swap(fVertex1, fVertex2);
      std::swap(iTria1, iTria2);
    }
    bool belongsToLocalTriangle = (trianglesProcIds[fVertex1] == me);

    if (belongsToLocalTriangle) {
      if (isboundary) {
        fEdgeToEdge[i] = edgeToFEdge.size();
        edgeToFEdge.push_back(i);
        trianglesOnEdge.push_back(iTria1);
        trianglesOnEdge.push_back(iTria2);
        isBoundaryEdge.push_back(true);
      }
    }
  }

  numBoundaryEdges = edgeToFEdge.size();

  //then, we compute the other edges
  for (int i = 0; i < nEdges_F; i++) {

    ID fVertex1(verticesOnEdge_F[2 * i] - 1), fVertex2(
        verticesOnEdge_F[2 * i + 1] - 1);

    // skip the (shared) edge when the associated MPAS vertices are not  valid
    if((fVertex1>=nVertices_F) || (fVertex2>=nVertices_F))
      continue;

    ID iTria1 = fVertexToTriangle[fVertex1];
    ID iTria2 = fVertexToTriangle[fVertex2];

    ID triaId_1 = fVertexToTriangleID[fVertex1]; //global Triangle
    ID triaId_2 = fVertexToTriangleID[fVertex2]; //global Triangle

    if (trianglesProcIds[fVertex1] != me) {
      std::swap(fVertex1, fVertex2);
      std::swap(iTria1, iTria2);
    }

    bool belongsToAnyTriangle = (triaId_1 != NotAnId) || (triaId_2 != NotAnId);
    bool isboundary = (triaId_1 == NotAnId) || (triaId_2 == NotAnId);
    bool belongsToLocalTriangle = (trianglesProcIds[fVertex1] == me);
    bool isMine = i < nEdgesSolve_F;

    if (belongsToLocalTriangle && !isboundary) {
      fEdgeToEdge[i] = edgeToFEdge.size();
      edgeToFEdge.push_back(i);
      trianglesOnEdge.push_back(iTria1);
      trianglesOnEdge.push_back(iTria2);
      isBoundaryEdge.push_back(false);
    }
  }

  for (int fEdge = 0; fEdge < nEdges_F; fEdge++)
    fEdgeToEdgeID[fEdge] = indexToEdgeID_F[fEdge];

  nEdges = edgeToFEdge.size();
  indexToEdgeID.resize(nEdges);
  iceMarginEdgesLIds.clear();
  iceMarginEdgesLIds.reserve(numBoundaryEdges);
  int maxEdgeID=std::numeric_limits<int>::min(), minEdgeID=std::numeric_limits<int>::max(), maxGlobalEdgeID, minGlobalEdgeID;
  for (int index = 0; index < nEdges; index++) {
    int fEdge = edgeToFEdge[index];
    indexToEdgeID[index] = fEdgeToEdgeID[fEdge];
    maxEdgeID = (indexToEdgeID[index] > maxEdgeID) ? indexToEdgeID[index] : maxEdgeID;
    minEdgeID = (indexToEdgeID[index] < minEdgeID) ? indexToEdgeID[index] : minEdgeID;

    if(index<numBoundaryEdges){
      int fCell0 = cellsOnEdge_F[2 * fEdge] - 1;
      int fCell1 = cellsOnEdge_F[2 * fEdge + 1] - 1;
      bool isCell0OnMargin = !(cellsMask_F[fCell0] & dynamic_ice_bit_value) &&
          (dirichletCellsMask_F[(nLayers+1)*fCell0] == 0);
      bool isCell1OnMargin = !(cellsMask_F[fCell1] & dynamic_ice_bit_value) &&
          (dirichletCellsMask_F[(nLayers+1)*fCell1] == 0);
      if(isCell0OnMargin || isCell1OnMargin)
        iceMarginEdgesLIds.push_back(index);
    }
  }

  MPI_Allreduce(&maxEdgeID, &maxGlobalEdgeID, 1, MPI_INT, MPI_MAX, comm);
  MPI_Allreduce(&minEdgeID, &minGlobalEdgeID, 1, MPI_INT, MPI_MIN, comm);
  globalEdgeStride = maxGlobalEdgeID - minGlobalEdgeID + 1;

  // Third, we compute the FE vertices belonging to the FE triangles owned by this processor.
  // We need to make sure that an FE vertex is owned by a proc that owns a FE triangle that contain that vertex
  // Otherwise we might end up with weird situation where a vertex could belong to a process with no associated triangle.

  vertexToFCell.clear();
  vertexToFCell.reserve(nCells_F);

  fCellToVertex.assign(nCells_F, NotAnId);
  std::vector<int> fCellToVertexID(nCells_F, NotAnId);

  vertexProcIDs.clear();

  std::vector<int> verticesProcIds(nCells_F, NotAnId);

  vertexProcIDs.reserve(nCells_F);
  for (int i = 0; i < nCells_F; i++) {
    bool isMine = i < nCellsSolve_F;
    bool belongsToLocalTriangle = false;
    bool belongsToAnyTriangle = false;
    int nEdg = nEdgesOnCells_F[i];
    bool invalidVertex = false;
    int minTriangleProcId = std::numeric_limits<int>::max();
    bool nodeOwnedByATriaProc = false;
    for (int j = 0; j < nEdg; j++) {
      ID fVertex(verticesOnCell_F[maxNEdgesOnCell_F * i + j] - 1);
      if(fVertex >= nVertices_F) {
        invalidVertex = true;
        break;
      }

      ID triaId = fVertexToTriangleID[fVertex];
      belongsToLocalTriangle = belongsToLocalTriangle || (trianglesProcIds[fVertex] == me);
      belongsToAnyTriangle = belongsToAnyTriangle || (triaId != NotAnId);

      if(triaId != NotAnId) {
        nodeOwnedByATriaProc = nodeOwnedByATriaProc || (trianglesProcIds[fVertex] == fCellsProcIds[i]);
        minTriangleProcId = (trianglesProcIds[fVertex] < minTriangleProcId) ? trianglesProcIds[fVertex] : minTriangleProcId;
      }
    }

    if(invalidVertex) continue;

    if(belongsToAnyTriangle)
      verticesProcIds[i] = nodeOwnedByATriaProc ? fCellsProcIds[i] : minTriangleProcId;

    if (belongsToLocalTriangle) {
      fCellToVertex[i] = vertexToFCell.size();
      vertexToFCell.push_back(i);
      vertexProcIDs.push_back(reduced_ranks[verticesProcIds[i]]);
    }
  }
  nVertices = vertexToFCell.size();

  for (int fcell = 0; fcell < nCells_F; fcell++)
    fCellToVertexID[fcell] = indexToCellID_F[fcell];

  int maxVertexID=std::numeric_limits<int>::min(), minVertexID=std::numeric_limits<int>::max(), maxGlobalVertexID, minGlobalVertexID;
  indexToVertexID.resize(nVertices);
  for (int index = 0; index < nVertices; index++) {
    int fCell = vertexToFCell[index];
    indexToVertexID[index] = fCellToVertexID[fCell];
    maxVertexID = (indexToVertexID[index] > maxVertexID) ? indexToVertexID[index] : maxVertexID;
    minVertexID = (indexToVertexID[index] < minVertexID) ? indexToVertexID[index] : minVertexID;
  }

  MPI_Allreduce(&maxVertexID, &maxGlobalVertexID, 1, MPI_INT, MPI_MAX, comm);
  MPI_Allreduce(&minVertexID, &minGlobalVertexID, 1, MPI_INT, MPI_MIN, comm);
  globalVertexStride = maxGlobalVertexID - minGlobalVertexID + 1;


  int vertexColumnShift = (Ordering == 1) ? 1 : globalVertexStride;
  int vertexLayerShift = (Ordering == 0) ? 1 : nLayers + 1;
  dirichletNodesIDs.clear();
  dirichletNodesIDs.reserve(nVertices); //need to improve storage efficiency
  isVertexBoundary.assign(nVertices, false);
  for (int index = 0; index < nVertices; index++) {
    int fCell = vertexToFCell[index];
    for(int il=0; il< nLayers+1; ++il)
    {
      int imask_F = il+(nLayers+1)*fCell;
      if(dirichletCellsMask_F[imask_F]!=0)
        dirichletNodesIDs.push_back((nLayers-il)*vertexColumnShift+indexToVertexID[index]*vertexLayerShift);
    }

    int nEdg = nEdgesOnCells_F[fCell];
    int j = 0;
    bool isBoundary;
    do {
      int fVertex = verticesOnCell_F[maxNEdgesOnCell_F * fCell + j++] - 1;
      isBoundary = !(verticesMask_F[fVertex] & dynamic_ice_bit_value);
    } while ((j < nEdg) && (!isBoundary));
    isVertexBoundary[index] = isBoundary;
  }

  // because we change the ownership of some vertices, we need to first communicate back to the processors that used to own those vertices
  // the data of the newly owned vertices. We do this by defining "reversed" send and receive lists, communicate back using that list, and
  // then communicate forward with the usual send and receive lists.
  // We could join these two step in one communication, but for the moment we do that separately.
  // We need to communicate info about the vertices when we get the ice velocity on vertices form the velocity solver/

  createReverseExchangeLists(sendCellsListReversed, recvCellsListReversed, verticesProcIds, indexToCellID_F, recvCellsList_F);

  //construct the local vector vertices on triangles making sure the area is positive
  verticesOnTria.resize(nTriangles * 3);
  double x[3], y[3], z[3];
  for (int index = 0; index < nTriangles; index++) {
    int iTria = triangleToFVertex[index];

    for (int j = 0; j < 3; j++) {
      int iCell = cellsOnVertex_F[3 * iTria + j] - 1;
      verticesOnTria[3 * index + j] = fCellToVertex[iCell];
      x[j] = xCell_F[iCell];
      y[j] = yCell_F[iCell];
      //  z[j] = zCell_F[iCell];
    }
    if (signedTriangleArea(x, y) < 0)
      std::swap(verticesOnTria[3 * index + 1], verticesOnTria[3 * index + 2]);
  }

  verticesOnEdge.resize(2 * nEdges);
  for (int index = 0; index < nEdges; index++) {
    int fEdge = edgeToFEdge[index];
    int fCell1 = cellsOnEdge_F[2 * fEdge] - 1;
    int fCell2 = cellsOnEdge_F[2 * fEdge + 1] - 1;
    verticesOnEdge[2 * index] = fCellToVertex[fCell1];
    verticesOnEdge[2 * index + 1] = fCellToVertex[fCell2];
  }


  //call the velocity solver only on procs owning some FE triangles
  if (isDomainEmpty)
    return;

  velocity_solver_compute_2d_grid__(reducedComm);
}

void velocity_solver_extrude_3d_grid(double const* levelsRatio_F) {

  if (isDomainEmpty)
    return;

  layersRatio.resize(nLayers);
  // !!Indexing of layers is reversed
  for (int i = 0; i < nLayers; i++)
    layersRatio[i] = levelsRatio_F[nLayers - 1 - i];

  levelsNormalizedThickness.resize(nLayers + 1);

  levelsNormalizedThickness[0] = 0;
  for (int i = 0; i < nLayers; i++)
    levelsNormalizedThickness[i + 1] = levelsNormalizedThickness[i]
        + layersRatio[i];

  //construct the local vector of coordinates
  std::vector<double> verticesCoords(3 * nVertices);

  for (int index = 0; index < nVertices; index++) {
    int iCell = vertexToFCell[index];
    verticesCoords[index * 3] = xCell_F[iCell] / unit_length;
    verticesCoords[index * 3 + 1] = yCell_F[iCell] / unit_length;
    verticesCoords[index * 3 + 2] = zCell_F[iCell] / unit_length;
  }

  std::vector<std::vector<int>> procsSharingVertices(nVertices);
  for(int i=0; i< nVertices; i++)
    procsSharingVertex(i, procsSharingVertices[i]);


  std::vector<int> iceMarginEdgesGIds(iceMarginEdgesLIds.size());
  for (int i=0; i< iceMarginEdgesGIds.size(); ++i) {
    iceMarginEdgesGIds[i] = indexToEdgeID[iceMarginEdgesLIds[i]];
  }

  velocity_solver_extrude_3d_grid__(nLayers, globalTriangleStride, globalVertexStride,
      globalEdgeStride, Ordering, reducedComm, indexToVertexID, vertexProcIDs,
      verticesCoords, verticesOnTria, procsSharingVertices, isBoundaryEdge,
      trianglesOnEdge, verticesOnEdge, indexToEdgeID,
      indexToTriangleID, dirichletNodesIDs, iceMarginEdgesGIds);
  }

// Function to set up how the MPAS log file will be used by Albany
void interface_init_log(){
  int me, tasks;
  MPI_Comm_rank(comm, &me);
  MPI_Comm_size(comm, &tasks);

  char fname[128];
  char oformat[128];

  strcpy(oformat, "log.albany.");

  if(tasks < 1E4) {
     strcat(oformat, "%4.4i");
  }
  else if (tasks < 1E5) {
     strcat(oformat, "%5.5i");
  }
  else if (tasks < 1E6) {
     strcat(oformat, "%6.6i");
  }
  else if (tasks < 1E7) {
     strcat(oformat, "%7.7i");;
  }
  else if (tasks < 1E8) {
     strcat(oformat, "%8.8i");
  }
  else if (tasks < 1E9) {
     strcat(oformat, "%9.9i");
  }
  else {
     if(me == 0){
       fprintf(stderr, "Error opening Albany stdout for 1E9 tasks or more\n");
     }
     return;
  }
  strcat(oformat, ".out");

  sprintf(albany_log_filename, oformat, me);

  Interface_stdout = open(albany_log_filename, O_CREAT|O_WRONLY|O_TRUNC,0644);
  if(Interface_stdout >=  0) {
     write(Interface_stdout, "-- Beginning log file for output from Albany velocity solver --", 63);
     write(Interface_stdout, " ", 1);
     fsync(Interface_stdout);
  } else {
    std::cerr << "Error opening Albany stdout file." << std::endl;
  }

}

// Function to redirect Albany stdout to the MPAS log file
void interface_redirect_stdout(int const* iTimestep) {
  /* Save current stdout for use later */
  //fsync(Interface_stdout);
  fflush(stdout);
  fflush(stderr);
  original_stdout = dup(1);
  original_stderr = dup(2);
  dup2(Interface_stdout, 1);
  dup2(Interface_stdout, 2);
  if (*iTimestep >= 0) {
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "--- Beginning Albany velocity solve for timestep " << *iTimestep << " ---" << std::endl;
    std::cout << std::endl;
  }
}


// Function to return to stdout to its previous location
void interface_reset_stdout() {
  /* Restore stdout */
     //fsync(Interface_stdout);
  fflush(stdout);
  fflush(stderr);
  dup2(original_stdout, 1);
  dup2(original_stderr, 2);
  close(original_stdout);
  close(original_stderr);
}


}

//This function computes the normal velocity on the edges midpoints of MPAS 2d cells.
//The function rely on the assumption that the locations of the edges midpoints are inside the union of the triangles that share the dual of that edge.
//The normal velocity is computed assuming to have linear/bilinear finite elements on the prisms.
//The normal velocity will be second order accurate for other finite elements (e.g. linear or quadratic on Tetrahedra).

void get_prism_velocity_on_FEdges(double * uNormal,
    const std::vector<double>& velocityOnCells,
    const std::vector<int>& edgeToFEdge) {

  //using layout of velocityOnCells
  int columnShift = 1;
  int layerShift = (nLayers + 1);

  UInt nPoints3D = nCells_F * (nLayers + 1);

  for (int iEdge = 0; iEdge < nEdgesSolve_F; iEdge++) {
    int iCell0 = cellsOnEdge_F[2 * iEdge] - 1;
    int iCell1 = cellsOnEdge_F[2 * iEdge + 1] - 1;

    //computing normal to the cell edge (dual of triangular edge)
    double nx = xCell_F[iCell1] - xCell_F[iCell0];
    double ny = yCell_F[iCell1] - yCell_F[iCell0];
    double n = sqrt(nx * nx + ny * ny);
    nx /= n;
    ny /= n;

    //identifying triangles that shares the edge
    ID fVertex0 = verticesOnEdge_F[2 * iEdge] - 1;
    ID fVertex1 = verticesOnEdge_F[2 * iEdge + 1] - 1;


    int iTria0 = fVertexToTriangleID[fVertex0];
    int iTria1 = fVertexToTriangleID[fVertex1];
    if((iTria0 == NotAnId) && (iTria1 == NotAnId)) continue;

    double t0[2*3], t1[2*3]; //t0[0] contains the x-coords of vertices of triangle 0 and t0[1] its y-coords.
    for (int j = 0; j < 3; j++) {
      int iCell = cellsOnVertex_F[3 * fVertex0 + j] - 1;
      t0[0 + 2 * j] = xCell_F[iCell];
      t0[1 + 2 * j] = yCell_F[iCell];
      iCell = cellsOnVertex_F[3 * fVertex1 + j] - 1;
      t1[0 + 2 * j] = xCell_F[iCell];
      t1[1 + 2 * j] = yCell_F[iCell];
    }

    //getting triangle circumcenters (vertices of MPAS cells).
    double bcoords[3];

   ID iCells[3]; //iCells[k] is the array of cells indexes of triangle k on iEdge

   //computing midpoint of fortran edge
   double e_mid[2];
   e_mid[0] =  0.5*(xVertex_F[fVertex0] + xVertex_F[fVertex1]);
   e_mid[1] =  0.5*(yVertex_F[fVertex0] + yVertex_F[fVertex1]);

   if((verticesMask_F[fVertex0] & dynamic_ice_bit_value) && belongToTria(e_mid, t0, bcoords)) {
      // triangle1 is in the mesh  AND midpoint is in triangle1
      for (int j = 0; j < 3; j++)
        iCells[j] = cellsOnVertex_F[3 * fVertex0 + j] - 1;
      }
   else if((verticesMask_F[fVertex1] & dynamic_ice_bit_value) && belongToTria(e_mid, t1, bcoords)) {
      //triangle2 is in the mesh  AND midpoint is in triangle2
      for (int j = 0; j < 3; j++)
        iCells[j] = cellsOnVertex_F[3 * fVertex1 + j] - 1;
      }
   else if((iTria0 == NotAnId) || (iTria1 == NotAnId)) { //is on Boundary
      //edge is on boundary and wasn't found by previous two cases
      //For boundary edges one of the two triangles sharing the edge won't be part of the velocity solver's mesh and the dynamic_ice_bit_value will be 0.

      //Compute iCells containing the vertices of the triangle that is part of the mesh and bcoords the corresponding barycentric coordinates
      if(verticesMask_F[fVertex0] & dynamic_ice_bit_value) { belongToTria(e_mid, t0, bcoords);
        for (int j = 0; j < 3; j++)
          iCells[j] = cellsOnVertex_F[3 * fVertex0 + j] - 1;
        }
      else {
        belongToTria(e_mid, t1, bcoords);
        for (int j = 0; j < 3; j++)
          iCells[j] = cellsOnVertex_F[3 * fVertex1 + j] - 1;
        }

      //We modify the barycentric coordinates so that they will be all non-negative (corresponding to a point on the triangle edge).
      double sum(0);
      for(int j=0; j<3; j++) {
        bcoords[j] = std::max(0.,bcoords[j]);
        sum += bcoords[j];
        }
      //Scale the coordinates so that they sum to 1.
      for(int j=0; j<3; j++)
        bcoords[j] /= sum;
      }
   else { //error, edge midpont does not belong to either triangle
      std::cout << "Error, edge midpont does not belong to either triangle" << std::endl;

      for (int j = 0; j < 3; j++)
        std::cout << "("<<t0[0 + 2 * j]<<","<<t0[1 + 2 * j]<<") ";
      std::cout <<std::endl;
      for (int j = 0; j < 3; j++)
        std::cout << "("<<t1[0 + 2 * j]<<","<<t1[1 + 2 * j]<<") ";
      std::cout <<"\n midpoint: ("<<e_mid[0]<<","<<e_mid[1]<<")"<<std::endl;
      exit(1);
   }

    for (int il = 0; il < nLayers+1; il++) { //loop over layers
      int ilReversed = nLayers - il;
      int index = iEdge * (nLayers+1) + ilReversed;
      uNormal[index] = 0;
      for(int j=0; j<3; ++j) {
        ID iCell3D = il * columnShift + layerShift * iCells[j];
        //contribution of x component of velocity at the edge midpoint
        uNormal[index] += bcoords[j] * nx * velocityOnCells[iCell3D];
        //getting IDs for y component of velocities
        iCell3D += nPoints3D;
        //contribution of y component of velocity at the edge midpoint
        uNormal[index] += bcoords[j] * ny * velocityOnCells[iCell3D];
      }
    } //end loop over layers
  } //loop over edges
}

void mapVerticesToCells(const std::vector<double>& velocityOnVertices,
    double* velocityOnCells, int fieldDim, int numLayers, int ordering) {
  int lVertexColumnShift = (ordering == 1) ? 1 : nVertices;
  int vertexLayerShift = (ordering == 0) ? 1 : numLayers + 1;

  int nVertices3D = nVertices * (numLayers + 1);

  // 0 entire field so no values from previous solve are left behind
  // if the ice extent has retreated.
  std::fill(velocityOnCells, velocityOnCells + nCells_F * (numLayers+1) * fieldDim, 0.);

  for (UInt j = 0; j < nVertices3D; ++j) {
    int ib = (ordering == 0) * (j % lVertexColumnShift)
        + (ordering == 1) * (j / vertexLayerShift);
    int il = (ordering == 0) * (j / lVertexColumnShift)
        + (ordering == 1) * (j % vertexLayerShift);

    int iCell = vertexToFCell[ib];
    int cellIndex = iCell * (numLayers + 1) + il;
    int vertexIndex = j;
    for (int dim = 0; dim < fieldDim; dim++) {
      velocityOnCells[cellIndex] = velocityOnVertices[vertexIndex];
      cellIndex += nCells_F * (numLayers + 1);
      vertexIndex += nVertices3D;
    }
  }

  for (int dim = 0; dim < fieldDim; dim++) {
    allToAll(&velocityOnCells[dim * nCells_F * (numLayers + 1)],
        &sendCellsListReversed, &recvCellsListReversed, (numLayers + 1));
    allToAll(&velocityOnCells[dim * nCells_F * (numLayers + 1)],
        sendCellsList_F, recvCellsList_F, (numLayers + 1));
  }
}


void createReverseExchangeLists(exchangeList_Type& sendListReverse_F,
    exchangeList_Type& receiveListReverse_F,
    const std::vector<int>& newProcIds, const int* indexToID_F, exchangeList_Type const * recvList_F) {
  sendListReverse_F.clear();
  receiveListReverse_F.clear();
  std::map<int, std::map<int, int> > sendMap, receiveMap;
  int nFEntities = newProcIds.size();
  std::vector<int> procIds(nFEntities);
  getProcIds(procIds, recvList_F);
  int me;
  MPI_Comm_rank(comm, &me);
  for (int fEntity = 0; fEntity < nFEntities; fEntity++) {
    if ((procIds[fEntity] != me) && (newProcIds[fEntity] == me)) {
      sendMap[procIds[fEntity]].insert(
        std::make_pair(indexToID_F[fEntity], fEntity));
    }
  }

  for (std::map<int, std::map<int, int> >::const_iterator it = sendMap.begin();
      it != sendMap.end(); it++) {
    std::vector<int> sendVec(it->second.size());
    int i = 0;
    for (std::map<int, int>::const_iterator iter = it->second.begin();
        iter != it->second.end(); iter++)
      sendVec[i++] = iter->second;
    sendListReverse_F.push_back(
        exchange(it->first, &sendVec[0], &sendVec[0] + sendVec.size()));
  }

  for (int fEntity = 0; fEntity < nFEntities; fEntity++) {
    if((procIds[fEntity] == me) && (newProcIds[fEntity] != NotAnId) && (newProcIds[fEntity] != me)) {
      receiveMap[newProcIds[fEntity]].insert(
        std::make_pair(indexToID_F[fEntity], fEntity));
    }
  }

  for (std::map<int, std::map<int, int> >::const_iterator it =
      receiveMap.begin(); it != receiveMap.end(); it++) {
    std::vector<int> receiveVec(it->second.size());
    int i = 0;
    for (std::map<int, int>::const_iterator iter = it->second.begin();
        iter != it->second.end(); iter++)
      receiveVec[i++] = iter->second;
    receiveListReverse_F.push_back(
        exchange(it->first, &receiveVec[0],
            &receiveVec[0] + receiveVec.size()));
  }
}

void mapCellsToVertices(const std::vector<double>& velocityOnCells,
    std::vector<double>& velocityOnVertices, int fieldDim, int numLayers,
    int ordering) {
  int lVertexColumnShift = (ordering == 1) ? 1 : nVertices;
  int vertexLayerShift = (ordering == 0) ? 1 : numLayers + 1;

  int nVertices3D = nVertices * (numLayers + 1);
  for (UInt j = 0; j < nVertices3D; ++j) {
    int ib = (ordering == 0) * (j % lVertexColumnShift)
        + (ordering == 1) * (j / vertexLayerShift);
    int il = (ordering == 0) * (j / lVertexColumnShift)
        + (ordering == 1) * (j % vertexLayerShift);

    int iCell = vertexToFCell[ib];
    int cellIndex = iCell * (numLayers + 1) + il;
    int vertexIndex = j;
    for (int dim = 0; dim < fieldDim; dim++) {
      velocityOnVertices[vertexIndex] = velocityOnCells[cellIndex];
      cellIndex += nCells_F * (numLayers + 1);
      vertexIndex += nVertices3D;
    }
  }
}

double signedTriangleArea(const double* x, const double* y) {
  double u[2] = { x[1] - x[0], y[1] - y[0] };
  double v[2] = { x[2] - x[0], y[2] - y[0] };

  return 0.5 * (u[0] * v[1] - u[1] * v[0]);
}

double signedTriangleAreaOnSphere(const double* x, const double* y,
    const double *z) {
  double u[3] = { x[1] - x[0], y[1] - y[0], z[1] - z[0] };
  double v[3] = { x[2] - x[0], y[2] - y[0], z[2] - z[0] };

  double crossProduct[3] = { u[1] * v[2] - u[2] * v[1], u[2] * v[0]
      - u[0] * v[2], u[0] * v[1] - u[1] * v[0] };
  double area = 0.5
      * std::sqrt(
          crossProduct[0] * crossProduct[0] + crossProduct[1] * crossProduct[1]
              + crossProduct[2] * crossProduct[2]);
  return
      (crossProduct[0] * x[0] + crossProduct[1] * y[0] + crossProduct[2] * z[0]
          > 0) ? area : -area;
}


void importFields(std::vector<std::pair<int, int> >& marineBdyExtensionMap,  double const* bedTopography_F, double const * lowerSurface_F, double const * thickness_F,
    double const * beta_F, double const* stiffnessFactor_F, double const* effecPress_F, double const* muFriction_F,
    double const * temperature_F, double const * smb_F, double eps) {

  int vertexLayerShift = (Ordering == 0) ? 1 : nLayers + 1;
  elevationData.assign(nVertices, 1e10);
  thicknessData.assign(nVertices, 1e10);
  bedTopographyData.assign(nVertices, 1e10);
  if (beta_F != 0)
    betaData.assign(nVertices, 1e10);
  if (stiffnessFactor_F != 0)
    stiffnessFactorData.assign(nVertices, 1.0);
  if (effecPress_F != 0)
    effecPressData.assign(nVertices, 1e10);
  if (muFriction_F!= 0)
    muFrictionData.assign(nVertices, 1e10);
  if(temperature_F != 0)
    temperatureDataOnPrisms.assign(nLayers * nTriangles, 1e10);
  if (smb_F != 0)
    smbData.assign(nVertices, 1e10);


  //import fields
  for (int index = 0; index < nVertices; index++) {
    int iCell = vertexToFCell[index];
    thicknessData[index] = std::max(thickness_F[iCell] / unit_length, eps);
    bedTopographyData[index] = bedTopography_F[iCell] / unit_length;
    elevationData[index] = lowerSurface_F[iCell] / unit_length + thicknessData[index];
    if (beta_F != 0)
      betaData[index] = beta_F[iCell] / unit_length;
    if (smb_F != 0)
      smbData[index] = smb_F[iCell] * secondsInAYear/rho_ice;
    if (stiffnessFactor_F != 0)
      stiffnessFactorData[index] = stiffnessFactor_F[iCell];
    if (effecPress_F != 0)
      effecPressData[index] = effecPress_F[iCell] / unit_length;  
    if (muFriction_F != 0)
      muFrictionData[index] = muFriction_F[iCell];
  }

  int lElemColumnShift = (Ordering == 1) ? 1 : nTriangles;
  int elemLayerShift = (Ordering == 0) ? 1 : nLayers;
  if(temperature_F != 0) {
    for (int index = 0; index < nTriangles; index++) {
      for (int il = 0; il < nLayers; il++) {
        double temperature = 0;
        int ilReversed = nLayers - il - 1;
        int nPoints = 0;
        for (int iVertex = 0; iVertex < 3; iVertex++) {
      int v = verticesOnTria[iVertex + 3 * index];
      int iCell = vertexToFCell[v];
      //compute temperature by averaging temperature values of triangles vertices where ice is present
      // Note that thermal_thickness_limit was imported in Albany units (km) but thickness_F is still in 
      // MPAS units (m), so thermal_thickness_limit is being scaled back to MPAS units for this comparison.
      if (thickness_F[iCell] > thermal_thickness_limit * unit_length) {
        temperature += temperature_F[iCell * nLayers + ilReversed];
        nPoints++;
         }
        }
        if (nPoints == 0)  //if triangle is in an ice-free area, set the temperature to T0
      temperatureDataOnPrisms[index*elemLayerShift + il*lElemColumnShift] = T0;
        else
      temperatureDataOnPrisms[index*elemLayerShift + il*lElemColumnShift] = temperature / nPoints;
      }
    }
  }

  //extend thickness and elevation data to the border for marine vertices
  marineBdyExtensionMap.clear();
  marineBdyExtensionMap.reserve(nVertices);

  for (int iV = 0; iV < nVertices; iV++) {
    int fCell = vertexToFCell[iV];
    bool isDynamicIceVertex = cellsMask_F[fCell] & dynamic_ice_bit_value;

    // Loop over boundary vertices to correctly extend relevant fields
    // Note, isDynamicIceVertex is typically true in the interior of the FE mesh and false at the boundary
    // However, isDynamicIceVertex can be false at interior points. This happens if there is narrow (one MPAS cell width) tongue
    // without ice surrounded by ice. After the extension the tongue will be part of the FE mesh,
    // and we need to properly extend the fields there as well.
    if (!isDynamicIceVertex) {
      if(bedTopographyData[iV]<0) {
        // -- marine margin --
        // Identify the lowest elevation neighboring cell with ice
        // Scalar values will be mapped from that location to here.
        double elevTemp =1e10;
        int nEdg = nEdgesOnCells_F[fCell];
        int neighbor_cell = -1;
        for (int j = 0; j < nEdg; j++) {
          int fEdge = edgesOnCell_F[maxNEdgesOnCell_F * fCell + j] - 1;
          int c0 = cellsOnEdge_F[2 * fEdge] - 1;
          int c1 = cellsOnEdge_F[2 * fEdge + 1] - 1;
          int c = (fCellToVertex[c0] == iV) ? c1 : c0;
          if((cellsMask_F[c] & dynamic_ice_bit_value)) {
            double elev = thickness_F[c] + lowerSurface_F[c]; // - 1e-8*std::sqrt(pow(xCell_F[c0],2)+std::pow(yCell_F[c0],2));
            if (elev < elevTemp) {
              elevTemp = elev;
              neighbor_cell = c;
            }
          }
        }
        // check if we didn't assign anything.  This occurs if this node has no neighbors with (dynamic) ice.
        // This should not ever occur, but including the check just to be safe.
        if (neighbor_cell != -1) {
          marineBdyExtensionMap.push_back(std::make_pair(iV,neighbor_cell));
        } else {
          std::cout << "WARNING: vertex with ID " << indexToVertexID[iV] <<
              " is not connected to active ice." << std::endl;
          // Check that this margin location is not below sea level!
          thicknessData[iV] = eps*3.0; // insert special value here to make identifying these points easier in exo output
          elevationData[iV] = (1.0 - rho_ice / rho_ocean) * thicknessData[iV];  // floating surface
        }
      } else {
        // -- nonmarine margin -- (extend to zero thickness)

        thicknessData[iV] = eps;
        elevationData[iV] = bedTopographyData[iV]+eps;
      }
    } // is on margin
  }  // vertex loop

  // Apply extension on marine margin
  for (std::vector<std::pair<int, int> >::iterator it = marineBdyExtensionMap.begin();
      it != marineBdyExtensionMap.end(); ++it) {
    int iv = it->first;
    int ic = it->second;

    double bed = bedTopographyData[iv];
    double elev = (lowerSurface_F[ic]+thickness_F[ic]) / unit_length;
    double thick = thickness_F[ic] / unit_length;

    //assume elevation as given and adjust thickness to avoid unphysical situations
    thick = std::min(thick, elev - bed);
    thick = std::min(thick, rho_ocean/(rho_ocean-rho_ice)*elev);

    if(thick < eps) { //thickness needs to be greater than eps
      thicknessData[iv] = eps;
      elevationData[iv] = std::max(bed+eps, (1.0 - rho_ice/rho_ocean)*eps);
    } else {
      thicknessData[iv] = thick;
      elevationData[iv] = elev;
    }

  }
}

void import2DFieldsObservations(std::vector<std::pair<int, int> >& marineBdyExtensionMap,
            double const * thicknessUncertainty_F,
            double const * smbUncertainty_F,
            double const * bmb_F, double const * bmbUncertainty_F,
            double const * observedSurfaceVelocityX_F, double const * observedSurfaceVelocityY_F,
            double const * observedSurfaceVelocityUncertainty_F,
            double const * observedThicknessTendency_F, double const * observedThicknessTendencyUncertainty_F,
            // these last 2 thermal fields are not necessarily observations but
            // they are not needed except for exporting to ascii format, so including
            // them in this function.
            double const* surfaceAirTemperature_F, double const* basalHeatFlux_F,
            int const * indexToCellID_F) {

  thicknessUncertaintyData.assign(nVertices, 1e10);
  smbUncertaintyData.assign(nVertices, 1e10);
  bmbData.assign(nVertices, 1e10);
  bmbUncertaintyData.assign(nVertices, 1e10);
  observedVeloXData.assign(nVertices, 1e10);
  observedVeloYData.assign(nVertices, 1e10);
  observedVeloUncertaintyData.assign(nVertices, 1e10);
  observedDHDtData.assign(nVertices, 1e10);
  observedDHDtUncertaintyData.assign(nVertices, 1e10);
  surfaceAirTemperatureData.assign(nVertices, 1e10);
  basalHeatFluxData.assign(nVertices, 1e10);
  indexToCellIDData.assign(nVertices, -1);


  //import fields
  for (int index = 0; index < nVertices; index++) {
    int iCell = vertexToFCell[index];

    thicknessUncertaintyData[index] = thicknessUncertainty_F[iCell] / unit_length;
    smbUncertaintyData[index] = smbUncertainty_F[iCell] * secondsInAYear / rho_ice;
    bmbData[index] = bmb_F[iCell] * secondsInAYear / rho_ice;
    bmbUncertaintyData[index] = bmbUncertainty_F[iCell] * secondsInAYear / rho_ice;

    observedVeloXData[index] = observedSurfaceVelocityX_F[iCell] * secondsInAYear;
    observedVeloYData[index] = observedSurfaceVelocityY_F[iCell] * secondsInAYear;
    observedVeloUncertaintyData[index] = observedSurfaceVelocityUncertainty_F[iCell] * secondsInAYear;

    observedDHDtData[index] = observedThicknessTendency_F[iCell] * secondsInAYear;
    observedDHDtUncertaintyData[index] = observedThicknessTendencyUncertainty_F[iCell] * secondsInAYear;

    surfaceAirTemperatureData[index] = surfaceAirTemperature_F[iCell];
    basalHeatFluxData[index] = basalHeatFlux_F[iCell];

    indexToCellIDData[index] = indexToCellID_F[iCell];
  }
}

void exportDissipationHeat(double * dissipationHeat_F) {
  std::fill(dissipationHeat_F, dissipationHeat_F + nVertices_F * (nLayers), 0.);
  int lElemColumnShift = (Ordering == 1) ? 1 : indexToTriangleID.size();
  int elemLayerShift = (Ordering == 0) ? 1 : nLayers;
  for (int index = 0; index < nTriangles; index++) {
    for (int il = 0; il < nLayers; il++) {
      int ilReversed = nLayers - il - 1;
      int fVertex = triangleToFVertex[index];
      dissipationHeat_F[fVertex * nLayers + ilReversed] = dissipationHeatOnPrisms[index * elemLayerShift + il * lElemColumnShift];
    }
  }
#ifdef  changeTrianglesOwnership
  allToAll (dissipationHeat_F,  &sendVerticesListReversed, &recvVerticesListReversed, nLayers);
#endif
  allToAll (dissipationHeat_F,  sendVerticesList_F, recvVerticesList_F, nLayers);
}

void exportBodyForce(double * bodyForce_F) {
  std::fill(bodyForce_F, bodyForce_F + nVertices_F, 0.);
  for (int index = 0; index < nTriangles; index++) {
    int fVertex = triangleToFVertex[index];
    bodyForce_F[fVertex] = bodyForceOnBasalCell[index];
  }
#ifdef  changeTrianglesOwnership
  allToAll (bodyForce_F,  &sendVerticesListReversed, &recvVerticesListReversed, 1);
#endif
  allToAll (bodyForce_F,  sendVerticesList_F, recvVerticesList_F, 1);
}

void exportBeta(double * beta_F) {
  std::fill(beta_F, beta_F + nCells_F, 0.);
  for (int index = 0; index < nVertices; index++) {
    int fCell = vertexToFCell[index];
    beta_F[fCell] = betaData[index] * unit_length;
  }
#ifdef  changeTrianglesOwnership
  allToAll (beta_F,  &sendCellsListReversed, &recvCellsListReversed, 1);
#endif
  allToAll (beta_F,  sendCellsList_F, recvCellsList_F, 1);
}

void createReducedMPI(int nLocalEntities, MPI_Comm& reduced_comm_id) {
  int numProcs, me;
  if (reduced_comm_id != MPI_COMM_NULL)
    MPI_Comm_free(&reduced_comm_id);
  MPI_Group world_group_id, reduced_group_id;
  MPI_Comm_size(comm, &numProcs);
  MPI_Comm_rank(comm, &me);
  std::vector<int> haveElements(numProcs);
  int nonEmpty = int(nLocalEntities > 0);
  MPI_Allgather(&nonEmpty, 1, MPI_INT, &haveElements[0], 1, MPI_INT, comm);
  std::vector<int> ranks;
  reduced_ranks.assign(numProcs,-1);
  for (int i = 0; i < numProcs; i++) {
    if (haveElements[i]) {
      reduced_ranks[i] = ranks.size();
      ranks.push_back(i);
    }
  }

  MPI_Comm_group(comm, &world_group_id);
  MPI_Group_incl(world_group_id, ranks.size(), &ranks[0], &reduced_group_id);
  MPI_Comm_create(comm, reduced_group_id, &reduced_comm_id);
}

void getProcIds(std::vector<int>& field, int const * recvArray) {
  int me;
  MPI_Comm_rank(comm, &me);
  field.assign(field.size(), me);

  //unpack recvArray and set the proc rank into filed
  for (int i(1), procID, size; i < recvArray[0]; i += size) {
    procID = recvArray[i++];
    size = recvArray[i++];
    if (procID == me)
      continue;
    for (int k = i; k < i + size; k++)
      field[recvArray[k]] = procID;
  }
}

void getProcIds(std::vector<int>& field, exchangeList_Type const * recvList) {
  int me;
  MPI_Comm_rank(comm, &me);
  field.assign(field.size(), me);
  exchangeList_Type::const_iterator it;

  for (it = recvList->begin(); it != recvList->end(); ++it) {
    if (it->procID == me)
      continue;
    for (int k = 0; k < (int) it->vec.size(); k++)
      field[it->vec[k]] = it->procID;
  }
}

exchangeList_Type unpackMpiArray(int const * array) {
  exchangeList_Type list;
  for (int i(1), procID, size; i < array[0]; i += size) {
    procID = array[i++];
    size = array[i++];
    list.push_back(exchange(procID, &array[i], &array[i + size]));
  }
  return list;
}

void allToAll(std::vector<int>& field, int const * sendArray,
    int const * recvArray, int fieldDim) {
  exchangeList_Type sendList, recvList;

  //unpack sendArray and build the sendList class
  for (int i(1), procID, size; i < sendArray[0]; i += size) {
    procID = sendArray[i++];
    size = sendArray[i++];
    sendList.push_back(
        exchange(procID, &sendArray[i], &sendArray[i + size], fieldDim));
  }

  //unpack recvArray and build the recvList class
  for (int i(1), procID, size; i < recvArray[0]; i += size) {
    procID = recvArray[i++];
    size = recvArray[i++];
    recvList.push_back(
        exchange(procID, &recvArray[i], &recvArray[i + size], fieldDim));
  }

  int me;
  MPI_Comm_rank(comm, &me);

  exchangeList_Type::iterator it;
  for (it = recvList.begin(); it != recvList.end(); ++it) {
    if (it->procID == me)
      continue;
    MPI_Irecv(&(it->buffer[0]), it->buffer.size(), MPI_INT, it->procID,
        it->procID, comm, &it->reqID);
  }

  for (it = sendList.begin(); it != sendList.end(); ++it) {
    if (it->procID == me)
      continue;
    for (ID i = 0; i < it->vec.size(); i++)
      for (int iComp = 0; iComp < fieldDim; iComp++)
        it->buffer[fieldDim * i + iComp] = field[fieldDim * it->vec[i] + iComp];

    MPI_Isend(&(it->buffer[0]), it->buffer.size(), MPI_INT, it->procID, me,
        comm, &it->reqID);
  }

  for (it = recvList.begin(); it != recvList.end(); ++it) {
    if (it->procID == me)
      continue;
    MPI_Wait(&it->reqID, MPI_STATUS_IGNORE);

    for (int i = 0; i < int(it->vec.size()); i++)
      for (int iComp = 0; iComp < fieldDim; iComp++)
        field[fieldDim * it->vec[i] + iComp] = it->buffer[fieldDim * i + iComp];
  }

  for (it = sendList.begin(); it != sendList.end(); ++it) {
    if (it->procID == me)
      continue;
    MPI_Wait(&it->reqID, MPI_STATUS_IGNORE);
  }
}

void allToAll(std::vector<int>& field, exchangeList_Type const * sendList,
    exchangeList_Type const * recvList, int fieldDim) {
  int me;
  MPI_Comm_rank(comm, &me);

  for (int iComp = 0; iComp < fieldDim; iComp++) {
    exchangeList_Type::const_iterator it;
    for (it = recvList->begin(); it != recvList->end(); ++it) {
      if (it->procID == me)
        continue;
      MPI_Irecv(&(it->buffer[0]), it->buffer.size(), MPI_INT, it->procID,
          it->procID, comm, &it->reqID);
    }

    for (it = sendList->begin(); it != sendList->end(); ++it) {
      if (it->procID == me)
        continue;
      for (ID i = 0; i < it->vec.size(); i++)
        it->buffer[i] = field[fieldDim * it->vec[i] + iComp];

      MPI_Isend(&(it->buffer[0]), it->buffer.size(), MPI_INT, it->procID, me,
          comm, &it->reqID);
    }

    for (it = recvList->begin(); it != recvList->end(); ++it) {
      if (it->procID == me)
        continue;
      MPI_Wait(&it->reqID, MPI_STATUS_IGNORE);

      for (int i = 0; i < int(it->vec.size()); i++)
        field[fieldDim * it->vec[i] + iComp] = it->buffer[i];
    }

    for (it = sendList->begin(); it != sendList->end(); ++it) {
      if (it->procID == me)
        continue;
      MPI_Wait(&it->reqID, MPI_STATUS_IGNORE);
    }
  }
}

void allToAll(double* field, exchangeList_Type const * sendList,
    exchangeList_Type const * recvList, int fieldDim) {
  int me;
  MPI_Comm_rank(comm, &me);

  for (int iComp = 0; iComp < fieldDim; iComp++) {
    exchangeList_Type::const_iterator it;

    for (it = recvList->begin(); it != recvList->end(); ++it) {
      if (it->procID == me)
        continue;
      MPI_Irecv(&(it->doubleBuffer[0]), it->doubleBuffer.size(), MPI_DOUBLE,
          it->procID, it->procID, comm, &it->reqID);
    }

    for (it = sendList->begin(); it != sendList->end(); ++it) {
      if (it->procID == me)
        continue;
      for (ID i = 0; i < it->vec.size(); i++)
        it->doubleBuffer[i] = field[fieldDim * it->vec[i] + iComp];

      MPI_Isend(&(it->doubleBuffer[0]), it->doubleBuffer.size(), MPI_DOUBLE,
          it->procID, me, comm, &it->reqID);
    }

    for (it = recvList->begin(); it != recvList->end(); ++it) {
      if (it->procID == me)
        continue;
      MPI_Wait(&it->reqID, MPI_STATUS_IGNORE);

      for (int i = 0; i < int(it->vec.size()); i++)
        field[fieldDim * it->vec[i] + iComp] = it->doubleBuffer[i];
    }

    for (it = sendList->begin(); it != sendList->end(); ++it) {
      if (it->procID == me)
        continue;
      MPI_Wait(&it->reqID, MPI_STATUS_IGNORE);
    }
  }
}

int initialize_iceProblem(int nTriangles) {
  bool keep_proc = nTriangles > 0;

  createReducedMPI(keep_proc, reducedComm);

  isDomainEmpty = !keep_proc;

  // initialize ice problem pointer
  if (keep_proc) {
    std::cout << nTriangles
        << " elements of the triangular grid are stored on this processor"
        << std::endl;
  } else {
    std::cout
        << "No elements of the triangular grid are stored on this processor"
        << std::endl;
  }

  return 0;
}

//barycentric coordinates bcoords are properly updated only when functions return true.
bool belongToTria(double const* x, double const* t, double bcoords[3], double eps) {
  double v1[2],v2[2],v3[2];
  for(int i=0; i<2; ++i) {
  v1[i] = t[i + 2 * 0];
  v2[i] = t[i + 2 * 1];
  v3[i] = t[i + 2 * 2];
  }
  double det = (v3[1]-v2[1])*(v3[0]-v1[0]) - (v3[0]-v2[0])*(v3[1]-v1[1]);
  double c1,c2;

  bcoords[0] = ((v3[1]-v2[1])*(v3[0]-x[0]) - (v3[0]-v2[0])*(v3[1]-x[1]))/det;
  bcoords[1] = (-(v3[1]-v1[1])*(v3[0]-x[0]) + (v3[0]-v1[0])*(v3[1]-x[1]))/det;
  bcoords[2] = 1.0 - bcoords[0] - bcoords[1];

  return ( bcoords[0] > -eps) &&
         ( bcoords[1] > -eps) &&
         ( bcoords[2] > -eps );
}

  void procsSharingVertex(const int vertex, std::vector<int>& procIds) {
    int fCell = vertexToFCell[vertex];
    procIds.clear();
    int nEdg = nEdgesOnCells_F[fCell];
    int me;
    MPI_Comm_rank(comm, &me);
    for(int i=0; i<nEdg; ++i) {
      int fVertex = verticesOnCell_F[maxNEdgesOnCell_F * fCell + i]-1;
      if (verticesMask_F[fVertex] & dynamic_ice_bit_value) {
        int proc = trianglesProcIds[fVertex];
        if(proc != me)
          procIds.push_back(reduced_ranks[proc]);
      }
    }
  }

  void write_ascii_mesh(int const* indexToCellID_F,
    double const* bedTopography_F, double const* lowerSurface_F,
    double const* beta_F, double const* temperature_F,
    double const* surfaceAirTemperature_F, double const* basalHeatFlux_F,
    double const* stiffnessFactor_F,
    double const* effecPress_F, double const* muFriction_F,
    double const* thickness_F, double const* thicknessUncertainty_F,
    double const* smb_F, double const* smbUncertainty_F,
    double const* bmb_F, double const* bmbUncertainty_F,
    double const* observedSurfaceVelocityX_F, double const* observedSurfaceVelocityY_F,
    double const* observedSurfaceVelocityUncertainty_F, 
    double const* observedThicknessTendency_F, double const * observedThicknessTendencyUncertainty_F) {


    // Write out ASCII format

    std::cout << "Writing mesh to albany.msh." << std::endl;
    // msh file
    std::ofstream outfile;
    outfile.precision(15);
    std::stringstream name;
    int me;
    MPI_Comm_rank(comm, &me);
    name <<  "albany.msh." <<  me;
    outfile.open (name.str(), std::ios::out | std::ios::trunc);
    if (outfile.is_open()) {

      //creating set from vector so that we can find elements in the set in log time
      //boundary edges labels: 2 if ice margin, 1 otherwise
       std::vector<int> bdEdgesLabels(numBoundaryEdges,1);
       for (int i = 0; i < iceMarginEdgesLIds.size(); i++)
         bdEdgesLabels[iceMarginEdgesLIds[i]] = 2;

       std::vector<int> boundaryEdges; //list of edge vertices and edge label
       boundaryEdges.resize(3 * numBoundaryEdges);
       for (int index = 0; index < numBoundaryEdges; index++) {
         boundaryEdges[0 + 3 * index] = verticesOnEdge[0 + 2 * index];
         boundaryEdges[1 + 3 * index] = verticesOnEdge[1 + 2 * index];
         boundaryEdges[2 + 3 * index] = bdEdgesLabels[index];
       }

       outfile << "Format: " << 1 << "\n";  // first line stating the format (we are going to provide global ids)
       outfile << "Triangle " << 3 << "\n";  // second line saying it is a mesh of triangles
       outfile << nVertices << " " << nTriangles << " " << numBoundaryEdges << "\n";  // second line

       if(!std::is_sorted(indexToVertexID.begin(), indexToVertexID.end())) {
         std::cout << "ERROR: Global vertices IDs need to be sorted in the mesh" << std::endl;
         exit(1);
         //note, if not sorted, one has to sort the vertices before storing them in the mesh,
         //also, one has to reorder accordingly all the fields defined on vertices
         //(similarly to what is done for temperature, defined on triangles)
       }

       for (int index = 0; index < nVertices; index++) { //coordinates lines
          int iCell = vertexToFCell[index];
          //setting boundary vertices labels, 2 for dirichlet nodes, 1 otherwise
          bool isDirichletVertex = dirichletCellsMask_F[(nLayers+1)*iCell] != 0;
          int vertexLabel = (!isVertexBoundary[index]) ? 0 :
                            isDirichletVertex ? 3 : 1;

          outfile << indexToVertexID[index] << " " << xCell_F[iCell] / unit_length << " " << yCell_F[iCell] / unit_length << " " << vertexLabel << "\n"  ;
       }

       // sort triangle IDs (needed by Albany)
       std::vector<int> sortingIndex;
       computeSortingIndices(sortingIndex, indexToTriangleID, nTriangles);

       for (int iTria = 0; iTria < nTriangles; iTria++) {//triangles lines
         int index = sortingIndex[iTria];
         outfile << indexToTriangleID[index] << " " << verticesOnTria[0 + 3 * index] + 1 << " " << verticesOnTria[1 + 3 * index] + 1 << " " << verticesOnTria[2 + 3 * index] + 1 << " " << 1 << "\n"; // last digit can be used to specify a 'material'.  Not used by Albany LandIce, so giving dummy value
       }

       // sort edges IDs (needed by Albany)
       computeSortingIndices(sortingIndex, indexToEdgeID, numBoundaryEdges);

       for (int iEdge = 0; iEdge < numBoundaryEdges; iEdge++) { // boundary edges lines
         int index = sortingIndex[iEdge];
         outfile <<  indexToEdgeID[index] << " " << boundaryEdges[0 + 3 * index] + 1 << " " << boundaryEdges[1 + 3 * index] + 1 << " " << boundaryEdges[2 + 3 * index] << "\n"; //last digit can be used to tell whether it's floating or not.. but let's worry about this later.
       }

       outfile.close();
       }
    else {
       std::cout << "Failed to open mshfile!"<< std::endl;
    }

    // individual field values
    // Call needed functions to process MPAS fields to Albany units/format

    std::vector<std::pair<int, int> > marineBdyExtensionMap;  // local map to be created by importFields
    importFields(marineBdyExtensionMap, bedTopography_F, lowerSurface_F, thickness_F, beta_F,
                   stiffnessFactor_F, effecPress_F, muFriction_F, temperature_F, smb_F, minThickness);

    import2DFieldsObservations(marineBdyExtensionMap,
                    thicknessUncertainty_F,
                    smbUncertainty_F,
                    bmb_F, bmbUncertainty_F,
                    observedSurfaceVelocityX_F, observedSurfaceVelocityY_F, observedSurfaceVelocityUncertainty_F,
                    observedThicknessTendency_F, observedThicknessTendencyUncertainty_F,
                    surfaceAirTemperature_F, basalHeatFlux_F,
                    indexToCellID_F);

    // apparent mass balance
    std::vector<double> appMbData(smbData.size()),
                        appMbUncertaintyData(smbData.size());
 
    for (int i=0; i<smbData.size(); ++i) {
      appMbData[i] = smbData[i]+bmbData[i]- observedDHDtData[i];
      //assuming fields are uncorrelated
      double variance = std::pow(smbUncertaintyData[i],2)+std::pow(bmbUncertaintyData[i],2)+std::pow(observedDHDtUncertaintyData[i],2);
      appMbUncertaintyData[i] =std::sqrt(variance);
    }
   

    // Write out individual fields
    write_ascii_mesh_field_int(indexToCellIDData, "mpas_cellID");
    write_ascii_mesh_field(thicknessData, "thickness");
    write_ascii_mesh_field(thicknessUncertaintyData, "thickness_uncertainty");

    write_ascii_mesh_field(elevationData, "surface_height");
    write_ascii_mesh_field(bedTopographyData, "bed_topography");

    write_ascii_mesh_field(surfaceAirTemperatureData, "surface_air_temperature");
    write_ascii_mesh_field(basalHeatFluxData, "basal_heat_flux");

    write_ascii_mesh_field(betaData, "basal_friction");

    write_ascii_mesh_field(smbData, "surface_mass_balance");
    write_ascii_mesh_field(smbUncertaintyData, "surface_mass_balance_uncertainty");

    write_ascii_mesh_field(bmbData, "basal_mass_balance");
    write_ascii_mesh_field(bmbUncertaintyData, "basal_mass_balance_uncertainty");
    
    write_ascii_mesh_field(observedDHDtData, "dhdt");
    write_ascii_mesh_field(observedDHDtUncertaintyData, "dhdt_uncertainty");
    
    write_ascii_mesh_field(appMbData, "apparent_mass_balance");
    write_ascii_mesh_field(appMbUncertaintyData, "apparent_mass_balance_uncertainty");
    
    write_ascii_mesh_field(observedVeloUncertaintyData, "surface_velocity_uncertainty");

    write_ascii_mesh_field(effecPressData, "effective_pressure");

    // These two fields are more complicated than the others so cannot use the function to write out

    std::cout << "Writing temperature.ascii." << std::endl;
    outfile.open ("temperature.ascii", std::ios::out | std::ios::trunc);
    if (outfile.is_open()) {
      outfile << nTriangles << " " << nLayers << "\n";  //number of triangles and number of layers on first line

      double midLayer = 0;
      for (int il = 0; il < nLayers; il++) { //sigma coordinates for temperature
        midLayer += layersRatio[il]/2.0;
        outfile << midLayer << "\n";
        midLayer += layersRatio[il]/2.0;
      }

      // sort triangle IDs (needed by Albany)
      std::vector<int> sortingIndex;
      computeSortingIndices(sortingIndex,indexToTriangleID,nTriangles);

      int lElemColumnShift = (Ordering == 1) ? 1 : nTriangles;
      int elemLayerShift = (Ordering == 0) ? 1 : nLayers;
      for(int il = 0; il<nLayers; ++il)
        for(int i = 0; i<nTriangles; ++i)  {//temperature values layer by layer
          int index = sortingIndex[i];
          outfile << temperatureDataOnPrisms[index*elemLayerShift + il*lElemColumnShift]<<"\n";
        }

      outfile.close();
    }
    else {
       std::cout << "Failed to open temperature ascii file!"<< std::endl;
    }

    std::cout << "Writing surface_velocity.ascii." << std::endl;
    outfile.open ("surface_velocity.ascii", std::ios::out | std::ios::trunc);
    if (outfile.is_open()) {
       outfile << nVertices << " " << 2 << "\n";  //number of vertices, number of components per vertex
       for(int i = 0; i<nVertices; ++i)
         outfile << observedVeloXData[i] << "\n";
       for(int i = 0; i<nVertices; ++i)
         outfile << observedVeloYData[i] << "\n";
       outfile.close();
    }
    else {
       std::cout << "Failed to open surface velocity ascii file!"<< std::endl;
    }

    //here we save the surface velocity as an extruded field,
    //having two levels at sigma coords 0 and 1
    //this enables Albany to prescribe surface velocity as Dirichlet BC
    std::cout << "Writing extruded_surface_velocity.ascii." << std::endl;
    outfile.open ("extruded_surface_velocity.ascii", std::ios::out | std::ios::trunc);
    if (outfile.is_open()) {
       outfile << nVertices << " " << 2 << " " << 2 << "\n";  //number of vertices, number of levels, number of components per vertex
       outfile << 0.0 << "\n";
       outfile << 1.0 << "\n";  // sigma coordinates for velocity
       for(int il=0; il<2; ++il) {
         for(int i = 0; i<nVertices; ++i)
           outfile << observedVeloXData[i] << "\n";
         for(int i = 0; i<nVertices; ++i)
           outfile << observedVeloYData[i] << "\n";
       }
       outfile.close();
    }
    else {
       std::cout << "Failed to open surface velocity ascii file!"<< std::endl;
    }


    std::cout << "\nWriting of all Albany fields complete." << std::endl;

  }


  // compute sortingIndices, so that vectorToSort[sortingIndices[i]], i=0, 1, ..., numIndices-1 is sorted
  void computeSortingIndices(std::vector<int>& sortingIndices, const std::vector<int>& vectorToSort, int numIndices) {
    // sort triangle IDs (needed by Albany)
    sortingIndices.resize(numIndices);

    // sortingIndices = [0,1,2,..,nTriangles-1];
    for (int i = 0; i < numIndices; i++)
      sortingIndices[i] = i;

    //compute sortingIndices that makes indexToTriangleID sorted
    std::sort(sortingIndices.begin(), sortingIndices.end(), [&](int il, int ir) {
            return (vectorToSort[il] < vectorToSort[ir]);
        });
  }


  void write_ascii_mesh_field(std::vector<double> fieldData, std::string filenamebase) {

    std::string filename = filenamebase+".ascii";
    std::cout << "Writing " << filename << std::endl;
    std::ofstream outfile;
    outfile.precision(15);
    outfile.open (filename.c_str(), std::ios::out | std::ios::trunc);
    if (outfile.is_open()) {
       outfile << nVertices << "\n";  //number of vertices on first line
       for(int i = 0; i < nVertices; ++i)
         outfile << fieldData[i] << "\n";
       outfile.close();
    }
    else {
       std::cout << "Error: Failed to open  "+filename << std::endl;
    }
  }


  void write_ascii_mesh_field_int(std::vector<int> fieldData, std::string filenamebase) {

    std::string filename = filenamebase+".ascii";
    std::cout << "Writing " << filename << std::endl;
    std::ofstream outfile;
    outfile.open (filename.c_str(), std::ios::out | std::ios::trunc);
    if (outfile.is_open()) {
       outfile << nVertices << "\n";  //number of vertices on first line
       for(int i = 0; i < nVertices; ++i)
         outfile << fieldData[i] << "\n";
       outfile.close();
    }
    else {
       std::cout << "Error: Failed to open  "+filename << std::endl;
    }
  }

