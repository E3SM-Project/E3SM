/*
Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
and the University Corporation for Atmospheric Research (UCAR).

Unless noted otherwise source code is licensed under the BSD license.
Additional copyright and license information can be found in the LICENSE file
distributed with this code, or at http://mpas-dev.github.io/license.html
*/

// ===================================================
//! Includes
// ===================================================
#include <cstring>
#include <vector>
#include <mpi.h>
#include <list>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <fcntl.h>
#include <limits>
#include <cmath>
#include <map>

#ifndef MPASLI_EXTERNAL_INTERFACE_DISABLE_MANGLING
#define velocity_solver_init_mpi velocity_solver_init_mpi_
#define velocity_solver_finalize velocity_solver_finalize_
#define velocity_solver_init_l1l2 velocity_solver_init_l1l2_
#define velocity_solver_solve_l1l2 velocity_solver_solve_l1l2_
#define velocity_solver_init_fo velocity_solver_init_fo_
#define velocity_solver_solve_fo velocity_solver_solve_fo_
#define velocity_solver_init_stokes velocity_solver_init_stokes_
#define velocity_solver_solve_stokes velocity_solver_solve_stokes_
#define velocity_solver_compute_2d_grid velocity_solver_compute_2d_grid_
#define velocity_solver_set_grid_data velocity_solver_set_grid_data_
#define velocity_solver_extrude_3d_grid velocity_solver_extrude_3d_grid_
#define velocity_solver_export_l1l2_velocity velocity_solver_export_l1l2_velocity_
#define velocity_solver_export_2d_data velocity_solver_export_2d_data_
#define velocity_solver_export_fo_velocity velocity_solver_export_fo_velocity_
#define velocity_solver_estimate_SS_SMB velocity_solver_estimate_ss_smb_
#define interface_init_log interface_init_log_
#define interface_redirect_stdout interface_redirect_stdout_
#define interface_reset_stdout interface_reset_stdout_
#define write_ascii_mesh write_ascii_mesh_
#endif

struct exchange {
  const int procID;
  const std::vector<int> vec;
  mutable std::vector<int> buffer;
  mutable std::vector<double> doubleBuffer;
  mutable MPI_Request reqID;

  exchange(int _procID, int const* vec_first, int const* vec_last,
      int fieldDim = 1);
};

typedef std::list<exchange> exchangeList_Type;

typedef unsigned int ID;
typedef unsigned int UInt;
const ID NotAnId = std::numeric_limits<int>::max();

// ===================================================
//! Interface function
// ===================================================
extern "C" {

int velocity_solver_init_mpi(int* fComm);

void velocity_solver_finalize();

void velocity_solver_set_parameters(double const* gravity_F, double const* ice_density_F, double const* ocean_density_F,
                        double const* sea_level_F, double const* flowParamA_F,
                        double const* flowLawExponent_F, double const* dynamic_thickness_F,
                        double const* clausius_clapeyron_coeff,
                        double const* thermal_thickness_limit_F,
                        int const* li_mask_ValueDynamicIce, int const* li_mask_ValueIce,
                        bool const* use_GLP_F);

void velocity_solver_init_l1l2(double const* levelsRatio);

void velocity_solver_init_fo(double const* levelsRatio);

void velocity_solver_solve_l1l2(double const* lowerSurface_F,
    double const* thickness_F, double const* beta_F, double const* temperature_F,
    double* const dirichletVelocityXValue = 0, double* const dirichletVelocitYValue = 0,
    double* u_normal_F = 0,
    double* xVelocityOnCell = 0, double* yVelocityOnCell = 0);

void velocity_solver_solve_fo(double const* bedTopography_F, double const* lowerSurface_F,
    double const* thickness_F, double * beta_F, double const* smb_F, double const* temperature_F, double const* stiffnessFactor_F,
    double const* effecPress_F, double const* muFriction_F,
    double* const dirichletVelocityXValue = 0, double* const dirichletVelocitYValue = 0,
    double* u_normal_F = 0, double* bodyForce_F = 0, double* dissipation_heat_F = 0,
    double* xVelocityOnCell = 0, double* yVelocityOnCell = 0, double const * deltat = 0,
    int *error = 0 );


void velocity_solver_compute_2d_grid(int const* verticesMask_F, int const* _cellsMask_F, int const* dirichletNodesMask_F);

void velocity_solver_set_grid_data(int const* _nCells_F, int const* _nEdges_F,
    int const* _nVertices_F, int const* _nLayers, int const* _nCellsSolve_F,
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
    int const* sendVerticesArray_F, int const* recvVerticesArray_F);

void velocity_solver_extrude_3d_grid(double const* levelsRatio_F);

void velocity_solver_export_l1l2_velocity();

void velocity_solver_export_fo_velocity();

void interface_init_log();

void interface_redirect_stdout(int const* iTimestep);

void interface_reset_stdout();

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
    double const* observedVelocityUncertainty_F, 
    double const* observedThicknessTendency_F, double const * observedThicknessTendencyUncertainty_F);

} // extern "C"

extern int velocity_solver_init_mpi__(MPI_Comm comm);
extern void velocity_solver_finalize__();

extern void velocity_solver_set_physical_parameters__(double const& gravity, double const& ice_density, double const& ocean_density, double const& sea_level, double const& flowParamA, 
                        double const& flowLawExponent, double const& dynamic_thickness, bool const& useGLP, double const& clausiusClapeyronCoeff); 

extern void velocity_solver_solve_fo__(int nLayers, int nGlobalVertices,
    int nGlobalTriangles, bool ordering, bool first_time_step,
    const std::vector<int>& indexToVertexID,
    const std::vector<int>& indexToTriangleID, double minBeta,
    const std::vector<double>& regulThk,
    const std::vector<double>& levelsNormalizedThickness,
    const std::vector<double>& elevationData,
    const std::vector<double>& thicknessData,
          std::vector<double>& betaData,
    const std::vector<double>& bedTopographyData,
    const std::vector<double>& smbData,
    const std::vector<double>& stiffnessFactorData,
    const std::vector<double>& effecPressData,
    const std::vector<double>& muFrictionData,
    const std::vector<double>& temperatureDataOnPrisms,
    std::vector<double>& bodyForceOnBasalCell,
    std::vector<double>& dissipationHeatOnPrisms,
    std::vector<double>& velocityOnVertices,
    int& error,
    const double& deltat = 0.0);

extern void velocity_solver_compute_2d_grid__(MPI_Comm);


extern void velocity_solver_export_2d_data__(MPI_Comm reducedComm,
    const std::vector<double>& elevationData,
    const std::vector<double>& thicknessData,
    const std::vector<double>& betaData,
    const std::vector<int>& indexToVertexID);

extern void velocity_solver_extrude_3d_grid__(
    int nLayers, int globalTriangleStride, int globalVertexStride, int globalEdgeStride,
    int Ordering, MPI_Comm reducedComm,
    const std::vector<int>& indexToVertexID,
    const std::vector<int>& vertexProcIDs,
    const std::vector<double>& verticesCoords,
    const std::vector<int>& verticesOnTria,
    const std::vector<std::vector<int>> procsSharingVertices,
    const std::vector<bool>& isBoundaryEdge,
    const std::vector<int>& trianglesOnEdge,
    const std::vector<int>& verticesOnEdge,
    const std::vector<int>& indexToEdgeID,
    const std::vector<int>& indexToTriangleID,
    const std::vector<int>& dirichletNodes,
    const std::vector<int>& iceMarginEdgesID);

extern void velocity_solver_export_fo_velocity__(MPI_Comm reducedComm);

exchangeList_Type unpackMpiArray(int const* array);

double signedTriangleArea(const double* x, const double* y);

double signedTriangleArea(const double* x, const double* y, const double* z);

void createReducedMPI(int nLocalEntities, MPI_Comm& reduced_comm_id);

void importFields(std::vector<std::pair<int, int> >& marineBdyExtensionMap,
                double const* bedTopography_F, double const* lowerSurface_F, double const* thickness_F,
    double const* beta_F = 0, double const* stiffnessFactor_F = 0, double const* effecPress_F = 0, double const* muFriction_F = 0, double const* temperature_F = 0, double const* smb_F = 0, double eps = 0);

void import2DFieldsObservations(std::vector<std::pair<int, int> >& marineBdyExtensionMap,
            double const * lowerSurface_F, 
            double const * thickness_F, double const * thicknessUncertainty_F,
            double const * smbUncertainty_F,
            double const * bmb_F, double const * bmbUncertainty_F,
            double const * observedSurfaceVelocityX_F, double const * observedSurfaceVelocityY_F,
            double const * observedSurfaceVelocityUncertainty_F,
            double const * observedThicknessTendency_F, double const * observedThicknessTendencyUncertainty_F,
            double const* surfaceAirTemperature_F, double const* basalHeatFlux_F,
            int const * indexToCellID_F);

void computeSortingIndices(std::vector<int>& sortingIndices, const std::vector<int>& vectorToSort, int numIndices);

void write_ascii_mesh_field(std::vector<double> fieldData, std::string filenamebase);

void write_ascii_mesh_field_int(std::vector<int> fieldData, std::string filenamebase);

std::vector<int> extendMaskByOneLayer(int const* verticesMask_F);

void exportDissipationHeat(double * dissipationHeat_F);
void exportBodyForce(double * bodyForce_F);
void exportBeta(double * beta_F);

void get_prism_velocity_on_FEdges(double* uNormal,
    const std::vector<double>& velocityOnCells,
    const std::vector<int>& edgeToFEdge);

int initialize_iceProblem(int nTriangles);

void createReverseExchangeLists(exchangeList_Type& sendListReverse_F,
    exchangeList_Type& receiveListReverse_F,
    const std::vector<int>& newProcIds, const int* indexToID_F, exchangeList_Type const * recvList_F);

void mapCellsToVertices(const std::vector<double>& velocityOnCells,
    std::vector<double>& velocityOnVertices, int fieldDim, int numLayers,
    int ordering);

void mapVerticesToCells(const std::vector<double>& velocityOnVertices,
    double* velocityOnCells, int fieldDim, int numLayers, int ordering);

void getProcIds(std::vector<int>& field, int const* recvArray);

void getProcIds(std::vector<int>& field, exchangeList_Type const* recvList);

void allToAll(std::vector<int>& field, int const* sendArray,
    int const* recvArray, int fieldDim = 1);

void allToAll(std::vector<int>& field, exchangeList_Type const* sendList,
    exchangeList_Type const* recvList, int fieldDim = 1);

void allToAll(double* field, exchangeList_Type const* sendList,
    exchangeList_Type const* recvList, int fieldDim = 1);

void procsSharingVertex(const int vertex, std::vector<int>& procIds);

bool belongToTria(double const* x, double const* t, double bcoords[3], double eps = 1e-3);



