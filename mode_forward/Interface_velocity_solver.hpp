/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s):
 Date: 2009-03-24

 Copyright (C) 2009 EPFL

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA
 */

// ===================================================
//! Includes
// ===================================================
//#include <boost/program_options.hpp>
#include <vector>
#include <mpi.h>
#include <list>
#include <iostream>
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
#endif

//#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
//#include <lifev/core/algorithm/PreconditionerML.hpp>

//#include <lifev/ice_sheet/solver/IceProblem.hpp>

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

void velocity_solver_set_parameters(double const* gravity_F, double const* ice_density_F, double const* ocean_density_F, double const* sea_level_F, double const* flowParamA_F, 
                        double const* enhancementFactor_F, double const* flowLawExponent_F, double const* dynamic_thickness_F, int const* li_mask_ValueDynamicIce, int const* li_mask_ValueIce);

void velocity_solver_init_l1l2(double const* levelsRatio);

void velocity_solver_init_fo(double const* levelsRatio);

void velocity_solver_solve_l1l2(double const* lowerSurface_F,
    double const* thickness_F, double const* beta_F, double const* temperature_F,
    double* const dirichletVelocityXValue = 0, double* const dirichletVelocitYValue = 0,
    double* u_normal_F = 0,
    double* xVelocityOnCell = 0, double* yVelocityOnCell = 0);

void velocity_solver_solve_fo(double const* bedTopography_F, double const* lowerSurface_F,
    double const* thickness_F, double const* beta_F, double const* smb_F, double const* temperature_F,
    double* const dirichletVelocityXValue = 0, double* const dirichletVelocitYValue = 0,
    double* u_normal_F = 0, /*double* dissipation_heat_F = 0,*/
    double* xVelocityOnCell = 0, double* yVelocityOnCell = 0, double const * deltat = 0);


void velocity_solver_compute_2d_grid(int const* verticesMask_F, int const* _cellsMask_F, int const* dirichletNodesMask_F, int const* floatingEdgeMask_F);

void velocity_solver_set_grid_data(int const* _nCells_F, int const* _nEdges_F,
    int const* _nVertices_F, int const* _nLayers, int const* _nCellsSolve_F,
    int const* _nEdgesSolve_F, int const* _nVerticesSolve_F,
    int const* _maxNEdgesOnCell_F, double const* radius_F,
    int const* _cellsOnEdge_F, int const* _cellsOnVertex_F,
    int const* _verticesOnCell_F, int const* _verticesOnEdge_F,
    int const* _edgesOnCell_F, int const* _nEdgesOnCells_F,
    int const* _indexToCellID_F,
    double const* _xCell_F, double const* _yCell_F, double const* _zCell_F,
    double const* _xVertex_F, double const* _yVertex_F, double const* _zVertex_F,
    double const* _areaTriangle_F,
    int const* sendCellsArray_F, int const* recvCellsArray_F,
    int const* sendEdgesArray_F, int const* recvEdgesArray_F,
    int const* sendVerticesArray_F, int const* recvVerticesArray_F);

void velocity_solver_extrude_3d_grid(double const* levelsRatio_F,
    double const* lowerSurface_F, double const* thickness_F);

void velocity_solver_export_l1l2_velocity();

void velocity_solver_export_fo_velocity();

//void velocity_solver_estimate_SS_SMB (const double* u_normal_F, double* sfcMassBal);

}

extern void velocity_solver_finalize__();

#ifdef LIFEV
extern void velocity_solver_init_l1l2__(const std::vector<double>& layersRatio, const std::vector<double>& velocityOnVertices, bool initialize_velocity);

extern void velocity_solver_solve_l1l2__(const std::vector<double>& elevationData,
    const std::vector<double>& thicknessData, const std::vector<double>& betaData,
    const std::vector<double>& temperatureData, const std::vector<int>& indexToVertexID,
    std::vector<double>& velocityOnVertices);

extern void velocity_solver_init_fo__(const std::vector<double>& layersRatio, const std::vector<double>& velocityOnVertices, const std::vector<int>& indexToVertexID, bool initialize_velocity);

extern void velocity_solver_export_l1l2_velocity__(const std::vector<double>& layersRatio, const std::vector<double>& elevationData, const std::vector<double>& regulThk,
    const std::vector<int>& mpasIndexToVertexID, MPI_Comm reducedComm);

#endif

extern void velocity_solver_set_physical_parameters__(double const& gravity, double const& ice_density, double const& ocean_density, double const& sea_level, double const& flowParamA, 
                        double const& enhancementFactor, double const& flowLawExponent, double const& dynamic_thickness); 

extern void velocity_solver_solve_fo__(int nLayers, int nGlobalVertices,
    int nGlobalTriangles, bool ordering, bool first_time_step,
    const std::vector<int>& indexToVertexID,
    const std::vector<int>& indexToTriangleID, double minBeta,
    const std::vector<double>& regulThk,
    const std::vector<double>& levelsNormalizedThickness,
    const std::vector<double>& elevationData,
    const std::vector<double>& thicknessData,
    const std::vector<double>& betaData,
    const std::vector<double>& bedTopographyData,
    const std::vector<double>& smbData,
    const std::vector<double>& temperatureOnTetra,
    std::vector<double>& dissipationHeatOnTetra,
    std::vector<double>& velocityOnVertices,
    const double& deltat = 0.0);


#ifdef LIFEV
extern void velocity_solver_compute_2d_grid__(int nGlobalTriangles,
     int nGlobalVertices, int nGlobalEdges,
     const std::vector<int>& indexToVertexID,
     const std::vector<double>& verticesCoords,
     const std::vector<bool>& isVertexBoundary,
     const std::vector<int>& verticesOnTria,
     const std::vector<bool>& isBoundaryEdge,
     const std::vector<int>& trianglesOnEdge,
     const std::vector<int>& trianglesPositionsOnEdge,
     const std::vector<int>& verticesOnEdge,
     const std::vector<int>& indexToEdgeID,
     const std::vector<int>& indexToTriangleID,
     const std::vector < std::pair<int, int> >& procOnInterfaceEdge);

#else
extern void velocity_solver_compute_2d_grid__(MPI_Comm);
#endif


extern void velocity_solver_export_2d_data__(MPI_Comm reducedComm,
    const std::vector<double>& elevationData,
    const std::vector<double>& thicknessData,
    const std::vector<double>& betaData,
    const std::vector<int>& indexToVertexID);

extern void velocity_solver_extrude_3d_grid__(int nLayers, int nGlobalTriangles,
    int nGlobalVertices, int nGlobalEdges, int Ordering, MPI_Comm reducedComm,
    const std::vector<int>& indexToVertexID,
    const std::vector<int>& mpasIndexToVertexID,
    const std::vector<double>& verticesCoords,
    const std::vector<bool>& isVertexBoundary,
    const std::vector<int>& verticesOnTria,
    const std::vector<bool>& isBoundaryEdge,
    const std::vector<int>& trianglesOnEdge,
    const std::vector<int>& trianglesPositionsOnEdge,
    const std::vector<int>& verticesOnEdge,
    const std::vector<int>& indexToEdgeID,
    const std::vector<int>& indexToTriangleID,
    const std::vector<int>& dirichletNodes,
    const std::vector<int>&floatingEdges);

//extern void velocity_solver_export_l1l2_velocity__();

extern void velocity_solver_export_fo_velocity__(MPI_Comm reducedComm);



#ifdef LIFEV
extern int  velocity_solver_initialize_iceProblem__(bool keep_proc, MPI_Comm reducedComm);
#endif

//extern void velocity_solver_estimate_SS_SMB__ (const double* u_normal_F, double* sfcMassBal);

exchangeList_Type unpackMpiArray(int const* array);

bool isGhostTriangle(int i, double relTol = 1e-1);

double signedTriangleArea(const double* x, const double* y);

double signedTriangleArea(const double* x, const double* y, const double* z);

void createReducedMPI(int nLocalEntities, MPI_Comm& reduced_comm_id);

void import2DFields(double const* bedTopography_F, double const* lowerSurface_F, double const* thickness_F,
    double const* beta_F = 0, double const* smb_F = 0, double eps = 0);

std::vector<int> extendMaskByOneLayer(int const* verticesMask_F);

void extendMaskByOneLayer(int const* verticesMask_F,
    std::vector<int>& extendedFVerticesMask);

void importP0Temperature(double const* temperature_F);

void exportDissipationHeat(double * dissipationHeat_F);

void get_prism_velocity_on_FEdges(double* uNormal,
    const std::vector<double>& velocityOnCells,
    const std::vector<int>& edgeToFEdge);

int initialize_iceProblem(int nTriangles);

void createReverseCellsExchangeLists(exchangeList_Type& sendListReverse_F,
    exchangeList_Type& receiveListReverse_F,
    const std::vector<int>& fVertexToTriangleID,
    const std::vector<int>& fCellToVertexID);

void createReverseEdgesExchangeLists(exchangeList_Type& sendListReverse_F,
    exchangeList_Type& receiveListReverse_F,
    const std::vector<int>& fVertexToTriangleID,
    const std::vector<int>& fEdgeToEdgeID);

void mapCellsToVertices(const std::vector<double>& velocityOnCells,
    std::vector<double>& velocityOnVertices, int fieldDim, int numLayers,
    int ordering);

void mapVerticesToCells(const std::vector<double>& velocityOnVertices,
    double* velocityOnCells, int fieldDim, int numLayers, int ordering);

void computeLocalOffset(int nLocalEntities, int& localOffset,
    int& nGlobalEntities);

void getProcIds(std::vector<int>& field, int const* recvArray);

void getProcIds(std::vector<int>& field, exchangeList_Type const* recvList);

void allToAll(std::vector<int>& field, int const* sendArray,
    int const* recvArray, int fieldDim = 1);

void allToAll(std::vector<int>& field, exchangeList_Type const* sendList,
    exchangeList_Type const* recvList, int fieldDim = 1);

void allToAll(double* field, exchangeList_Type const* sendList,
    exchangeList_Type const* recvList, int fieldDim = 1);

int prismType(long long int const* prismVertexMpasIds, int& minIndex);
void tetrasFromPrismStructured (long long int const* prismVertexMpasIds, long long int const* prismVertexGIds, long long int tetrasIdsOnPrism[][4]);
void computeMap();

void setBdFacesOnPrism (const std::vector<std::vector<std::vector<int> > >& prismStruct, const std::vector<int>& prismFaceIds, std::vector<int>& tetraPos, std::vector<int>& facePos);
void tetrasFromPrismStructured (int const* prismVertexMpasIds, int const* prismVertexGIds, int tetrasIdsOnPrism[][4]);
void procsSharingVertex(const int vertex, std::vector<int>& procIds);

bool belongToTria(double const* x, double const* t, double bcoords[3], double eps = 1e-3);



