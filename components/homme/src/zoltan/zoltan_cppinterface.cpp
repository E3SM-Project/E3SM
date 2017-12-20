
#ifdef HAVE_CONFIG_H
#include "config.h.c"
#endif

#include "zoltan_cppinterface.hpp"
//#define VISUALIZE_MAPPING

#if TRILINOS_HAVE_ZOLTAN2
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_TaskMapping.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_Map.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_GlobalMPISession.hpp>

#include <Zoltan2_MappingProblem.hpp>
#include <Zoltan2_MappingSolution.hpp>
#include <Zoltan2_EvaluatePartition.hpp>
#include <Zoltan2_EvaluateMapping.hpp>
#include <iostream>
#include <vector>
#include <algorithm>

struct SortItem{
  int id;
  int weight;
  bool operator<(const SortItem& a) const
  {
      return this->id < a.id;
  }
};
void sort_graph(
    int *nelem,
    int *xadj,
    int *adjncy,
    double *adjwgt,
    double *vwgt){

  std::vector <SortItem> sorter(*nelem);

  for (int i =0; i < *nelem; ++i){
    int b = xadj[i];
    int e = xadj[i + 1];
    for (int j = b; j < e; ++j){
      sorter[j-b].id = adjncy[j];
      sorter[j-b].weight = adjwgt[j];
    }
    std::sort(sorter.begin(), sorter.begin() + e - b);
    for (int j = b; j < e; ++j){
      adjncy[j] = sorter[j-b].id;
      adjwgt[j] = sorter[j-b].weight;
    }
  }

}

void zoltan_partition_problem(
    int *nelem,
    int *xadj,
    int *adjncy,
    double *adjwgt,
    double *vwgt,
    int *nparts,
    MPI_Comm comm,
    double *xcoord,
    double *ycoord,
    double *zcoord, int *coord_dimension,
    int *result_parts,
    int *partmethod,
    int *mapmethod){
  using namespace Teuchos;
  typedef int zlno_t;
  typedef int zgno_t;
  //typedef int part_t;
  typedef double zscalar_t;

  typedef Tpetra::Map<>::node_type znode_t;
  typedef Tpetra::Map<zlno_t, zgno_t, znode_t> map_t;
  size_t numGlobalCoords = *nelem;


  Teuchos::RCP<const Teuchos::Comm<int> > tcomm =
      Teuchos::RCP<const Teuchos::Comm<int> > (new Teuchos::MpiComm<int> (comm));

  RCP<const map_t> map = rcp (new map_t (numGlobalCoords, 0, tcomm));

  typedef Tpetra::CrsGraph<zlno_t, zgno_t, znode_t> tcrsGraph_t;
  RCP<tcrsGraph_t> TpetraCrsGraph(new tcrsGraph_t (map, 0));

  const zlno_t numMyElements = map->getNodeNumElements ();
  const zgno_t myBegin = map->getGlobalElement (0);

  for (zlno_t lclRow = 0; lclRow < numMyElements; ++lclRow) {
    const zgno_t gblRow = map->getGlobalElement (lclRow);
    zgno_t begin = xadj[gblRow];
    zgno_t end = xadj[gblRow + 1];
    const ArrayView< const zgno_t > indices(adjncy+begin, end-begin);
    TpetraCrsGraph->insertGlobalIndices(gblRow, indices);
  }
  TpetraCrsGraph->fillComplete ();

  RCP<const tcrsGraph_t> const_data = rcp_const_cast<const tcrsGraph_t>(TpetraCrsGraph);
  typedef Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> tMVector_t;
  typedef Zoltan2::XpetraCrsGraphAdapter<tcrsGraph_t, tMVector_t> adapter_t;
  RCP<adapter_t> ia (new adapter_t(const_data/*,(int)vtx_weights.size(),(int)edge_weights.size()*/, 1, 1));

  //for now no edge weights, and no vertex weights.
  //ia->setVertexWeights(vtx_weights[i],vtx_weightStride[i],i);
  //ia->setEdgeWeights(edge_weights[i],edge_weightStride[i],i);


  /***********************************SET COORDINATES*********************/
  const int coord_dim = *coord_dimension;
  // make an array of array views containing the coordinate data
  Teuchos::Array<Teuchos::ArrayView<const zscalar_t> > coordView(coord_dim);

  if(numMyElements > 0){
    Teuchos::ArrayView<const zscalar_t> a(xcoord + myBegin, numMyElements);
    coordView[0] = a;
    Teuchos::ArrayView<const zscalar_t> b(ycoord + myBegin, numMyElements);
    coordView[1] = b;
    if (coord_dim == 3){
      Teuchos::ArrayView<const zscalar_t> c(zcoord + myBegin, numMyElements);
      coordView[2] = c;
    }
  }
  else {
    Teuchos::ArrayView<const zscalar_t> a;
    coordView[0] = a;
    coordView[1] = a;

    if (coord_dim == 3){
      coordView[2] = a;
    }
  }

  RCP<tMVector_t> coords(new tMVector_t(map, coordView.view(0, coord_dim), coord_dim));//= set multivector;
  RCP<const tMVector_t> const_coords = rcp_const_cast<const tMVector_t>(coords);
  Zoltan2::XpetraMultiVectorAdapter<tMVector_t> *adapter = (new Zoltan2::XpetraMultiVectorAdapter<tMVector_t>(const_coords));

  ia->setCoordinateInput(adapter);
  ia->setEdgeWeights(adjwgt, 1, 0);
  ia->setVertexWeights(vwgt, 1, 0);
  /***********************************SET COORDINATES*********************/


  typedef Zoltan2::PartitioningProblem<adapter_t> xcrsGraph_problem_t; // xpetra_graph problem type
  ParameterList zoltan2_parameters;
  zoltan2_parameters.set("compute_metrics", true);
  zoltan2_parameters.set("imbalance_tolerance", "1.0");
  zoltan2_parameters.set("num_global_parts", tcomm->getSize());
  switch (*partmethod){
  case 6:
    zoltan2_parameters.set("algorithm", "rcb");
    break;
  case 7:
    zoltan2_parameters.set("algorithm", "multijagged");
    zoltan2_parameters.set("mj_recursion_depth", "3");
    break;
  case 8:
    zoltan2_parameters.set("algorithm", "multijagged");
    zoltan2_parameters.set("mj_enable_rcb", true);
    break;
  case 9:
    zoltan2_parameters.set("algorithm", "rib");
    break;
  case 10:
    zoltan2_parameters.set("algorithm", "hsfc");
    break;
  case 11:
    zoltan2_parameters.set("algorithm", "patoh");
    break;
  case 12:
    zoltan2_parameters.set("algorithm", "zoltan");
    {
      Teuchos::ParameterList &zparams =
        zoltan2_parameters.sublist("zoltan_parameters",false);
      zparams.set("LB_METHOD", "PHG");
      zparams.set("LB_APPROACH", "PARTITION");
    }
    //zoltan2_parameters.set("algorithm", "phg");
    break;
  case 13:
    zoltan2_parameters.set("algorithm", "metis");
    break;
  case 14:
    zoltan2_parameters.set("algorithm", "parmetis");
    break;
  case 15:
    zoltan2_parameters.set("algorithm", "parma");
    break;
  case 16:
    zoltan2_parameters.set("algorithm", "scotch");
    break;
  case 17:
    zoltan2_parameters.set("algorithm", "ptscotch");
    break;
  case 18:
    zoltan2_parameters.set("algorithm", "block");
    break;
  case 19:
    zoltan2_parameters.set("algorithm", "cyclic");
    break;
  case 20:
    zoltan2_parameters.set("algorithm", "random");
    break;
  case 21:
    zoltan2_parameters.set("algorithm", "zoltan");
    break;
  case 22:
    zoltan2_parameters.set("algorithm", "nd");
    break;
  default :
    zoltan2_parameters.set("algorithm", "rcb");

  }
  zoltan2_parameters.set("mj_keep_part_boxes", false);


  RCP<xcrsGraph_problem_t> homme_partition_problem (new xcrsGraph_problem_t(ia.getRawPtr(),&zoltan2_parameters,tcomm));

  homme_partition_problem->solve();
  tcomm->barrier();

  std::vector<int> tmp_result_parts(numGlobalCoords, 0);
  int *parts =  (int *)homme_partition_problem->getSolution().getPartListView();

  int fortran_shift = 1;
  if (*mapmethod > 1){
    fortran_shift = 0;
  }

  for (zlno_t lclRow = 0; lclRow < numMyElements; ++lclRow) {
    const zgno_t gblRow = map->getGlobalElement (lclRow);
    tmp_result_parts[gblRow] = parts[lclRow] + fortran_shift;
  }

  Teuchos::reduceAll<int, int>(
      *(tcomm),
      Teuchos::REDUCE_SUM,
      numGlobalCoords,
      &(tmp_result_parts[0]),
      result_parts);
}


void zoltan_map_problem(
    int *nelem,
    int *xadj,
    int *adjncy,
    double *adjwgt,
    double *vwgt,
    int *nparts,
    MPI_Comm comm,
    double *xcoord,
    double *ycoord,
    double *zcoord, int *coord_dimension,
    int *result_parts,
    int *partmethod,
    int *mapmethod){

  using namespace Teuchos;

  typedef int zlno_t;
  typedef int zgno_t;
  typedef int part_t;
  typedef double zscalar_t;



  typedef Tpetra::Map<>::node_type znode_t;
  typedef Tpetra::Map<zlno_t, zgno_t, znode_t> map_t;
  size_t numGlobalCoords = *nelem;


  Teuchos::RCP<const Teuchos::Comm<int> > tcomm =  Teuchos::createSerialComm<int>();
  Teuchos::RCP<const Teuchos::Comm<int> > global_comm =
      Teuchos::RCP<const Teuchos::Comm<int> > (new Teuchos::MpiComm<int> (comm));


  Teuchos::ParameterList problemParams;
  Teuchos::RCP<Zoltan2::Environment> env (new Zoltan2::Environment(problemParams, global_comm));
  RCP<Zoltan2::TimerManager> timer(new Zoltan2::TimerManager(global_comm, &std::cout, Zoltan2::MACRO_TIMERS));
  env->setTimer(timer);


  env->timerStart(Zoltan2::MACRO_TIMERS, "TpetraGraphCreate");
  RCP<const map_t> map = rcp (new map_t (numGlobalCoords, 0, tcomm));

  typedef Tpetra::CrsGraph<zlno_t, zgno_t, znode_t> tcrsGraph_t;
  RCP<tcrsGraph_t> TpetraCrsGraph(new tcrsGraph_t (map, 0));

  const zlno_t numMyElements = map->getNodeNumElements ();
  const zgno_t myBegin = map->getGlobalElement (0);

  for (zlno_t lclRow = 0; lclRow < numMyElements; ++lclRow) {
    const zgno_t gblRow = map->getGlobalElement (lclRow);
    zgno_t begin = xadj[gblRow];
    zgno_t end = xadj[gblRow + 1];
    const ArrayView< const zgno_t > indices(adjncy+begin, end-begin);
    TpetraCrsGraph->insertGlobalIndices(gblRow, indices);
  }
  TpetraCrsGraph->fillComplete ();
  env->timerStop(Zoltan2::MACRO_TIMERS, "TpetraGraphCreate");

  env->timerStart(Zoltan2::MACRO_TIMERS, "AdapterCreate");
  RCP<const tcrsGraph_t> const_data = rcp_const_cast<const tcrsGraph_t>(TpetraCrsGraph);
  typedef Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> tMVector_t;
  typedef Zoltan2::XpetraCrsGraphAdapter<tcrsGraph_t, tMVector_t> adapter_t;
  RCP<adapter_t> ia (new adapter_t(const_data/*,(int)vtx_weights.size(),(int)edge_weights.size()*/, 1, 1));

  /***********************************SET COORDINATES*********************/
  const int coord_dim = *coord_dimension;
  //const int coord_dim = 2;
  // make an array of array views containing the coordinate data
  Teuchos::Array<Teuchos::ArrayView<const zscalar_t> > coordView(coord_dim);

  if(numMyElements > 0){
    Teuchos::ArrayView<const zscalar_t> a(xcoord + myBegin, numMyElements);
    coordView[0] = a;
    Teuchos::ArrayView<const zscalar_t> b(ycoord + myBegin, numMyElements);
    coordView[1] = b;
    if (coord_dim == 3){
      Teuchos::ArrayView<const zscalar_t> c(zcoord + myBegin, numMyElements);
      coordView[2] = c;
    }
  }
  else {
    Teuchos::ArrayView<const zscalar_t> a;
    coordView[0] = a;
    coordView[1] = a;
    if (coord_dim == 3){
      coordView[2] = a;
    }
  }

  RCP<tMVector_t> coords(new tMVector_t(map, coordView.view(0, coord_dim), coord_dim));//= set multivector;
  RCP<const tMVector_t> const_coords = rcp_const_cast<const tMVector_t>(coords);
  Zoltan2::XpetraMultiVectorAdapter<tMVector_t> *adapter = (new Zoltan2::XpetraMultiVectorAdapter<tMVector_t>(const_coords));

  ia->setCoordinateInput(adapter);
  ia->setEdgeWeights(adjwgt, 1, 0);
  ia->setVertexWeights(vwgt, 1, 0);

  env->timerStop(Zoltan2::MACRO_TIMERS, "AdapterCreate");
  /***********************************SET COORDINATES*********************/

  //int *parts = result_parts;





  zgno_t num_map_task = global_comm->getSize();
  //partitioning performed using sf curve, or metis and so on.
  //this was handled in zoltan part, therefore, for mapping we need to shift the part numbers by 1.
  if (*partmethod <= 4){
    for (int i = 0; i <  *nelem; ++i){
      result_parts[i] = result_parts[i] - 1;
    }
  }
  if (*partmethod == 5){
    //in this case there is no partitioning performed before mapping.
    num_map_task = *nelem;
    for (int i = 0; i <  *nelem; ++i){
      result_parts[i] = i;
    }
  }
  //otherwise part method is > 5 then part numbers are the output of zoltan2 partitioning mehtods.
  //we do not need to shift them.


  if (*mapmethod == 3){
    //the optimization parameters that are specific to architecture is set here.

    //this is for BGQ
    problemParams.set("Machine_Optimization_Level", 10);
  }

  env->timerStart(Zoltan2::MACRO_TIMERS, "MachineCreate");

  Zoltan2::MachineRepresentation<zscalar_t, part_t> mach(*global_comm, problemParams);

  env->timerStop(Zoltan2::MACRO_TIMERS, "MachineCreate");

  env->timerStart(Zoltan2::MACRO_TIMERS, "CoordinateTaskMapper");

  Zoltan2::CoordinateTaskMapper<adapter_t, part_t> ctm (
      global_comm,
      Teuchos::rcpFromRef(mach),
      ia,
      num_map_task,
      result_parts,
      env,
      false);
  env->timerStop(Zoltan2::MACRO_TIMERS, "CoordinateTaskMapper");
  //timer->printAndResetToZero();

  //then shift the results for fortran base.
  const int fortran_shift = 1;
  for (zlno_t lclRow = 0; lclRow < numMyElements; ++lclRow) {
    int proc = ctm.getAssignedProcForTask(result_parts[lclRow]);
    result_parts[lclRow] = proc + fortran_shift;
  }
}

void zoltan_mapping_problem(
    int *nelem,
    int *xadj,
    int *adjncy,
    double *adjwgt,
    double *vwgt,
    int *nparts,
    MPI_Comm comm,
    double *xcoord,
    double *ycoord,
    double *zcoord, int *coord_dimension,
    int *result_parts,
    int *partmethod,
    int *mapmethod){

  using namespace Teuchos;

  typedef int zlno_t;
  typedef int zgno_t;
  //typedef int part_t;
  typedef double zscalar_t;



  typedef Tpetra::Map<>::node_type znode_t;
  typedef Tpetra::Map<zlno_t, zgno_t, znode_t> map_t;
  size_t numGlobalCoords = *nelem;


  Teuchos::RCP<const Teuchos::Comm<int> > tcomm =  Teuchos::createSerialComm<int>();
  Teuchos::RCP<const Teuchos::Comm<int> > global_comm =
      Teuchos::RCP<const Teuchos::Comm<int> > (new Teuchos::MpiComm<int> (comm));


  Teuchos::ParameterList problemParams;
  Teuchos::RCP<Zoltan2::Environment> env (new Zoltan2::Environment(problemParams, global_comm));
  RCP<Zoltan2::TimerManager> timer(new Zoltan2::TimerManager(global_comm, &std::cout, Zoltan2::MACRO_TIMERS));
  env->setTimer(timer);


  env->timerStart(Zoltan2::MACRO_TIMERS, "TpetraGraphCreate");
  RCP<const map_t> map = rcp (new map_t (numGlobalCoords, 0, tcomm));

  typedef Tpetra::CrsGraph<zlno_t, zgno_t, znode_t> tcrsGraph_t;
  RCP<tcrsGraph_t> TpetraCrsGraph(new tcrsGraph_t (map, 0));

  const zlno_t numMyElements = map->getNodeNumElements ();
  const zgno_t myBegin = map->getGlobalElement (0);

  for (zlno_t lclRow = 0; lclRow < numMyElements; ++lclRow) {
    const zgno_t gblRow = map->getGlobalElement (lclRow);
    zgno_t begin = xadj[gblRow];
    zgno_t end = xadj[gblRow + 1];
    const ArrayView< const zgno_t > indices(adjncy+begin, end-begin);
    TpetraCrsGraph->insertGlobalIndices(gblRow, indices);
  }
  TpetraCrsGraph->fillComplete ();
  env->timerStop(Zoltan2::MACRO_TIMERS, "TpetraGraphCreate");

  env->timerStart(Zoltan2::MACRO_TIMERS, "AdapterCreate");
  RCP<const tcrsGraph_t> const_data = rcp_const_cast<const tcrsGraph_t>(TpetraCrsGraph);
  typedef Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> tMVector_t;
  typedef Zoltan2::XpetraCrsGraphAdapter<tcrsGraph_t, tMVector_t> adapter_t;
  RCP<adapter_t> ia (new adapter_t(const_data/*,(int)vtx_weights.size(),(int)edge_weights.size()*/, 1, 1));

  /***********************************SET COORDINATES*********************/
  const int coord_dim = *coord_dimension;
  //const int coord_dim = 2;
  // make an array of array views containing the coordinate data
  Teuchos::Array<Teuchos::ArrayView<const zscalar_t> > coordView(coord_dim);

  if(numMyElements > 0){
    Teuchos::ArrayView<const zscalar_t> a(xcoord + myBegin, numMyElements);
    coordView[0] = a;
    Teuchos::ArrayView<const zscalar_t> b(ycoord + myBegin, numMyElements);
    coordView[1] = b;
    if (coord_dim == 3){
      Teuchos::ArrayView<const zscalar_t> c(zcoord + myBegin, numMyElements);
      coordView[2] = c;
    }
  }
  else {
    Teuchos::ArrayView<const zscalar_t> a;
    coordView[0] = a;
    coordView[1] = a;
    if (coord_dim == 3){
      coordView[2] = a;
    }
  }

  RCP<tMVector_t> coords(new tMVector_t(map, coordView.view(0, coord_dim), coord_dim));//= set multivector;
  RCP<const tMVector_t> const_coords = rcp_const_cast<const tMVector_t>(coords);
  Zoltan2::XpetraMultiVectorAdapter<tMVector_t> *adapter = (new Zoltan2::XpetraMultiVectorAdapter<tMVector_t>(const_coords));

  ia->setCoordinateInput(adapter);
  ia->setEdgeWeights(adjwgt, 1, 0);
  ia->setVertexWeights(vwgt, 1, 0);

  env->timerStop(Zoltan2::MACRO_TIMERS, "AdapterCreate");
  /***********************************SET COORDINATES*********************/

  //int *parts = result_parts;


  //partitioning performed using sf curve, or metis and so on.
  //this was not handled in zoltan part, therefore, for mapping we need to shift the part numbers by 1.
  if (*partmethod <= 4){
    for (int i = 0; i <  *nelem; ++i){
      result_parts[i] = result_parts[i] - 1;
    }
  }
  if (*partmethod == 5){
    //in this case there is no partitioning performed before mapping.
    for (int i = 0; i <  *nelem; ++i){
      result_parts[i] = i;
    }
  }

  ArrayRCP<int> initial_part_ids(result_parts, 0,*nelem, false);
  Zoltan2::PartitioningSolution<adapter_t> partitioning_solution(env, global_comm, 0);
  partitioning_solution.setParts(initial_part_ids);

  //otherwise part method is > 5 then part numbers are the output of zoltan2 partitioning mehtods.
  //we do not need to shift them.


  if (*mapmethod >= 2){
    //the optimization parameters that are specific to architecture is set here.
    //this is for BGQ
    problemParams.set("Machine_Optimization_Level", *mapmethod - 2);
    //problemParams.set("machine_coord_transformation", "EIGNORE");
  }

  problemParams.set("mapping_algorithm", "geometric");
  problemParams.set("distributed_input_adapter", false);

#ifdef VISUALIZE_MAPPING
  problemParams.set("Input_RCA_Machine_Coords", "titan_alloc.txt");
#endif
  char* dpf=NULL;
  dpf = getenv ("DIVIDEPRIMEFIRST");
  if (dpf != NULL){
    int val = atoi(dpf);
    if (val) {
      problemParams.set("divide_prime_first", true);
    }
    else {
      problemParams.set("divide_prime_first", false);
    }
  }
  else {
    problemParams.set("divide_prime_first", true);
  }

  dpf = getenv ("REDUCEBESTMAPPING");
  if (dpf != NULL){
    int val = atoi(dpf);
    if (val) {
      problemParams.set("reduce_best_mapping", true);
    }
    else {
      problemParams.set("reduce_best_mapping", false);
    }
  }
  env->timerStart(Zoltan2::MACRO_TIMERS, "MappingProblemCreate");
  Zoltan2::MappingProblem<adapter_t> serial_map_problem(
      ia.getRawPtr(),
      &problemParams,
      global_comm,
      &partitioning_solution);



  env->timerStop(Zoltan2::MACRO_TIMERS, "MappingProblemCreate");
  serial_map_problem.solve(true);
  Zoltan2::MappingSolution<adapter_t> *msoln3 = serial_map_problem.getSolution();
  int *parts =  (int *)msoln3->getPartListView();

  //timer->printAndResetToZero();

  //then shift the results for fortran base.
  const int fortran_shift = 1;
  for (zlno_t lclRow = 0; lclRow < numMyElements; ++lclRow) {
    int proc = parts[lclRow];
    result_parts[lclRow] = proc + fortran_shift;
  }
/*
  typedef Zoltan2::EvaluateMapping<adapter_t> quality_t;
  Teuchos::RCP<quality_t> metricObject_1 = rcp(new quality_t(ia.getRawPtr(),&problemParams,tcomm,msoln3,
      serial_map_problem.getMachine().getRawPtr()));
  if (global_comm->getRank() == 0){
    metricObject_1->printMetrics(std::cout);
  }
*/
}


void zoltan2_print_metrics(
    int *nelem,
    int *xadj,
    int *adjncy,
    double *adjwgt,
    double *vwgt,
    int *nparts,
    MPI_Comm comm,
    int *result_parts){

  sort_graph(nelem,xadj,adjncy,adjwgt,vwgt);



  int nv = *nelem;
  int np = *nparts;
  std::vector <double> part_edge_cuts(nv, 0);
  std::vector <double> part_vertex_weights(np, 0);
  int num_messages = 0;
  double myEdgeCut = 0;
  double weighted_hops = 0;

  Teuchos::RCP<const Teuchos::Comm<int> > tcomm =
      Teuchos::RCP<const Teuchos::Comm<int> > (new Teuchos::MpiComm<int> (comm));
  int myRank = tcomm->getRank();

  Zoltan2::MachineRepresentation<double, int> mach(*tcomm);
  double **proc_coords;
  int mach_coord_dim = mach.getMachineDim();
  int *machine_extent = new int [mach_coord_dim];
  bool *machine_extent_wrap_around = new bool[mach_coord_dim];

  mach.getAllMachineCoordinatesView(proc_coords);
  mach.getMachineExtent(machine_extent);
  mach.getMachineExtentWrapArounds(machine_extent_wrap_around);

  for (int i = 0; i < nv; ++i){
    int part_of_i = result_parts[i] - 1;
    part_vertex_weights[part_of_i] += vwgt[i];
    if (part_of_i != myRank){
      continue;
    }
    const int adj_begin = xadj[i];
    const int adj_end = xadj[i + 1];
    for (int j = adj_begin; j < adj_end; ++j){
      int neighbor_vertex = adjncy[j];
      double neighbor_conn = adjwgt[j];
      int neighbor_part = result_parts[neighbor_vertex] - 1;
      if (neighbor_part != myRank){
        if (part_edge_cuts[neighbor_part] < 0.00001){
          num_messages += 1;
        }
        part_edge_cuts[neighbor_part] += neighbor_conn;
        myEdgeCut += neighbor_conn;
        double hops = 0;
        mach.getHopCount(part_of_i, neighbor_part, hops);
        weighted_hops += hops * neighbor_conn;
      }
    }
  }

  int global_num_messages = 0;

  Teuchos::reduceAll<int, int>(
      *(tcomm),
      Teuchos::REDUCE_SUM,
      1,
      &(num_messages),
      &(global_num_messages));

  int global_max_messages = 0;
  Teuchos::reduceAll<int, int>(
      *(tcomm),
      Teuchos::REDUCE_MAX,
      1,
      &(num_messages),
      &(global_max_messages));

  double global_edge_cut = 0;
  Teuchos::reduceAll<int, double>(
      *(tcomm),
      Teuchos::REDUCE_SUM,
      1,
      &(myEdgeCut),
      &(global_edge_cut));

  double global_max_edge_cut = 0;
  Teuchos::reduceAll<int, double>(
      *(tcomm),
      Teuchos::REDUCE_MAX,
      1,
      &(myEdgeCut),
      &(global_max_edge_cut));

  double total_weighted_hops = 0;
  Teuchos::reduceAll<int, double>(
      *(tcomm),
      Teuchos::REDUCE_SUM,
      1,
      &(weighted_hops),
      &(total_weighted_hops));

  if (myRank == 0){
    std::cout << "\tGLOBAL NUM MESSAGES:" << global_num_messages << std::endl
              << "\tMAX MESSAGES:       " << global_max_messages << std::endl
              << "\tGLOBAL EDGE CUT:    " << global_edge_cut << std::endl
              << "\tMAX EDGE CUT:       " << global_max_edge_cut << std::endl
              << "\tTOTAL WEIGHTED HOPS:" << total_weighted_hops << std::endl;

  }
  delete [] machine_extent_wrap_around;
  delete [] machine_extent;


}


void zoltan2_print_metrics2(
    int *nelem,
    int *xadj,
    int *adjncy,
    double *adjwgt,
    double *vwgt,
    int *nparts,
    MPI_Comm comm,
    int *result_parts){


  using namespace Teuchos;

  typedef int zlno_t;
  typedef int zgno_t;
  //typedef int part_t;
  typedef double zscalar_t;



  typedef Tpetra::Map<>::node_type znode_t;
  typedef Tpetra::Map<zlno_t, zgno_t, znode_t> map_t;
  size_t numGlobalCoords = *nelem;


  Teuchos::RCP<const Teuchos::Comm<int> > tcomm =  Teuchos::createSerialComm<int>();
  Teuchos::RCP<const Teuchos::Comm<int> > global_comm =
      Teuchos::RCP<const Teuchos::Comm<int> > (new Teuchos::MpiComm<int> (comm));


  Teuchos::ParameterList problemParams;
  Teuchos::RCP<Zoltan2::Environment> env (new Zoltan2::Environment(problemParams, global_comm));
  RCP<Zoltan2::TimerManager> timer(new Zoltan2::TimerManager(global_comm, &std::cout, Zoltan2::MACRO_TIMERS));
  env->setTimer(timer);


  env->timerStart(Zoltan2::MACRO_TIMERS, "TpetraGraphCreate");
  RCP<const map_t> map = rcp (new map_t (numGlobalCoords, 0, tcomm));

  typedef Tpetra::CrsGraph<zlno_t, zgno_t, znode_t> tcrsGraph_t;
  RCP<tcrsGraph_t> TpetraCrsGraph(new tcrsGraph_t (map, 0));

  const zlno_t numMyElements = map->getNodeNumElements ();
  //const zgno_t myBegin = map->getGlobalElement (0);

  for (zlno_t lclRow = 0; lclRow < numMyElements; ++lclRow) {
    const zgno_t gblRow = map->getGlobalElement (lclRow);
    zgno_t begin = xadj[gblRow];
    zgno_t end = xadj[gblRow + 1];
    const ArrayView< const zgno_t > indices(adjncy+begin, end-begin);
    TpetraCrsGraph->insertGlobalIndices(gblRow, indices);
  }
  TpetraCrsGraph->fillComplete ();
  env->timerStop(Zoltan2::MACRO_TIMERS, "TpetraGraphCreate");

  env->timerStart(Zoltan2::MACRO_TIMERS, "AdapterCreate");
  RCP<const tcrsGraph_t> const_data = rcp_const_cast<const tcrsGraph_t>(TpetraCrsGraph);
  typedef Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> tMVector_t;
  typedef Zoltan2::XpetraCrsGraphAdapter<tcrsGraph_t, tMVector_t> adapter_t;
  RCP<adapter_t> ia (new adapter_t(const_data/*,(int)vtx_weights.size(),(int)edge_weights.size()*/, 1, 1));

  /***********************************SET COORDINATES*********************/
  ia->setEdgeWeights(adjwgt, 1, 0);
  ia->setVertexWeights(vwgt, 1, 0);
  env->timerStop(Zoltan2::MACRO_TIMERS, "AdapterCreate");
  /***********************************SET COORDINATES*********************/

  Zoltan2::MappingSolution<adapter_t> single_phase_mapping_solution(env, global_comm);
  ArrayRCP<int> initial_part_ids(numMyElements);

  for (zgno_t i = 0; i < numMyElements; ++i){
    initial_part_ids[i] = result_parts[i] - 1;
  }
  single_phase_mapping_solution.setParts(initial_part_ids);
  typedef Zoltan2::EvaluateMapping<adapter_t> quality_t;
  Teuchos::ParameterList distributed_problemParams;

#ifdef VISUALIZE_MAPPING
  Teuchos::ParameterList serial_problemParams_2;
  serial_problemParams_2.set("Input_RCA_Machine_Coords", "titan_alloc.txt");
  Zoltan2::MachineRepresentation<double, int> mach(*global_comm, serial_problemParams_2);
#else
  Zoltan2::MachineRepresentation<double, int> mach(*global_comm);
#endif

  RCP<quality_t> metricObject_1 = rcp(new quality_t(ia.getRawPtr(),&distributed_problemParams,tcomm,
      &single_phase_mapping_solution, &mach));
  if (global_comm->getRank() == 0){
    metricObject_1->printMetrics(std::cout);

#ifdef VISUALIZE_MAPPING
    double ** coords;
    mach.getAllMachineCoordinatesView(coords);

    Zoltan2::visualize_mapping<double, int> (
        0, mach.getMachineDim(), mach.getNumRanks(), coords,
        *nelem, xadj,adjncy, result_parts);
#endif
  }
}

#else //HAVE_TRILINOS
void zoltan_partition_problem(
    int *nelem,
    int *xadj,
    int *adjncy,
    int *adjwgt,
    int *vwgt,
    int *nparts,
    MPI_Comm comm,
    double *xcoord,
    double *ycoord,
    double *zcoord,
    int *result_parts,
    int *mapmethod){
#if HAVE_TRILINOS
  std::cerr << "Trilinos is not compiled with Zoltan2!!" << std::endl;
#else
  std::cerr << "Homme is not compiled with Trilinos!!" << std::endl;
#endif
}
void zoltan_map_problem(
    int *nelem,
    int *xadj,
    int *adjncy,
    int *adjwgt,
    int *vwgt,
    int *nparts,
    MPI_Comm comm,
    double *xcoord,
    double *ycoord,
    double *zcoord,
    int *result_parts,
    int *mapmethod){
#if HAVE_TRILINOS
  std::cerr << "Trilinos is not compiled with Zoltan2!!" << std::endl;
#else
  std::cerr << "Homme is not compiled with Trilinos!!" << std::endl;
#endif}
void zoltan2_print_metrics(
    int *nelem,
    int *xadj,
    int *adjncy,
    double *adjwgt,
    double *vwgt,
    int *nparts,
    MPI_Comm comm,
    int *result_parts){
#if HAVE_TRILINOS
  std::cerr << "Trilinos is not compiled with Zoltan2!!" << std::endl;
#else
  std::cerr << "Homme is not compiled with Trilinos!!" << std::endl;
#endif}

void sort_graph(
    int *nelem,
    int *xadj,
    int *adjncy,
    double *adjwgt,
    double *vwgt){
#if HAVE_TRILINOS
  std::cerr << "Trilinos is not compiled with Zoltan2!!" << std::endl;
#else
  std::cerr << "Homme is not compiled with Trilinos!!" << std::endl;
#endif
  }



#endif // HAVE_TRILINOS

//#endif
