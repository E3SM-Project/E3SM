
#ifdef HAVE_CONFIG_H
#include "config.h.c"
#endif

#include "zoltan_cppinterface.hpp"

#if HAVE_TRILINOS
#if TRILINOS_HAVE_ZOLTAN2
//#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_Map.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_GlobalMPISession.hpp>

#include <iostream>

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
    double *zcoord,
    int *result_parts,
    int *partmethod){
  using namespace Teuchos;

  typedef int zlno_t;
  typedef int zgno_t;
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
  const int coord_dim = 3;
  // make an array of array views containing the coordinate data
  Teuchos::Array<Teuchos::ArrayView<const zscalar_t> > coordView(coord_dim);

  if(numMyElements > 0){
    Teuchos::ArrayView<const zscalar_t> a(xcoord + myBegin, numMyElements);
    coordView[0] = a;
    Teuchos::ArrayView<const zscalar_t> b(ycoord + myBegin, numMyElements);
    coordView[1] = b;
    Teuchos::ArrayView<const zscalar_t> c(zcoord + myBegin, numMyElements);
    coordView[2] = c;
  }
  else {
    Teuchos::ArrayView<const zscalar_t> a;
    coordView[0] = a;
    coordView[1] = a;
    coordView[2] = a;

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
  zoltan2_parameters.set("compute_metrics", "true");
  zoltan2_parameters.set("imbalance_tolerance", "1.0");
  zoltan2_parameters.set("num_global_parts", tcomm->getSize());
  switch (*partmethod){
  case 5:
    zoltan2_parameters.set("algorithm", "rcb");
    break;

  case 6:
    zoltan2_parameters.set("algorithm", "multijagged");
    break;
  case 7:
    zoltan2_parameters.set("algorithm", "rib");
    break;
  case 8:
    zoltan2_parameters.set("algorithm", "hsfc");
    break;
  case 9:
    zoltan2_parameters.set("algorithm", "patoh");
    break;
  case 10:
    zoltan2_parameters.set("algorithm", "phg");
    break;
  case 11:
    zoltan2_parameters.set("algorithm", "metis");
    break;
  case 12:
    zoltan2_parameters.set("algorithm", "parmetis");
    break;
  case 13:
    zoltan2_parameters.set("algorithm", "parma");
    break;
  case 14:
    zoltan2_parameters.set("algorithm", "scotch");
    break;
  case 15:
    zoltan2_parameters.set("algorithm", "ptscotch");
    break;
  case 16:
    zoltan2_parameters.set("algorithm", "block");
    break;
  case 17:
    zoltan2_parameters.set("algorithm", "cyclic");
    break;
  case 18:
    zoltan2_parameters.set("algorithm", "random");
    break;
  case 19:
    zoltan2_parameters.set("algorithm", "zoltan");
    break;
  case 20:
    zoltan2_parameters.set("algorithm", "nd");
    break;
  default :
    zoltan2_parameters.set("algorithm", "multijagged");

  }

  zoltan2_parameters.set("mj_keep_part_boxes", false);
  zoltan2_parameters.set("mj_recursion_depth", "3");
  zoltan2_parameters.set("remap_parts", true);


  RCP<xcrsGraph_problem_t> homme_partition_problem (new xcrsGraph_problem_t(ia.getRawPtr(),&zoltan2_parameters,tcomm));



  homme_partition_problem->solve();
  tcomm->barrier();

  int *parts =  (int *)homme_partition_problem->getSolution().getPartListView();
  std::vector<int> tmp_result_parts(numGlobalCoords, 0);

  for (zlno_t lclRow = 0; lclRow < numMyElements; ++lclRow) {
    const zgno_t gblRow = map->getGlobalElement (lclRow);
    tmp_result_parts[gblRow] = parts[lclRow] + 1;
  }

  Teuchos::reduceAll<int, int>(
      *(tcomm),
      Teuchos::REDUCE_SUM,
      numGlobalCoords,
      &(tmp_result_parts[0]),
      result_parts);

}

#else //TRILINOS_HAVE_ZOLTAN2
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
    int *result_parts){
  std::cerr << "Trilinos is not compiled with Zoltan2!!" << std::endl;
}
#endif
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
    int *result_parts){
  std::cerr << "Homme is not compiled with Trilinos!!" << std::endl;
}
#endif // HAVE_TRILINOS

//#endif
