//===-- base/Halo.cpp - halo exchange methods -------------------*- C++ -*-===//
//
// The Halo class and its nested classes Neighbor and ExchList contain all
// the information and methods needed for exchanging the halo elements of
// supported Kokkos array types for a given machine environment (MachEnv)
// and parallel decomposition (Decomp). These exchanges are carried out
// via non-blocking MPI library routines. Constructor and private member
// functions are defined here. The Halo class public member function
// exchangeFullArrayHalo which is called by the user to perform halo
// exchanges on a given array is a template function and thus is defined
// in the associated header file, Halo.h.
//
//===----------------------------------------------------------------------===//

#include "Halo.h"
#include "mpi.h"
#include <algorithm>
#include <iterator>
#include <numeric>

namespace OMEGA {

// create the static class members
Halo *Halo::DefaultHalo = nullptr;
std::map<std::string, std::unique_ptr<Halo>> Halo::AllHalos;

//------------------------------------------------------------------------------
// Local routine that searches a std::vector<I4> for a particular entry and
// returns the index of that entry. If not found, the size is returned
// (corresponding to the last index + 1)

I4 searchVector(const std::vector<I4> &InVector, // vector to search
                I4 Value                         // value to search for
) {

   // first use the std::find routine to determine the iterator location
   auto it = std::find(InVector.begin(), InVector.end(), Value);
   // now translate the iterator into an actual vector index
   I4 LocIndx = std::distance(InVector.begin(), it);

   return LocIndx;

} // end function searchVector (std::vector)

//------------------------------------------------------------------------------
// Construct a new ExchList based on input 2D vector which contains a list
// of indices sorted by halo layer

Halo::ExchList::ExchList(
    const std::vector<std::vector<I4>> List // list of indices
) {

   // First dimension of List is the number of halo layers
   I4 HaloLayers = List.size();

   // Set member vector sizes to number of halo layers
   NList.resize(HaloLayers);
   Offsets.resize(HaloLayers);

   // Copy List into member 2D vector Ind which holds the indices
   Ind = List;

   // Count the total number of elements in the list and set the
   // number of elements in each halo layer
   NTot = 0;
   for (int I = 0; I < HaloLayers; ++I) {
      NList[I] = List[I].size();
      NTot += NList[I];
   }

   // Set the index offsets for each halo layer
   Offsets[0] = 0;
   for (int I = 0; I < HaloLayers - 1; ++I) {
      Offsets[I + 1] = Offsets[I] + NList[I];
   }

} // end ExchList constructor

// Empty constructor for ExchList class
Halo::ExchList::ExchList() = default;

//------------------------------------------------------------------------------
// Construct a new Neighbor given the send and receive lists for each index
// space of the neighboring task identified by NghbrID

Halo::Neighbor::Neighbor(
    const std::vector<std::vector<I4>> &SendCell, // list of cells to send
    const std::vector<std::vector<I4>> &SendEdge, // list of edges to send
    const std::vector<std::vector<I4>> &SendVrtx, // list of vertices to send
    const std::vector<std::vector<I4>> &RecvCell, // list of cells to receive
    const std::vector<std::vector<I4>> &RecvEdge, // list of edges to receive
    const std::vector<std::vector<I4>> &RecvVrtx, // list of vertices to receive
    const I4 NghbrID                              // ID of neighboring task
) {

   // Set the ID of the neighboring task
   TaskID = NghbrID;

   // Construct the member ExchList objects
   SendLists[0] = ExchList(SendCell);
   SendLists[1] = ExchList(SendEdge);
   SendLists[2] = ExchList(SendVrtx);

   RecvLists[0] = ExchList(RecvCell);
   RecvLists[1] = ExchList(RecvEdge);
   RecvLists[2] = ExchList(RecvVrtx);

} // end Neighbor constructor

//------------------------------------------------------------------------------
// Initialize and construct the default Halo. MachEnv and Decomp must already
// be initialized

int Halo::init() {

   I4 IErr{0}; // error code

   MachEnv *DefEnv   = MachEnv::getDefault();
   Decomp *DefDecomp = Decomp::getDefault();

   Halo::DefaultHalo = create("Default", DefEnv, DefDecomp);

   return IErr;

} // End Halo init

// Constructor
//------------------------------------------------------------------------------
// Construct a Halo for the input Name, MachEnv, and Decomp.

Halo::Halo(const std::string &Name, const MachEnv *InEnv,
           const Decomp *InDecomp) {

   I4 IErr{0}; // error code

   // Set pointer for the Decomp
   MyDecomp = InDecomp;

   // Set member variable for the halo width
   HaloWidth = MyDecomp->HaloWidth;

   // Set local task ID and the MPI communicator
   MyTask = InEnv->getMyTask();
   MyComm = InEnv->getComm();

   // Fetch the total number of tasks
   I4 NumTasks = InEnv->getNumTasks();

   // Declare 3D vectors to hold lists of indices generated below which are
   // used to construct a Neighbor for each neighboring task
   std::vector<std::vector<std::vector<I4>>> RecvCellLists;
   std::vector<std::vector<std::vector<I4>>> SendCellLists;
   std::vector<std::vector<std::vector<I4>>> RecvEdgeLists;
   std::vector<std::vector<std::vector<I4>>> SendEdgeLists;
   std::vector<std::vector<std::vector<I4>>> RecvVrtxLists;
   std::vector<std::vector<std::vector<I4>>> SendVrtxLists;

   // Determine which tasks are neighbors to the local task
   IErr = determineNeighbors(NumTasks);
   if (IErr != 0)
      LOG_ERROR("Halo: Error determining neighbors");

   // Generate the exchange lists for each neighboring task in each index space
   IErr = generateExchangeLists(SendCellLists, RecvCellLists, OnCell);
   if (IErr != 0)
      LOG_ERROR("Halo: Error generating exchange lists for Cells");
   IErr = generateExchangeLists(SendEdgeLists, RecvEdgeLists, OnEdge);
   if (IErr != 0)
      LOG_ERROR("Halo: Error generating exchange lists for Edges");
   IErr = generateExchangeLists(SendVrtxLists, RecvVrtxLists, OnVertex);
   if (IErr != 0)
      LOG_ERROR("Halo: Error generating exchange lists for Vertices");

   // Construct the Neighbor objects and save them in class member Neighbors
   for (int INghbr = 0; INghbr < NNghbr; ++INghbr) {
      Neighbors.push_back(Neighbor(SendCellLists[INghbr], SendEdgeLists[INghbr],
                                   SendVrtxLists[INghbr], RecvCellLists[INghbr],
                                   RecvEdgeLists[INghbr], RecvVrtxLists[INghbr],
                                   NeighborList[INghbr]));
   }

} // end Halo constructor

/// Creates a new halo by calling the constructor and puts it in the AllHalos
/// map
Halo *Halo::create(const std::string &Name, const MachEnv *Env,
                   const Decomp *Decomp) {
   // Check to see if a halo of the same name already exists and
   // if so, exit with an error
   if (AllHalos.find(Name) != AllHalos.end()) {
      LOG_ERROR("Attempted to create a Halo with name {} but a Halo of "
                "that name already exists",
                Name);
      return nullptr;
   }

   // create a new halo on the heap and put it in a map of
   // unique_ptrs, which will manage its lifetime
   auto *NewHalo = new Halo(Name, Env, Decomp);
   AllHalos.emplace(Name, NewHalo);

   return NewHalo;
} // end Halo create

// Destructor
//------------------------------------------------------------------------------
// Destroy instance of Halo

Halo::~Halo() {

   // No operations necessary

} // end destructor

//------------------------------------------------------------------------------
// Removes a Halo from AllHalos map and destroys it

void Halo::erase(std::string InName // name of Halo to remove
) {

   AllHalos.erase(InName); // removes the Halo from the map and in
                           // the process, calls the destructor

} // end Halo erase

//------------------------------------------------------------------------------
// Removes all Halos and destroys them

void Halo::clear() {

   AllHalos.clear(); // removes all Halos from the map and in the
                     // process, calls the destructor for each

} // end Halo clear

// Retrieval functions
//------------------------------------------------------------------------------
// Get default Halo

Halo *Halo::getDefault() { return Halo::DefaultHalo; }

//------------------------------------------------------------------------------
// Get Halo by name

Halo *Halo::get(const std::string Name // name of Halo to retrieve
) {

   // look for an instance of this name
   auto it = AllHalos.find(Name);

   // if found, return the Halo pointer
   if (it != AllHalos.end()) {
      return it->second.get();
   } else {
      // otherwise print an error and retrun a null pointer
      LOG_ERROR("Halo::get: Attempt to retrieve non-existent Halo:");
      LOG_ERROR(" {} has not been defined or has been removed", Name);
      return nullptr;
   }
} // end Halo get

//------------------------------------------------------------------------------
// Sets Halo class members NeighborList, NNghbr, SendFlags, and RecvFlags during
// Halo construction

int Halo::determineNeighbors(const I4 NumTasks) {

   I4 IErr{0}; // internal error code
   I4 Err{0};  // error code to return

   NeighborList.clear();

   std::vector<I4> CellTasks;
   std::vector<I4> EdgeTasks;
   std::vector<I4> VertexTasks;

   // create lists of tasks that own elements in the halo of the local task for
   // each index space
   generateListOfTasksInHalo(MyDecomp->NCellsOwned, MyDecomp->NCellsAll,
                             MyDecomp->CellLocH, CellTasks);
   generateListOfTasksInHalo(MyDecomp->NEdgesOwned, MyDecomp->NEdgesAll,
                             MyDecomp->EdgeLocH, EdgeTasks);
   generateListOfTasksInHalo(MyDecomp->NVerticesOwned, MyDecomp->NVerticesAll,
                             MyDecomp->VertexLocH, VertexTasks);

   std::vector<I4> UofCE;
   std::vector<I4> UofCEV;

   // union of the three sets results in a sorted list of IDs for all tasks that
   // own at least one halo element in at least one index space
   std::set_union(CellTasks.begin(), CellTasks.end(), EdgeTasks.begin(),
                  EdgeTasks.end(), std::back_inserter(UofCE));
   std::set_union(UofCE.begin(), UofCE.end(), VertexTasks.begin(),
                  VertexTasks.end(), std::back_inserter(UofCEV));

   std::vector<I4> HaloAll(NumTasks, 0);
   std::vector<I4> OwnedAll(NumTasks, 0);

   // set vector of flags for all tasks in MyComm which signals whether the
   // halo for the local task needs elements from each task
   for (int ITask = 0; ITask < UofCEV.size(); ++ITask) {
      HaloAll[UofCEV[ITask]] = 1;
   }

   // perform all to all with all tasks in MyComm in order to determine if there
   // are any tasks that need elements from the local task that the local task
   // does not already consider a neighbor
   IErr =
       MPI_Alltoall(&HaloAll[0], 1, MPI_INT, &OwnedAll[0], 1, MPI_INT, MyComm);
   if (IErr != 0) {
      LOG_ERROR("Halo: MPI_Alltoall error");
      Err = -1;
   }

   // set vector of IDs for all tasks that need locally owned elements for
   // their halos
   std::vector<I4> AddNeighbors;
   for (int ITask = 0; ITask < NumTasks; ++ITask) {
      if (OwnedAll[ITask])
         AddNeighbors.push_back(ITask);
   }

   // one final union results in a list of IDs for all tasks that need to
   // send elements to the local task or need locally owned elements during
   // a halo exchange in at least one index space, save in Halo member
   // vector NeighborList and set member variable NNghbr
   std::set_union(UofCEV.begin(), UofCEV.end(), AddNeighbors.begin(),
                  AddNeighbors.end(), std::back_inserter(NeighborList));
   NNghbr = NeighborList.size();

   // set SendFlags and RecvFlags for each index space
   setNeighborFlags(CellTasks, OnCell);
   setNeighborFlags(EdgeTasks, OnEdge);
   setNeighborFlags(VertexTasks, OnVertex);

   return Err;
}

//------------------------------------------------------------------------------
// Using input decomposition info for a particular index space, generate a
// sorted list of the tasks that own elements in the Halo of the local task

int Halo::generateListOfTasksInHalo(const I4 NOwned, const I4 NAll,
                                    HostArray2DI4 Loc,
                                    std::vector<I4> &ListOfTasks) {

   I4 Err{0}; // error code to return

   // search through halo elements in input Loc array to find each unique
   // task ID and save in ListOfTasks
   for (int Idx = NOwned; Idx < NAll; ++Idx) {
      I4 Value = Loc(Idx, 0);
      auto It  = std::find(ListOfTasks.begin(), ListOfTasks.end(), Value);
      if (It == ListOfTasks.end())
         ListOfTasks.push_back(Value);
   }

   std::sort(ListOfTasks.begin(), ListOfTasks.end());

   return Err;
}

//------------------------------------------------------------------------------
// For the input index space, set SendFlags and RecvFlags vectors that flag
// which Neighbors in NeighborList the local task needs to send elements to or
// receive elements from during a halo exchange

int Halo::setNeighborFlags(std::vector<I4> ListOfTasks,
                           const MeshElement IdxSpace) {

   I4 Err{0}; // error code to return

   // allocate size of SendFlags and RecvFlags
   SendFlags[IdxSpace].resize(NNghbr);
   RecvFlags[IdxSpace].resize(NNghbr);

   for (int INghbr = 0; INghbr < NNghbr; ++INghbr) {
      auto It = std::find(ListOfTasks.begin(), ListOfTasks.end(),
                          NeighborList[INghbr]);
      if (It != ListOfTasks.end()) {
         RecvFlags[IdxSpace][INghbr] = 1;
      } else {
         RecvFlags[IdxSpace][INghbr] = 0;
      }
   }

   // initialize vectors to track MPI errors for each MPI_Isend and MPI_IRecv
   std::vector<I4> SendErr(NNghbr, 0);
   std::vector<I4> RecvErr(NNghbr, 0);

   // initialize vectors of MPI_Request variables to control non-blocking
   // MPI communications
   std::vector<MPI_Request> RecvReqs(NNghbr);
   std::vector<MPI_Request> SendReqs(NNghbr);

   for (int INghbr = 0; INghbr < NNghbr; ++INghbr) {
      RecvErr[INghbr] = MPI_Irecv(&SendFlags[IdxSpace][INghbr], 1, MPI_INT,
                                  NeighborList[INghbr], MPI_ANY_TAG, MyComm,
                                  &RecvReqs[INghbr]);
      if (RecvErr[INghbr] != 0) {
         LOG_ERROR("MPI error {} on task {} receive from task {}",
                   RecvErr[INghbr], MyTask, NeighborList[INghbr]);
         Err = -1;
      }
   }

   for (int INghbr = 0; INghbr < NNghbr; ++INghbr) {
      SendErr[INghbr] =
          MPI_Isend(&RecvFlags[IdxSpace][INghbr], 1, MPI_INT,
                    NeighborList[INghbr], 0, MyComm, &SendReqs[INghbr]);
      if (SendErr[INghbr] != 0) {
         LOG_ERROR("MPI error {} on task {} send to task {}", SendErr[INghbr],
                   MyTask, NeighborList[INghbr]);
         Err = -1;
      }
   }

   MPI_Waitall(NNghbr, SendReqs.data(), MPI_STATUS_IGNORE);
   MPI_Waitall(NNghbr, RecvReqs.data(), MPI_STATUS_IGNORE);

   return Err;
}

//------------------------------------------------------------------------------
// Generate the lists of indices that are used to construct the ExchList
// objects of the input index space for each Neighbor, and save the lists in
// the input 3D vectors SendLists and RecvLists. The first dimension of these
// vectors represent the task in the order they appear in NeighborList, and the
// remaining 2D vector is passed to the Neighbor constructor for that
// neighboring task.

int Halo::generateExchangeLists(
    std::vector<std::vector<std::vector<I4>>> &SendLists,
    std::vector<std::vector<std::vector<I4>>> &RecvLists,
    const MeshElement IndexSpace) {

   I4 IErr{0}; // error code

   // Pointers to the needed info from the Decomp for the input index space
   const I4 *NOwnedPtr{nullptr};
   const I4 *NAllPtr{nullptr};
   HostArray1DI4 NHaloPtr;
   HostArray2DI4 LocPtr;

   // Fetch the proper info for this index space
   switch (IndexSpace) {
   case OnCell:
      NOwnedPtr = &MyDecomp->NCellsOwned;
      NAllPtr   = &MyDecomp->NCellsAll;
      NHaloPtr  = MyDecomp->NCellsHaloH;
      LocPtr    = MyDecomp->CellLocH;
      NumLayers = HaloWidth;

      break;
   case OnEdge:
      NOwnedPtr = &MyDecomp->NEdgesOwned;
      NAllPtr   = &MyDecomp->NEdgesAll;
      NHaloPtr  = MyDecomp->NEdgesHaloH;
      LocPtr    = MyDecomp->EdgeLocH;
      NumLayers = HaloWidth + 1;

      break;
   case OnVertex:
      NOwnedPtr = &MyDecomp->NVerticesOwned;
      NAllPtr   = &MyDecomp->NVerticesAll;
      NHaloPtr  = MyDecomp->NVerticesHaloH;
      LocPtr    = MyDecomp->VertexLocH;
      NumLayers = HaloWidth + 1;

      break;
   }

   // Save the indices that define the bounds of each halo layer
   std::vector<I4> HaloBnds{*NOwnedPtr};
   for (int ILayer = 0; ILayer < HaloWidth; ++ILayer) {
      HaloBnds.push_back(NHaloPtr(ILayer));
   }
   if (IndexSpace != OnCell)
      HaloBnds.push_back(*NAllPtr);

   // Determine the number of halo elements owned by each neighbor in each
   // layer of the halo.
   std::vector<std::vector<I4>> NumNghbrHalo(NNghbr,
                                             std::vector<I4>(NumLayers, 0));

   for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
      for (int Idx = HaloBnds[ILayer]; Idx < HaloBnds[ILayer + 1]; ++Idx) {
         I4 NewVal = LocPtr(Idx, 0);
         I4 INghbr = searchVector(NeighborList, NewVal);
         ++NumNghbrHalo[INghbr][ILayer];
      }
   }

   // Allocate vectors to save halo index info extracted from Decomp now that
   // the number of halo elements owned by each neighboring task is known.
   std::vector<std::vector<I4>> HaloIdx;

   RecvLists.resize(NNghbr);
   HaloIdx.resize(NNghbr);

   for (int INghbr = 0; INghbr < NNghbr; ++INghbr) {
      HaloIdx[INghbr].resize(NumLayers);
      RecvLists[INghbr].resize(NumLayers);
      for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
         RecvLists[INghbr][ILayer].resize(NumNghbrHalo[INghbr][ILayer]);
      }
      I4 TotHalo = std::accumulate(NumNghbrHalo[INghbr].begin(),
                                   NumNghbrHalo[INghbr].end(), 0);
      HaloIdx[INghbr].resize(TotHalo);
   }

   // Save the needed halo info from Decomp in RecvLists and HaloIdx. RecvLists
   // will contain the local indices of halo elements that are owned by
   // neighbors, sorted by neighbor and halo layer, these are then ready to
   // construct the ExchList objects for each Neighbor for this index space.
   // HaloIdx will collect lists of indices for the same mesh elements as
   // defined on the neighboring tasks that own them, this array is collapsed
   // along the halo layer dimension to facilitate MPI communication since each
   // of these lists need to be sent to the corresponding neighboring task.
   std::vector<I4> IOffsets(NNghbr, 0);
   for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
      std::vector<I4> ListIdx(NNghbr, 0);
      for (int Idx = HaloBnds[ILayer]; Idx < HaloBnds[ILayer + 1]; ++Idx) {
         I4 INghbr = searchVector(NeighborList, LocPtr(Idx, 0));
         I4 IList  = ListIdx[INghbr]++;
         RecvLists[INghbr][ILayer][IList]          = Idx;
         HaloIdx[INghbr][IOffsets[INghbr] + IList] = LocPtr(Idx, 1);
      }
      for (int INghbr = 0; INghbr < NNghbr; ++INghbr) {
         IOffsets[INghbr] += NumNghbrHalo[INghbr][ILayer];
      }
   }

   // Allocate a vector of vectors for each neighbor that will receive the
   // number of locally owned elements that belong to each layer of the halos
   // of neighboring tasks.
   std::vector<std::vector<I4>> NumLocalHalo(NNghbr,
                                             std::vector<I4>(NumLayers, 0));

   // Exchange halo sizes with neighboring tasks.
   IErr = exchangeVectorInt(NumNghbrHalo, NumLocalHalo);

   // Now that the number of locally owned elements that belong to the halos
   // of neighboring tasks is known, allocate space to receive this info.
   std::vector<std::vector<I4>> OwnedIdx;

   OwnedIdx.resize(NNghbr);
   for (int INghbr = 0; INghbr < NNghbr; ++INghbr) {
      I4 TotHalo = std::accumulate(NumLocalHalo[INghbr].begin(),
                                   NumLocalHalo[INghbr].end(), 0);
      OwnedIdx[INghbr].resize(TotHalo);
   }

   // Send indices of elements needed by the local halo and receive
   // the indices of elements locally owned that are needed by the
   // halos of neighboring tasks.
   IErr = exchangeVectorInt(HaloIdx, OwnedIdx);

   // Sort out the received lists of indices by halo layer and save
   // in SendLists, these are now ready to construct the ExchList
   // objects for each Neighbor for this index space.
   SendLists.resize(NNghbr);
   for (int INghbr = 0; INghbr < NNghbr; ++INghbr) {
      SendLists[INghbr].resize(NumLayers);
      I4 IOffset = 0;
      for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
         SendLists[INghbr][ILayer].resize(NumLocalHalo[INghbr][ILayer]);
         for (int IList = 0; IList < NumLocalHalo[INghbr][ILayer]; ++IList) {
            SendLists[INghbr][ILayer][IList] =
                OwnedIdx[INghbr][IOffset + IList];
         }
         IOffset += NumLocalHalo[INghbr][ILayer];
      }
   }

   return IErr;
} // end generateExchangeLists

//------------------------------------------------------------------------------
// Exchange 1D integer vectors with each neighbor. Takes as input two 2D
// vectors, SendVec and RecvVec, where the first dimension is the neighboring
// task to send to or receive from, and the second dimension is the 1D vector
// to send to that task or memory space to receive a vector from that task.
// This communication is done using nonblocking MPI routines MPI_Isend
// and MPI_Irecv

int Halo::exchangeVectorInt(
    const std::vector<std::vector<I4>> &SendVec, // vector of vectors to send
    std::vector<std::vector<I4>> &RecvVec // space to receive sent vectors
) {

   // initialize vectors to track MPI errors for each MPI_Isend and MPI_Irecv
   std::vector<I4> SendErr(NNghbr, 0);
   std::vector<I4> RecvErr(NNghbr, 0);

   I4 Err{0}; // error code to return

   // initialize vectors of MPI_Request variables to control non-blocking
   // MPI communications
   std::vector<MPI_Request> RecvReqs(NNghbr);
   std::vector<MPI_Request> SendReqs(NNghbr);

   for (int INghbr = 0; INghbr < NNghbr; ++INghbr) {
      I4 DimLen = RecvVec[INghbr].size();
      RecvErr[INghbr] =
          MPI_Irecv(&RecvVec[INghbr][0], DimLen, MPI_INT, NeighborList[INghbr],
                    MPI_ANY_TAG, MyComm, &RecvReqs[INghbr]);
      if (RecvErr[INghbr] != 0) {
         LOG_ERROR("MPI error {} on task {} receive from task {}",
                   RecvErr[INghbr], MyTask, NeighborList[INghbr]);
         Err = -1;
      }
   }

   for (int INghbr = 0; INghbr < NNghbr; ++INghbr) {
      I4 DimLen = SendVec[INghbr].size();
      SendErr[INghbr] =
          MPI_Isend(&SendVec[INghbr][0], DimLen, MPI_INT, NeighborList[INghbr],
                    0, MyComm, &SendReqs[INghbr]);
      if (SendErr[INghbr] != 0) {
         LOG_ERROR("MPI error {} on task {} send to task {}", SendErr[INghbr],
                   MyTask, NeighborList[INghbr]);
         Err = -1;
      }
   }

   MPI_Waitall(NNghbr, SendReqs.data(), MPI_STATUS_IGNORE);
   MPI_Waitall(NNghbr, RecvReqs.data(), MPI_STATUS_IGNORE);

   return Err;
} // end exchangeVectorInt

//------------------------------------------------------------------------------
// Allocate RecvBuffer and prepare for MPI communication by calling MPI_Irecv
// for each Neighbor

int Halo::startReceives() {

   // Initialize vector to track MPI errors for each MPI_Irecv
   std::vector<I4> IErr(NNghbr, 0);

   I4 Err{0}; // Error code to return

   for (int INghbr = 0; INghbr < NNghbr; ++INghbr) {
      if (RecvFlags[MyElem][INghbr]) {
         MyNeighbor    = &Neighbors[INghbr];
         I4 BufferSize = TotSize * MyNeighbor->RecvLists[MyElem].NTot;
         MyNeighbor->RecvBuffer.resize(BufferSize);
         IErr[INghbr] = MPI_Irecv(&MyNeighbor->RecvBuffer[0], BufferSize,
                                  MPI_DOUBLE, MyNeighbor->TaskID, MPI_ANY_TAG,
                                  MyComm, &MyNeighbor->RReq);
         if (IErr[INghbr] != 0) {
            LOG_ERROR("MPI error {} on task {} receive from task {}",
                      IErr[INghbr], MyTask, MyNeighbor->TaskID);
            Err = -1;
         }
      }
   }

   return Err;
} // end startReceives

//------------------------------------------------------------------------------
// Initiate MPI communication by calling MPI_Isend for each Neighbor to send
// the packed buffers to each task

int Halo::startSends() {

   // Initialize vector to track MPI errors for each MPI_Isend
   std::vector<I4> IErr(NNghbr, 0);

   I4 Err{0}; // Error code to return

   for (int INghbr = 0; INghbr < NNghbr; ++INghbr) {
      if (SendFlags[MyElem][INghbr]) {
         MyNeighbor    = &Neighbors[INghbr];
         I4 BufferSize = TotSize * MyNeighbor->SendLists[MyElem].NTot;
         IErr[INghbr] =
             MPI_Isend(&MyNeighbor->SendBuffer[0], BufferSize, MPI_DOUBLE,
                       MyNeighbor->TaskID, 0, MyComm, &MyNeighbor->SReq);
         if (IErr[INghbr] != 0) {
            LOG_ERROR("MPI error {} on task {} send to task {}", IErr[INghbr],
                      MyTask, MyNeighbor->TaskID);
            Err = -1;
         }
      }
   }

   return Err;
} // end startSends

//------------------------------------------------------------------------------
// The packBuffer function is overloaded to all supported data types. First, the
// send buffer for the neighbor is allocated with enough space to send all the
// halo elements. Then the exchange list for the neighbor and index space
// is used to select the proper elements and pack them into the send buffer.
// In multidimensional arrays the second fastest index (second index from the
// right) is the mesh element dimension. For integer arrays, the value is
// recast as a Real in a bit-preserving manner using reinterpret_cast to pack
// into the buffer, which is of type std::vector<Real>.

int Halo::packBuffer(const HostArray1DI4 Array) {

   ExchList *MyList = &MyNeighbor->SendLists[MyElem];

   MyNeighbor->SendBuffer.resize(MyList->NTot);

   for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
      for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
         I4 IBuff = MyList->Offsets[ILayer] + IExch;
         MyNeighbor->SendBuffer[IBuff] =
             reinterpret_cast<R8 &>(Array(MyList->Ind[ILayer][IExch]));
      }
   }

   return 0;
} // end packBuffer HostArray1DI4

int Halo::packBuffer(const HostArray1DI8 Array) {

   ExchList *MyList = &MyNeighbor->SendLists[MyElem];

   MyNeighbor->SendBuffer.resize(MyList->NTot);

   for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
      for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
         I4 IBuff = MyList->Offsets[ILayer] + IExch;
         MyNeighbor->SendBuffer[IBuff] =
             reinterpret_cast<R8 &>(Array(MyList->Ind[ILayer][IExch]));
      }
   }

   return 0;
} // end packBuffer HostArray1DI8

int Halo::packBuffer(const HostArray1DR4 Array) {

   ExchList *MyList = &MyNeighbor->SendLists[MyElem];

   MyNeighbor->SendBuffer.resize(MyList->NTot);

   for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
      for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
         I4 IBuff                      = MyList->Offsets[ILayer] + IExch;
         MyNeighbor->SendBuffer[IBuff] = Array(MyList->Ind[ILayer][IExch]);
      }
   }

   return 0;
} // end packBuffer HostArray1DR4

int Halo::packBuffer(const HostArray1DR8 Array) {

   ExchList *MyList = &MyNeighbor->SendLists[MyElem];

   MyNeighbor->SendBuffer.resize(MyList->NTot);

   for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
      for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
         I4 IBuff                      = MyList->Offsets[ILayer] + IExch;
         MyNeighbor->SendBuffer[IBuff] = Array(MyList->Ind[ILayer][IExch]);
      }
   }

   return 0;
} // end packBuffer HostArray1DR8

int Halo::packBuffer(const HostArray2DI4 Array) {

   ExchList *MyList = &MyNeighbor->SendLists[MyElem];
   int NJ           = Array.extent(1);

   MyNeighbor->SendBuffer.resize(MyList->NTot * TotSize);

   for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
      for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
         for (int J = 0; J < NJ; ++J) {
            I4 IBuff = (MyList->Offsets[ILayer] + IExch) * NJ + J;
            MyNeighbor->SendBuffer[IBuff] =
                reinterpret_cast<R8 &>(Array(MyList->Ind[ILayer][IExch], J));
         }
      }
   }

   return 0;
} // end packBuffer HostArray2DI4

int Halo::packBuffer(const HostArray2DI8 Array) {

   ExchList *MyList = &MyNeighbor->SendLists[MyElem];
   int NJ           = Array.extent(1);

   MyNeighbor->SendBuffer.resize(MyList->NTot * TotSize);

   for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
      for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
         for (int J = 0; J < NJ; ++J) {
            I4 IBuff = (MyList->Offsets[ILayer] + IExch) * NJ + J;
            MyNeighbor->SendBuffer[IBuff] =
                reinterpret_cast<R8 &>(Array(MyList->Ind[ILayer][IExch], J));
         }
      }
   }

   return 0;
} // end packBuffer HostArray2DI8

int Halo::packBuffer(const HostArray2DR4 Array) {

   ExchList *MyList = &MyNeighbor->SendLists[MyElem];
   int NJ           = Array.extent(1);

   MyNeighbor->SendBuffer.resize(MyList->NTot * TotSize);

   for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
      for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
         for (int J = 0; J < NJ; ++J) {
            I4 IBuff = (MyList->Offsets[ILayer] + IExch) * NJ + J;
            MyNeighbor->SendBuffer[IBuff] =
                Array(MyList->Ind[ILayer][IExch], J);
         }
      }
   }

   return 0;
} // end packBuffer HostArray2DR4

int Halo::packBuffer(const HostArray2DR8 Array) {

   ExchList *MyList = &MyNeighbor->SendLists[MyElem];
   int NJ           = Array.extent(1);

   MyNeighbor->SendBuffer.resize(MyList->NTot * TotSize);

   for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
      for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
         for (int J = 0; J < NJ; ++J) {
            I4 IBuff = (MyList->Offsets[ILayer] + IExch) * NJ + J;
            MyNeighbor->SendBuffer[IBuff] =
                Array(MyList->Ind[ILayer][IExch], J);
         }
      }
   }

   return 0;
} // end packBuffer HostArray2DR8

int Halo::packBuffer(const HostArray3DI4 Array) {

   ExchList *MyList = &MyNeighbor->SendLists[MyElem];
   int NK           = Array.extent(0);
   int NJ           = Array.extent(2);

   MyNeighbor->SendBuffer.resize(MyList->NTot * TotSize);

   for (int K = 0; K < NK; ++K) {
      for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
         for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
            for (int J = 0; J < NJ; ++J) {
               I4 IBuff =
                   (K * MyList->NTot + MyList->Offsets[ILayer] + IExch) * NJ +
                   J;
               MyNeighbor->SendBuffer[IBuff] = reinterpret_cast<R8 &>(
                   Array(K, MyList->Ind[ILayer][IExch], J));
            }
         }
      }
   }

   return 0;
} // end packBuffer HostArray3DI4

int Halo::packBuffer(const HostArray3DI8 Array) {

   ExchList *MyList = &MyNeighbor->SendLists[MyElem];
   int NK           = Array.extent(0);
   int NJ           = Array.extent(2);

   MyNeighbor->SendBuffer.resize(MyList->NTot * TotSize);

   for (int K = 0; K < NK; ++K) {
      for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
         for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
            for (int J = 0; J < NJ; ++J) {
               I4 IBuff =
                   (K * MyList->NTot + MyList->Offsets[ILayer] + IExch) * NJ +
                   J;
               MyNeighbor->SendBuffer[IBuff] = reinterpret_cast<R8 &>(
                   Array(K, MyList->Ind[ILayer][IExch], J));
            }
         }
      }
   }

   return 0;
} // end packBuffer HostArray3DI8

int Halo::packBuffer(const HostArray3DR4 Array) {

   ExchList *MyList = &MyNeighbor->SendLists[MyElem];
   int NK           = Array.extent(0);
   int NJ           = Array.extent(2);

   MyNeighbor->SendBuffer.resize(MyList->NTot * TotSize);

   for (int K = 0; K < NK; ++K) {
      for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
         for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
            for (int J = 0; J < NJ; ++J) {
               I4 IBuff =
                   (K * MyList->NTot + MyList->Offsets[ILayer] + IExch) * NJ +
                   J;
               MyNeighbor->SendBuffer[IBuff] =
                   Array(K, MyList->Ind[ILayer][IExch], J);
            }
         }
      }
   }

   return 0;
} // end packBuffer HostArray3DR4

int Halo::packBuffer(const HostArray3DR8 Array) {

   ExchList *MyList = &MyNeighbor->SendLists[MyElem];
   int NK           = Array.extent(0);
   int NJ           = Array.extent(2);

   MyNeighbor->SendBuffer.resize(MyList->NTot * TotSize);

   for (int K = 0; K < NK; ++K) {
      for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
         for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
            for (int J = 0; J < NJ; ++J) {
               I4 IBuff =
                   (K * MyList->NTot + MyList->Offsets[ILayer] + IExch) * NJ +
                   J;
               MyNeighbor->SendBuffer[IBuff] =
                   Array(K, MyList->Ind[ILayer][IExch], J);
            }
         }
      }
   }

   return 0;
} // end packBuffer HostArray3DR8

int Halo::packBuffer(const HostArray4DI4 Array) {

   ExchList *MyList = &MyNeighbor->SendLists[MyElem];
   int NL           = Array.extent(0);
   int NK           = Array.extent(1);
   int NJ           = Array.extent(3);

   MyNeighbor->SendBuffer.resize(MyList->NTot * TotSize);

   for (int L = 0; L < NL; ++L) {
      for (int K = 0; K < NK; ++K) {
         for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
            for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
               for (int J = 0; J < NJ; ++J) {
                  I4 IBuff = ((L * NK + K) * MyList->NTot +
                              MyList->Offsets[ILayer] + IExch) *
                                 NJ +
                             J;
                  MyNeighbor->SendBuffer[IBuff] = reinterpret_cast<R8 &>(
                      Array(L, K, MyList->Ind[ILayer][IExch], J));
               }
            }
         }
      }
   }

   return 0;
} // end packBuffer HostArray4DI4

int Halo::packBuffer(const HostArray4DI8 Array) {

   ExchList *MyList = &MyNeighbor->SendLists[MyElem];
   int NL           = Array.extent(0);
   int NK           = Array.extent(1);
   int NJ           = Array.extent(3);

   MyNeighbor->SendBuffer.resize(MyList->NTot * TotSize);

   for (int L = 0; L < NL; ++L) {
      for (int K = 0; K < NK; ++K) {
         for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
            for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
               for (int J = 0; J < NJ; ++J) {
                  I4 IBuff = ((L * NK + K) * MyList->NTot +
                              MyList->Offsets[ILayer] + IExch) *
                                 NJ +
                             J;
                  MyNeighbor->SendBuffer[IBuff] = reinterpret_cast<R8 &>(
                      Array(L, K, MyList->Ind[ILayer][IExch], J));
               }
            }
         }
      }
   }

   return 0;
} // end packBuffer HostArray4DI8

int Halo::packBuffer(const HostArray4DR4 Array) {

   ExchList *MyList = &MyNeighbor->SendLists[MyElem];
   int NL           = Array.extent(0);
   int NK           = Array.extent(1);
   int NJ           = Array.extent(3);

   MyNeighbor->SendBuffer.resize(MyList->NTot * TotSize);

   for (int L = 0; L < NL; ++L) {
      for (int K = 0; K < NK; ++K) {
         for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
            for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
               for (int J = 0; J < NJ; ++J) {
                  I4 IBuff = ((L * NK + K) * MyList->NTot +
                              MyList->Offsets[ILayer] + IExch) *
                                 NJ +
                             J;
                  MyNeighbor->SendBuffer[IBuff] =
                      Array(L, K, MyList->Ind[ILayer][IExch], J);
               }
            }
         }
      }
   }

   return 0;
} // end packBuffer HostArray4DR4

int Halo::packBuffer(const HostArray4DR8 Array) {

   ExchList *MyList = &MyNeighbor->SendLists[MyElem];
   int NL           = Array.extent(0);
   int NK           = Array.extent(1);
   int NJ           = Array.extent(3);

   MyNeighbor->SendBuffer.resize(MyList->NTot * TotSize);

   for (int L = 0; L < NL; ++L) {
      for (int K = 0; K < NK; ++K) {
         for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
            for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
               for (int J = 0; J < NJ; ++J) {
                  I4 IBuff = ((L * NK + K) * MyList->NTot +
                              MyList->Offsets[ILayer] + IExch) *
                                 NJ +
                             J;
                  MyNeighbor->SendBuffer[IBuff] =
                      Array(L, K, MyList->Ind[ILayer][IExch], J);
               }
            }
         }
      }
   }

   return 0;
} // end packBuffer HostArray4DR8

int Halo::packBuffer(const HostArray5DI4 Array) {

   ExchList *MyList = &MyNeighbor->SendLists[MyElem];
   int NM           = Array.extent(0);
   int NL           = Array.extent(1);
   int NK           = Array.extent(2);
   int NJ           = Array.extent(4);

   MyNeighbor->SendBuffer.resize(MyList->NTot * TotSize);

   for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
      for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
         for (int J = 0; J < NJ; ++J) {
            for (int K = 0; K < NK; ++K) {
               for (int L = 0; L < NL; ++L) {
                  for (int M = 0; M < NM; ++M) {
                     I4 IBuff = (((M * NL + L) * NK + K) * MyList->NTot +
                                 MyList->Offsets[ILayer] + IExch) *
                                    NJ +
                                J;
                     MyNeighbor->SendBuffer[IBuff] = reinterpret_cast<R8 &>(
                         Array(M, L, K, MyList->Ind[ILayer][IExch], J));
                  }
               }
            }
         }
      }
   }

   return 0;
} // end packBuffer HostArray5DI4

int Halo::packBuffer(const HostArray5DI8 Array) {

   ExchList *MyList = &MyNeighbor->SendLists[MyElem];
   int NM           = Array.extent(0);
   int NL           = Array.extent(1);
   int NK           = Array.extent(2);
   int NJ           = Array.extent(4);

   MyNeighbor->SendBuffer.resize(MyList->NTot * TotSize);

   for (int M = 0; M < NM; ++M) {
      for (int L = 0; L < NL; ++L) {
         for (int K = 0; K < NK; ++K) {
            for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
               for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
                  for (int J = 0; J < NJ; ++J) {
                     I4 IBuff = (((M * NL + L) * NK + K) * MyList->NTot +
                                 MyList->Offsets[ILayer] + IExch) *
                                    NJ +
                                J;
                     MyNeighbor->SendBuffer[IBuff] = reinterpret_cast<R8 &>(
                         Array(M, L, K, MyList->Ind[ILayer][IExch], J));
                  }
               }
            }
         }
      }
   }

   return 0;
} // end packBuffer HostArray5DI8

int Halo::packBuffer(const HostArray5DR4 Array) {

   ExchList *MyList = &MyNeighbor->SendLists[MyElem];
   int NM           = Array.extent(0);
   int NL           = Array.extent(1);
   int NK           = Array.extent(2);
   int NJ           = Array.extent(4);

   MyNeighbor->SendBuffer.resize(MyList->NTot * TotSize);

   for (int M = 0; M < NM; ++M) {
      for (int L = 0; L < NL; ++L) {
         for (int K = 0; K < NK; ++K) {
            for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
               for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
                  for (int J = 0; J < NJ; ++J) {
                     I4 IBuff = (((M * NL + L) * NK + K) * MyList->NTot +
                                 MyList->Offsets[ILayer] + IExch) *
                                    NJ +
                                J;
                     MyNeighbor->SendBuffer[IBuff] =
                         Array(M, L, K, MyList->Ind[ILayer][IExch], J);
                  }
               }
            }
         }
      }
   }

   return 0;
} // end packBuffer HostArray5DR4

int Halo::packBuffer(const HostArray5DR8 Array) {

   ExchList *MyList = &MyNeighbor->SendLists[MyElem];
   int NM           = Array.extent(0);
   int NL           = Array.extent(1);
   int NK           = Array.extent(2);
   int NJ           = Array.extent(4);

   MyNeighbor->SendBuffer.resize(MyList->NTot * TotSize);

   for (int M = 0; M < NM; ++M) {
      for (int L = 0; L < NL; ++L) {
         for (int K = 0; K < NK; ++K) {
            for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
               for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
                  for (int J = 0; J < NJ; ++J) {
                     I4 IBuff = (((M * NL + L) * NK + K) * MyList->NTot +
                                 MyList->Offsets[ILayer] + IExch) *
                                    NJ +
                                J;
                     MyNeighbor->SendBuffer[IBuff] =
                         Array(M, L, K, MyList->Ind[ILayer][IExch], J);
                  }
               }
            }
         }
      }
   }

   return 0;
} // end packBuffer HostArray5DR8

//------------------------------------------------------------------------------
// The unpackBuffer function is overloaded to all supported data types. After
// a message has been received from a neighboring task, the RecvList for the
// corresponding Neighbor and index space is used to save the elements of the
// receive buffer in their proper locations in the input Array. In multi-
// dimensional arrays the second fastest index (second index from the
// right) is the mesh element dimension. For integer arrays, the value from
// the buffer is recast in a bit-preserving manner from a Real to the proper
// integer type (I4 or I8) using reinterpret_cast, and then saved in the
// input Array.

int Halo::unpackBuffer(HostArray1DI4 &Array) {

   ExchList *MyList = &MyNeighbor->RecvLists[MyElem];

   for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
      for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
         I4 IBuff = MyList->Offsets[ILayer] + IExch;
         Array(MyList->Ind[ILayer][IExch]) =
             reinterpret_cast<I4 &>(MyNeighbor->RecvBuffer[IBuff]);
      }
   }

   return 0;
} // end unpackBuffer HostArray1DI4

int Halo::unpackBuffer(HostArray1DI8 &Array) {

   ExchList *MyList = &MyNeighbor->RecvLists[MyElem];

   for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
      for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
         I4 IBuff = MyList->Offsets[ILayer] + IExch;
         Array(MyList->Ind[ILayer][IExch]) =
             reinterpret_cast<I8 &>(MyNeighbor->RecvBuffer[IBuff]);
      }
   }

   return 0;
} // end unpackBuffer HostArray1DI8

int Halo::unpackBuffer(HostArray1DR4 &Array) {

   ExchList *MyList = &MyNeighbor->RecvLists[MyElem];

   for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
      for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
         I4 IBuff                          = MyList->Offsets[ILayer] + IExch;
         Array(MyList->Ind[ILayer][IExch]) = MyNeighbor->RecvBuffer[IBuff];
      }
   }

   return 0;
} // end unpackBuffer HostArray1DR4

int Halo::unpackBuffer(HostArray1DR8 &Array) {

   ExchList *MyList = &MyNeighbor->RecvLists[MyElem];

   for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
      for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
         I4 IBuff                          = MyList->Offsets[ILayer] + IExch;
         Array(MyList->Ind[ILayer][IExch]) = MyNeighbor->RecvBuffer[IBuff];
      }
   }

   return 0;
} // end unpackBuffer HostArray1DR8

int Halo::unpackBuffer(HostArray2DI4 &Array) {

   ExchList *MyList = &MyNeighbor->RecvLists[MyElem];
   int NJ           = Array.extent(1);

   for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
      for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
         for (int J = 0; J < NJ; ++J) {
            I4 IBuff = (MyList->Offsets[ILayer] + IExch) * NJ + J;
            Array(MyList->Ind[ILayer][IExch], J) =
                reinterpret_cast<I4 &>(MyNeighbor->RecvBuffer[IBuff]);
         }
      }
   }

   return 0;
} // end unpackBuffer HostArray2DI4

int Halo::unpackBuffer(HostArray2DI8 &Array) {

   ExchList *MyList = &MyNeighbor->RecvLists[MyElem];
   int NJ           = Array.extent(1);

   for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
      for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
         for (int J = 0; J < NJ; ++J) {
            I4 IBuff = (MyList->Offsets[ILayer] + IExch) * NJ + J;
            Array(MyList->Ind[ILayer][IExch], J) =
                reinterpret_cast<I8 &>(MyNeighbor->RecvBuffer[IBuff]);
         }
      }
   }

   return 0;
} // end unpackBuffer HostArray2DI8

int Halo::unpackBuffer(HostArray2DR4 &Array) {

   ExchList *MyList = &MyNeighbor->RecvLists[MyElem];
   int NJ           = Array.extent(1);

   for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
      for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
         for (int J = 0; J < NJ; ++J) {
            I4 IBuff = (MyList->Offsets[ILayer] + IExch) * NJ + J;
            Array(MyList->Ind[ILayer][IExch], J) =
                MyNeighbor->RecvBuffer[IBuff];
         }
      }
   }

   return 0;
} // end unpackBuffer HostArray2DR4

int Halo::unpackBuffer(HostArray2DR8 &Array) {

   ExchList *MyList = &MyNeighbor->RecvLists[MyElem];
   int NJ           = Array.extent(1);

   for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
      for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
         for (int J = 0; J < NJ; ++J) {
            I4 IBuff = (MyList->Offsets[ILayer] + IExch) * NJ + J;
            Array(MyList->Ind[ILayer][IExch], J) =
                MyNeighbor->RecvBuffer[IBuff];
         }
      }
   }

   return 0;
} // end unpackBuffer HostArray2DR8

int Halo::unpackBuffer(HostArray3DI4 &Array) {

   ExchList *MyList = &MyNeighbor->RecvLists[MyElem];
   int NK           = Array.extent(0);
   int NJ           = Array.extent(2);

   for (int K = 0; K < NK; ++K) {
      for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
         for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
            for (int J = 0; J < NJ; ++J) {
               I4 IBuff =
                   (K * MyList->NTot + MyList->Offsets[ILayer] + IExch) * NJ +
                   J;
               Array(K, MyList->Ind[ILayer][IExch], J) =
                   reinterpret_cast<I4 &>(MyNeighbor->RecvBuffer[IBuff]);
            }
         }
      }
   }

   return 0;
} // end unpackBuffer HostArray3DI4

int Halo::unpackBuffer(HostArray3DI8 &Array) {

   ExchList *MyList = &MyNeighbor->RecvLists[MyElem];
   int NK           = Array.extent(0);
   int NJ           = Array.extent(2);

   for (int K = 0; K < NK; ++K) {
      for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
         for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
            for (int J = 0; J < NJ; ++J) {
               I4 IBuff =
                   (K * MyList->NTot + MyList->Offsets[ILayer] + IExch) * NJ +
                   J;
               Array(K, MyList->Ind[ILayer][IExch], J) =
                   reinterpret_cast<I8 &>(MyNeighbor->RecvBuffer[IBuff]);
            }
         }
      }
   }

   return 0;
} // end unpackBuffer HostArray3DI8

int Halo::unpackBuffer(HostArray3DR4 &Array) {

   ExchList *MyList = &MyNeighbor->RecvLists[MyElem];
   int NK           = Array.extent(0);
   int NJ           = Array.extent(2);

   for (int K = 0; K < NK; ++K) {
      for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
         for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
            for (int J = 0; J < NJ; ++J) {
               I4 IBuff =
                   (K * MyList->NTot + MyList->Offsets[ILayer] + IExch) * NJ +
                   J;
               Array(K, MyList->Ind[ILayer][IExch], J) =
                   MyNeighbor->RecvBuffer[IBuff];
            }
         }
      }
   }

   return 0;
} // end unpackBuffer HostArray3DR4

int Halo::unpackBuffer(HostArray3DR8 &Array) {

   ExchList *MyList = &MyNeighbor->RecvLists[MyElem];
   int NK           = Array.extent(0);
   int NJ           = Array.extent(2);

   for (int K = 0; K < NK; ++K) {
      for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
         for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
            for (int J = 0; J < NJ; ++J) {
               I4 IBuff =
                   (K * MyList->NTot + MyList->Offsets[ILayer] + IExch) * NJ +
                   J;
               Array(K, MyList->Ind[ILayer][IExch], J) =
                   MyNeighbor->RecvBuffer[IBuff];
            }
         }
      }
   }

   return 0;
} // end unpackBuffer HostArray3DR8

int Halo::unpackBuffer(HostArray4DI4 &Array) {

   ExchList *MyList = &MyNeighbor->RecvLists[MyElem];
   int NL           = Array.extent(0);
   int NK           = Array.extent(1);
   int NJ           = Array.extent(3);

   for (int L = 0; L < NL; ++L) {
      for (int K = 0; K < NK; ++K) {
         for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
            for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
               for (int J = 0; J < NJ; ++J) {
                  I4 IBuff = ((L * NK + K) * MyList->NTot +
                              MyList->Offsets[ILayer] + IExch) *
                                 NJ +
                             J;
                  Array(L, K, MyList->Ind[ILayer][IExch], J) =
                      reinterpret_cast<I4 &>(MyNeighbor->RecvBuffer[IBuff]);
               }
            }
         }
      }
   }

   return 0;
} // end unpackBuffer HostArray4DI4

int Halo::unpackBuffer(HostArray4DI8 &Array) {

   ExchList *MyList = &MyNeighbor->RecvLists[MyElem];
   int NL           = Array.extent(0);
   int NK           = Array.extent(1);
   int NJ           = Array.extent(3);

   for (int L = 0; L < NL; ++L) {
      for (int K = 0; K < NK; ++K) {
         for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
            for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
               for (int J = 0; J < NJ; ++J) {
                  I4 IBuff = ((L * NK + K) * MyList->NTot +
                              MyList->Offsets[ILayer] + IExch) *
                                 NJ +
                             J;
                  Array(L, K, MyList->Ind[ILayer][IExch], J) =
                      reinterpret_cast<I8 &>(MyNeighbor->RecvBuffer[IBuff]);
               }
            }
         }
      }
   }

   return 0;
} // end unpackBuffer HostArray4DI8

int Halo::unpackBuffer(HostArray4DR4 &Array) {

   ExchList *MyList = &MyNeighbor->RecvLists[MyElem];
   int NL           = Array.extent(0);
   int NK           = Array.extent(1);
   int NJ           = Array.extent(3);

   for (int L = 0; L < NL; ++L) {
      for (int K = 0; K < NK; ++K) {
         for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
            for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
               for (int J = 0; J < NJ; ++J) {
                  I4 IBuff = ((L * NK + K) * MyList->NTot +
                              MyList->Offsets[ILayer] + IExch) *
                                 NJ +
                             J;
                  Array(L, K, MyList->Ind[ILayer][IExch], J) =
                      MyNeighbor->RecvBuffer[IBuff];
               }
            }
         }
      }
   }

   return 0;
} // end unpackBuffer HostArray4DR4

int Halo::unpackBuffer(HostArray4DR8 &Array) {

   ExchList *MyList = &MyNeighbor->RecvLists[MyElem];
   int NL           = Array.extent(0);
   int NK           = Array.extent(1);
   int NJ           = Array.extent(3);

   for (int L = 0; L < NL; ++L) {
      for (int K = 0; K < NK; ++K) {
         for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
            for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
               for (int J = 0; J < NJ; ++J) {
                  I4 IBuff = ((L * NK + K) * MyList->NTot +
                              MyList->Offsets[ILayer] + IExch) *
                                 NJ +
                             J;
                  Array(L, K, MyList->Ind[ILayer][IExch], J) =
                      MyNeighbor->RecvBuffer[IBuff];
               }
            }
         }
      }
   }

   return 0;
} // end unpackBuffer HostArray4DR8

int Halo::unpackBuffer(HostArray5DI4 &Array) {

   ExchList *MyList = &MyNeighbor->RecvLists[MyElem];
   int NM           = Array.extent(0);
   int NL           = Array.extent(1);
   int NK           = Array.extent(2);
   int NJ           = Array.extent(4);

   for (int M = 0; M < NM; ++M) {
      for (int L = 0; L < NL; ++L) {
         for (int K = 0; K < NK; ++K) {
            for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
               for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
                  for (int J = 0; J < NJ; ++J) {
                     I4 IBuff = (((M * NL + L) * NK + K) * MyList->NTot +
                                 MyList->Offsets[ILayer] + IExch) *
                                    NJ +
                                J;
                     Array(M, L, K, MyList->Ind[ILayer][IExch], J) =
                         reinterpret_cast<I4 &>(MyNeighbor->RecvBuffer[IBuff]);
                  }
               }
            }
         }
      }
   }

   return 0;
} // end unpackBuffer HostArray5DI4

int Halo::unpackBuffer(HostArray5DI8 &Array) {

   ExchList *MyList = &MyNeighbor->RecvLists[MyElem];
   int NM           = Array.extent(0);
   int NL           = Array.extent(1);
   int NK           = Array.extent(2);
   int NJ           = Array.extent(4);

   for (int M = 0; M < NM; ++M) {
      for (int L = 0; L < NL; ++L) {
         for (int K = 0; K < NK; ++K) {
            for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
               for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
                  for (int J = 0; J < NJ; ++J) {
                     I4 IBuff = (((M * NL + L) * NK + K) * MyList->NTot +
                                 MyList->Offsets[ILayer] + IExch) *
                                    NJ +
                                J;
                     Array(M, L, K, MyList->Ind[ILayer][IExch], J) =
                         reinterpret_cast<I8 &>(MyNeighbor->RecvBuffer[IBuff]);
                  }
               }
            }
         }
      }
   }

   return 0;
} // end unpackBuffer HostArray5DI8

int Halo::unpackBuffer(HostArray5DR4 &Array) {

   ExchList *MyList = &MyNeighbor->RecvLists[MyElem];
   int NM           = Array.extent(0);
   int NL           = Array.extent(1);
   int NK           = Array.extent(2);
   int NJ           = Array.extent(4);

   for (int M = 0; M < NM; ++M) {
      for (int L = 0; L < NL; ++L) {
         for (int K = 0; K < NK; ++K) {
            for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
               for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
                  for (int J = 0; J < NJ; ++J) {
                     I4 IBuff = (((M * NL + L) * NK + K) * MyList->NTot +
                                 MyList->Offsets[ILayer] + IExch) *
                                    NJ +
                                J;
                     Array(M, L, K, MyList->Ind[ILayer][IExch], J) =
                         MyNeighbor->RecvBuffer[IBuff];
                  }
               }
            }
         }
      }
   }

   return 0;
} // end unpackBuffer HostArray5DR4

int Halo::unpackBuffer(HostArray5DR8 &Array) {

   ExchList *MyList = &MyNeighbor->RecvLists[MyElem];
   int NM           = Array.extent(0);
   int NL           = Array.extent(1);
   int NK           = Array.extent(2);
   int NJ           = Array.extent(4);

   for (int M = 0; M < NM; ++M) {
      for (int L = 0; L < NL; ++L) {
         for (int K = 0; K < NK; ++K) {
            for (int ILayer = 0; ILayer < NumLayers; ++ILayer) {
               for (int IExch = 0; IExch < MyList->NList[ILayer]; ++IExch) {
                  for (int J = 0; J < NJ; ++J) {
                     I4 IBuff = (((M * NL + L) * NK + K) * MyList->NTot +
                                 MyList->Offsets[ILayer] + IExch) *
                                    NJ +
                                J;
                     Array(M, L, K, MyList->Ind[ILayer][IExch], J) =
                         MyNeighbor->RecvBuffer[IBuff];
                  }
               }
            }
         }
      }
   }

   return 0;
} // end unpackBuffer HostArray5DR8

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
