#ifndef OMEGA_HALO_H
#define OMEGA_HALO_H
//===-- base/Halo.h - halo exchange definitions -----------------*- C++ -*-===//
//
/// \file
/// \brief Defines classes and methods for performing halo exchanges
///
/// This header defines classes and methods needed for exchanging the halo
/// elements of every supported array type in OMEGA with elements of the array
/// on neighboring mesh partitions. The Halo class constructor creates a Halo
/// object based on a given machine environment and parallel decomposition.
/// The Halo object contains all the info needed for performing a halo
/// exchange The halo exchanges are carried out via non-blocking MPI library
/// routines. The Halo class public member function exchangeFullArrayHalo
/// which is called by the user to perform halo exchanges is a template
/// function and thus is fully defined in this header.
///
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "Decomp.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "mpi.h"
#include <memory>
#include <numeric>
#include <type_traits>
#include <vector>

namespace OMEGA {

// Set the default MPI real data type as single or double precision based on
// the default real data type used by OMEGA.
#ifdef OMEGA_SINGLE_PRECISION
static const MPI_Datatype MPI_RealKind = MPI_FLOAT;
#else
static const MPI_Datatype MPI_RealKind = MPI_DOUBLE;
#endif

/// The MeshElement enum identifies the index space to use for a halo exchange.
enum MeshElement { OnCell, OnEdge, OnVertex };

/// The Halo class contains two nested classes, ExchList and Neighbor classes,
/// defined below. The Halo class holds all the Neighbor objects needed by a
/// task to perform a full halo exchange with each of its neighboring tasks for
/// any array defined on the mesh. The local task ID and the MPI communicator
/// handle are also stored here. NumLayers, CurElem, and TotSize are temporary
/// variables utilized by the current halo exchange which are stored here
/// for easy accessibility by the Halo methods.
class Halo {
 private:
   /// Flag to allow passing device arrays to MPI_Irecv and MPI_Isend in
   /// startReceives and startSends, determined by pre-processing parameter
#ifdef OMEGA_MPI_ON_DEVICE
   const bool ExchOnDev = true;
#else
   const bool ExchOnDev = false;
#endif

   /// The default Halo handles halo exchanges for arrays defined on the mesh
   /// with the default decomposition. A pointer is stored here for easy
   /// retrieval.
   static Halo *DefaultHalo;

   /// All halos are tracked/stored within the class as a map paired with a
   /// name for later retrieval.
   static std::map<std::string, std::unique_ptr<Halo>> AllHalos;

   const Decomp *MyDecomp{nullptr}; /// Pointer to decomposition object

   I4 NNghbr;           /// number of neighboring tasks
   I4 MyTask;           /// local MPI Task ID
   I4 HaloWidth;        /// cell width of halo
   I4 NumLayers;        /// number of halo layers for current exchange
   I4 TotSize;          /// Array size at each mesh element for current exchange
   MPI_Comm MyComm;     /// MPI communicator handle
   MeshElement CurElem; /// index space of current array

   /// Forward Declaration of Neighbor class, defined below
   class Neighbor;

   /// Vector of Neighbor class objects for all neighboring tasks, contains
   /// all the exchange lists and buffer memory sufficient for a full halo
   /// exchange in any index space.
   std::vector<Neighbor> Neighbors;

   /// Sorted list of Task IDs of all neighboring tasks in ascending order.
   /// Another task is considered a neighbor if a mesh element owned by that
   /// task is in the halo of the local task in at least one index space, or if
   /// an owned element is in the halo of that task in at least one index space.
   std::vector<I4> NeighborList;

   /// Flags to control which tasks in NeighborList the local task needs to
   /// send elements to and receive elements from for each index space.
   std::vector<I4> SendFlags[3], RecvFlags[3];

   /// The ExchList class contains the information needed to pack values from
   /// an array into a buffer or unpack values from a buffer into an array
   /// for a particular index space and neighbor
   class ExchList {
    private:
      /// Array containing number of mesh elements in each halo layer
      std::vector<I4> NHalo;
      /// Total number of elements in ExchList, sum of NList
      I4 NTot;
      /// Host and device arrays containing starting indices in the
      /// buffer for each halo layer
      HostArray1DI4 OffsetsH;
      Array1DI4 Offsets;
      /// Host and device arrays containing the list of local indices of array
      /// elements to be packed into the send buffer, or the local indices of
      /// elements unpacked from the receive buffer
      HostArray1DI4 IndexH;
      Array1DI4 Index;

      /// The constructor for the ExchList class takes as input an array of
      /// vectors, each containing a list of indices to be sent or received for
      /// each halo layer for a particular neighbor
      ExchList(const std::vector<std::vector<I4>> List);

      /// Empty constructor for ExchList class
      ExchList();

    public:
      /// Destructor
      ~ExchList() {}

      /// Halo and Neighbor are friend class to allow access to private
      /// members of the class
      friend class Halo;
      friend class Neighbor;
   }; // end class ExchList

   /// The Neighbor class contains all the information and buffer memory
   /// needed to carry out a halo exchange in each index space for one
   /// neighboring task.
   class Neighbor {
    private:
      I4 TaskID; /// ID of neighboring task

      /// Arrays of ExchList objects for sends and recieves for each
      /// index space. 0 = OnCell, 1 = OnEdge, 2 = OnVertex
      ExchList SendLists[3], RecvLists[3];
      /// Buffers for MPI communication on host and device
      HostArray1DR8 SendBufferH, RecvBufferH;
      Array1DR8 SendBuffer, RecvBuffer;
      /// MPI request handles for non-blocking MPI communication
      MPI_Request RReq, SReq;

      // Flags to track whether a message has been received and whether
      // the buffer has been unpacked yet. The MPI_Test routine requires an
      // integer, so Received is stored as an I4.
      I4 Received{false};
      bool Unpacked{false};

      /// Neighbor constructor takes takes as input six unique vectors of
      /// vectors containing the indices of array elements to send or receive
      /// for each type of mesh element to be stored in the member ExchList
      /// objects, as well as the ID of the neighboring task.
      Neighbor(const std::vector<std::vector<I4>> &SendCell,
               const std::vector<std::vector<I4>> &SendEdge,
               const std::vector<std::vector<I4>> &SendVert,
               const std::vector<std::vector<I4>> &RecvCell,
               const std::vector<std::vector<I4>> &RecvEdge,
               const std::vector<std::vector<I4>> &RecvVert, const I4 NghbrID);

    public:
      /// Destructor
      ~Neighbor() {}

      /// Halo is a friend class to allow access to private members
      /// of the class
      friend class Halo;

   }; // end class Neighbor

   // Private methods

   /// Uses info from Decomp to generate a sorted list of tasks that own
   /// elements in the the Halo of the local task for a particular index space.
   /// Utilized only during halo construction
   int generateListOfTasksInHalo(const I4 NOwned, const I4 NAll,
                                 HostArray2DI4 Locs,
                                 std::vector<I4> &ListOfTasks);

   /// Set SendFlags and RecvFlags for the input index space. Utilized only
   /// during halo construction
   int setNeighborFlags(std::vector<I4> NeighborElem,
                        const MeshElement IdxSpace);

   /// Uses info from Decomp to determine all tasks which own elements in the
   /// halo of the local task or need locally owned elements for their halo.
   /// Utilized only during halo construction
   int determineNeighbors(const I4 NumTasks);

   /// Send a vector of integers to each neighboring task and receive a vector
   /// of integers from each neighboring task. The first dimension of each
   /// input 2D vector represents the task in the order they appear in
   /// NeighborList. Utilized only during halo construction
   int exchangeVectorInt(const std::vector<std::vector<I4>> &SendVec,
                         std::vector<std::vector<I4>> &RecvVec);

   /// Generate the lists of indices to send to and receive from each
   /// neighboring task for the input IndexSpace and the Decomp pointed to
   /// by member variable MyDecomp, and save the lists in the input 3D
   /// vectors SendLists and RecvLists. The first dimension of these vectors
   /// represent the task in the order they appear in NeighborList, and the
   /// remaining 2D vector is used in constructing a Neighbor object for
   /// that task. Utilized only during halo construction
   int
   generateExchangeLists(std::vector<std::vector<std::vector<I4>>> &SendLists,
                         std::vector<std::vector<std::vector<I4>>> &RecvLists,
                         const MeshElement IndexSpace);

   /// Allocate the recieve buffers and call MPI_Irecv for each Neighbor.
   /// The input bool UseDevBuffer specifies whether or not the device buffer
   /// will be used in the unpackBuffer functionfor unpacking into the array.
   int startReceives(bool UseDevBuffer);

   /// Call MPI_Isend for each Neighbor to send the packed buffers to the
   /// neighboring tasks. The input bool UseDevBuffer specifies whether or not
   /// the device buffer was packed in the packBuffer function.
   int startSends(bool UseDevBuffer);

   /// Function template that returns a bool that is true if the Array is
   /// on the device, or if the device and host memory spaces are the same
   /// space. Used to determine if buffer packs and unpacks are done within
   /// Kokkos parallelFor kernels.
   template <typename T> bool devBufferPUP(const T &Array) {
      bool OnDev = Impl::findArrayMemLoc<T>() == ArrayMemLoc::Both ||
                   Impl::findArrayMemLoc<T>() == ArrayMemLoc::Device;
      return OnDev;
   }

   /// Construct a new halo labeled Name for the input MachEnv and Decomp
   Halo(const std::string &Name, const MachEnv *InEnv, const Decomp *InDecomp);

   // Forbid copy and move construction
   Halo(const Halo &) = delete;
   Halo(Halo &&)      = delete;

 public:
   // Methods

   /// initialize default Halo
   static int init();

   /// Creates a new halo by calling the constructor and puts it in the AllHalos
   /// map
   static Halo *create(const std::string &Name, const MachEnv *Env,
                       const Decomp *Decomp);

   /// Destructor
   ~Halo();

   /// Erase - removes Halo by name
   static void erase(std::string InName ///< [in] name of Halo to remove
   );

   /// Clear - removes all defined Halo instances
   static void clear();

   /// Retrieves a pointer to the default Halo object.
   static Halo *getDefault();

   /// Retrieves a pointer to a Halo object by Name
   static Halo *get(std::string Name);

   /// Buffer pack specialized function templates for supported Kokkos array
   /// ranks. Select out the proper elements from the input Array to send to a
   /// neighboring task and pack them into the proper send buffer for
   /// that Neighbor.
   template <typename T>
   std::enable_if_t<ArrayRank<T>::Is1D>
   packBuffer(const T &Array,      // 1D Kokkos array of any type
              const I4 CurNeighbor // current neighbor
   ) {

      OMEGA_SCOPE(LocList, Neighbors[CurNeighbor].SendLists[CurElem]);
      OMEGA_SCOPE(LocNeighbor, Neighbors[CurNeighbor]);

      if (devBufferPUP(Array)) {
         OMEGA_SCOPE(LocIndex, LocList.Index);
         Kokkos::resize(LocNeighbor.SendBuffer, LocList.NTot);
         OMEGA_SCOPE(LocBuff, Neighbors[CurNeighbor].SendBuffer);
         parallelFor(
             {LocList.NTot}, KOKKOS_LAMBDA(int IExch) {
                auto Val       = Array(LocIndex(IExch));
                const R8 RVal  = reinterpret_cast<R8 &>(Val);
                LocBuff(IExch) = RVal;
             });
      } else {
         OMEGA_SCOPE(LocIndexH, LocList.IndexH);
         Kokkos::resize(LocNeighbor.SendBufferH, LocList.NTot);
         OMEGA_SCOPE(LocBuffH, Neighbors[CurNeighbor].SendBufferH);
         for (int IExch = 0; IExch < LocList.NTot; ++IExch) {
            const R8 RVal   = reinterpret_cast<R8 &>(Array(LocIndexH(IExch)));
            LocBuffH(IExch) = RVal;
         }
      }
   }

   template <typename T>
   std::enable_if_t<ArrayRank<T>::Is2D>
   packBuffer(const T &Array,      // 2D Kokkos array of any type
              const I4 CurNeighbor // current neighbor
   ) {

      OMEGA_SCOPE(LocList, Neighbors[CurNeighbor].SendLists[CurElem]);
      OMEGA_SCOPE(LocNeighbor, Neighbors[CurNeighbor]);

      const I4 NJ = Array.extent(1);

      if (devBufferPUP(Array)) {
         OMEGA_SCOPE(LocIndex, LocList.Index);
         Kokkos::resize(LocNeighbor.SendBuffer, LocList.NTot * TotSize);
         OMEGA_SCOPE(LocBuff, Neighbors[CurNeighbor].SendBuffer);

         parallelFor(
             {LocList.NTot, NJ}, KOKKOS_LAMBDA(int IExch, int J) {
                auto Val       = Array(LocIndex(IExch), J);
                const R8 RVal  = reinterpret_cast<R8 &>(Val);
                const I4 IBuff = IExch * NJ + J;
                LocBuff(IBuff) = RVal;
             });
      } else {
         OMEGA_SCOPE(LocIndexH, LocList.IndexH);
         Kokkos::resize(LocNeighbor.SendBufferH, LocList.NTot * TotSize);
         OMEGA_SCOPE(LocBuffH, Neighbors[CurNeighbor].SendBufferH);
         for (int IExch = 0; IExch < LocList.NTot; ++IExch) {
            for (int J = 0; J < NJ; ++J) {
               const I4 IBuff = IExch * NJ + J;
               const R8 RVal =
                   reinterpret_cast<R8 &>(Array(LocIndexH(IExch), J));
               LocBuffH(IBuff) = RVal;
            }
         }
      }
   }

   template <typename T>
   std::enable_if_t<ArrayRank<T>::Is3D>
   packBuffer(const T &Array,      // 3D Kokkos array of any type
              const I4 CurNeighbor // current neighbor
   ) {

      OMEGA_SCOPE(LocList, Neighbors[CurNeighbor].SendLists[CurElem]);
      OMEGA_SCOPE(LocNeighbor, Neighbors[CurNeighbor]);

      const I4 NJ = Array.extent(2);
      const I4 NK = Array.extent(0);

      const I4 NTotList = LocList.NTot;

      if (devBufferPUP(Array)) {
         OMEGA_SCOPE(LocIndex, LocList.Index);
         Kokkos::resize(LocNeighbor.SendBuffer, LocList.NTot * TotSize);
         OMEGA_SCOPE(LocBuff, Neighbors[CurNeighbor].SendBuffer);

         parallelFor(
             {NK, NTotList, NJ}, KOKKOS_LAMBDA(int K, int IExch, int J) {
                auto Val       = Array(K, LocIndex(IExch), J);
                const R8 RVal  = reinterpret_cast<R8 &>(Val);
                const I4 IBuff = (K * NTotList + IExch) * NJ + J;
                LocBuff(IBuff) = RVal;
             });
      } else {
         OMEGA_SCOPE(LocIndexH, LocList.IndexH);
         Kokkos::resize(LocNeighbor.SendBufferH, LocList.NTot * TotSize);
         OMEGA_SCOPE(LocBuffH, Neighbors[CurNeighbor].SendBufferH);
         for (int K = 0; K < NK; ++K) {
            for (int IExch = 0; IExch < LocList.NTot; ++IExch) {
               for (int J = 0; J < NJ; ++J) {
                  const I4 IBuff = (K * LocList.NTot + IExch) * NJ + J;
                  const R8 RVal =
                      reinterpret_cast<R8 &>(Array(K, LocIndexH(IExch), J));
                  LocBuffH(IBuff) = RVal;
               }
            }
         }
      }
   }

   template <typename T>
   std::enable_if_t<ArrayRank<T>::Is4D>
   packBuffer(const T &Array,      // 4D Kokkos array of any type
              const I4 CurNeighbor // current neighbor
   ) {

      OMEGA_SCOPE(LocList, Neighbors[CurNeighbor].SendLists[CurElem]);
      OMEGA_SCOPE(LocNeighbor, Neighbors[CurNeighbor]);

      const I4 NJ = Array.extent(3);
      const I4 NK = Array.extent(1);
      const I4 NL = Array.extent(0);

      const I4 NTotList = LocList.NTot;

      if (devBufferPUP(Array)) {
         OMEGA_SCOPE(LocIndex, LocList.Index);
         Kokkos::resize(LocNeighbor.SendBuffer, LocList.NTot * TotSize);
         OMEGA_SCOPE(LocBuff, Neighbors[CurNeighbor].SendBuffer);

         parallelFor(
             {NL, NK, NTotList, NJ},
             KOKKOS_LAMBDA(int L, int K, int IExch, int J) {
                auto Val       = Array(L, K, LocIndex(IExch), J);
                const R8 RVal  = reinterpret_cast<R8 &>(Val);
                const I4 IBuff = ((L * NK + K) * NTotList + IExch) * NJ + J;
                LocBuff(IBuff) = RVal;
             });
      } else {
         OMEGA_SCOPE(LocIndexH, LocList.IndexH);
         Kokkos::resize(LocNeighbor.SendBufferH, LocList.NTot * TotSize);
         OMEGA_SCOPE(LocBuffH, Neighbors[CurNeighbor].SendBufferH);
         for (int L = 0; L < NL; ++L) {
            for (int K = 0; K < NK; ++K) {
               for (int IExch = 0; IExch < LocList.NTot; ++IExch) {
                  for (int J = 0; J < NJ; ++J) {
                     const I4 IBuff =
                         ((L * NK + K) * NTotList + IExch) * NJ + J;
                     const R8 RVal = reinterpret_cast<R8 &>(
                         Array(L, K, LocIndexH(IExch), J));
                     LocBuffH(IBuff) = RVal;
                  }
               }
            }
         }
      }
   }

   template <typename T>
   std::enable_if_t<ArrayRank<T>::Is5D>
   packBuffer(const T &Array,      // 5D Kokkos array of any type
              const I4 CurNeighbor // current neighbor
   ) {

      OMEGA_SCOPE(LocList, Neighbors[CurNeighbor].SendLists[CurElem]);
      OMEGA_SCOPE(LocNeighbor, Neighbors[CurNeighbor]);

      const I4 NJ = Array.extent(4);
      const I4 NK = Array.extent(2);
      const I4 NL = Array.extent(1);
      const I4 NM = Array.extent(0);

      const I4 NTotList = LocList.NTot;

      if (devBufferPUP(Array)) {
         OMEGA_SCOPE(LocIndex, LocList.Index);
         Kokkos::resize(LocNeighbor.SendBuffer, LocList.NTot * TotSize);
         OMEGA_SCOPE(LocBuff, Neighbors[CurNeighbor].SendBuffer);

         parallelFor(
             {NM, NL, NK, NTotList, NJ},
             KOKKOS_LAMBDA(int M, int L, int K, int IExch, int J) {
                auto Val      = Array(M, L, K, LocIndex(IExch), J);
                const R8 RVal = reinterpret_cast<R8 &>(Val);
                const I4 IBuff =
                    (((M * NL + L) * NK + K) * NTotList + IExch) * NJ + J;
                LocBuff(IBuff) = RVal;
             });
      } else {
         OMEGA_SCOPE(LocIndexH, LocList.IndexH);
         Kokkos::resize(LocNeighbor.SendBufferH, LocList.NTot * TotSize);
         OMEGA_SCOPE(LocBuffH, Neighbors[CurNeighbor].SendBufferH);
         for (int M = 0; M < NM; ++M) {
            for (int L = 0; L < NL; ++L) {
               for (int K = 0; K < NK; ++K) {
                  for (int IExch = 0; IExch < LocList.NTot; ++IExch) {
                     for (int J = 0; J < NJ; ++J) {
                        const I4 IBuff =
                            (((M * NL + L) * NK + K) * NTotList + IExch) * NJ +
                            J;
                        const R8 RVal = reinterpret_cast<R8 &>(
                            Array(M, L, K, LocIndexH(IExch), J));
                        LocBuffH(IBuff) = RVal;
                     }
                  }
               }
            }
         }
      }
   }

   /// Buffer unpack specialized function templates for supported Kokkos array
   /// ranks. After receiving a message from a neighboring task, save the
   /// elements of the proper receive buffer for that Neighbor into the
   /// corresponding halo elements of the input Array
   template <typename T>
   std::enable_if_t<ArrayRank<T>::Is1D>
   unpackBuffer(const T &Array,      // 1D Kokkos array of any type
                const I4 CurNeighbor // current neighbor
   ) {

      using ValType = typename T::non_const_value_type;

      OMEGA_SCOPE(LocList, Neighbors[CurNeighbor].RecvLists[CurElem]);
      OMEGA_SCOPE(LocNeighbor, Neighbors[CurNeighbor]);

      if (devBufferPUP(Array)) {
         OMEGA_SCOPE(LocIndex, LocList.Index);
         OMEGA_SCOPE(LocBuff, Neighbors[CurNeighbor].RecvBuffer);
         parallelFor(
             {LocList.NTot}, KOKKOS_LAMBDA(int IExch) {
                const I4 IArr = LocIndex(IExch);
                Array(IArr)   = reinterpret_cast<ValType &>(LocBuff(IExch));
             });
      } else {
         OMEGA_SCOPE(LocIndexH, LocList.IndexH);
         OMEGA_SCOPE(LocBuffH, Neighbors[CurNeighbor].RecvBufferH);
         for (int IExch = 0; IExch < LocList.NTot; ++IExch) {
            const I4 IArr = LocIndexH(IExch);
            Array(IArr)   = reinterpret_cast<ValType &>(LocBuffH(IExch));
         }
      }
   }

   template <typename T>
   std::enable_if_t<ArrayRank<T>::Is2D>
   unpackBuffer(const T &Array,      // 2D Kokkos array of any type
                const I4 CurNeighbor // current neighbor
   ) {

      using ValType = typename T::non_const_value_type;

      OMEGA_SCOPE(LocList, Neighbors[CurNeighbor].RecvLists[CurElem]);
      OMEGA_SCOPE(LocNeighbor, Neighbors[CurNeighbor]);

      const I4 NJ = Array.extent(1);

      if (devBufferPUP(Array)) {
         OMEGA_SCOPE(LocIndex, LocList.Index);
         OMEGA_SCOPE(LocBuff, Neighbors[CurNeighbor].RecvBuffer);

         parallelFor(
             {LocList.NTot, NJ}, KOKKOS_LAMBDA(int IExch, int J) {
                const I4 IBuff = IExch * NJ + J;
                const I4 IArr  = LocIndex(IExch);
                Array(IArr, J) = reinterpret_cast<ValType &>(LocBuff(IBuff));
             });
      } else {
         OMEGA_SCOPE(LocIndexH, LocList.IndexH);
         OMEGA_SCOPE(LocBuffH, Neighbors[CurNeighbor].RecvBufferH);
         for (int IExch = 0; IExch < LocList.NTot; ++IExch) {
            for (int J = 0; J < NJ; ++J) {
               const I4 IBuff = IExch * NJ + J;
               const I4 IArr  = LocIndexH(IExch);
               Array(IArr, J) = reinterpret_cast<ValType &>(LocBuffH(IBuff));
            }
         }
      }
   }

   template <typename T>
   std::enable_if_t<ArrayRank<T>::Is3D>
   unpackBuffer(const T &Array,      // 3D Kokkos array of any type
                const I4 CurNeighbor // current neighbor
   ) {

      using ValType = typename T::non_const_value_type;

      OMEGA_SCOPE(LocList, Neighbors[CurNeighbor].RecvLists[CurElem]);
      OMEGA_SCOPE(LocNeighbor, Neighbors[CurNeighbor]);

      const I4 NJ = Array.extent(2);
      const I4 NK = Array.extent(0);

      const I4 NTotList = LocList.NTot;

      if (devBufferPUP(Array)) {
         OMEGA_SCOPE(LocIndex, LocList.Index);
         OMEGA_SCOPE(LocBuff, Neighbors[CurNeighbor].RecvBuffer);

         parallelFor(
             {NK, NTotList, NJ}, KOKKOS_LAMBDA(int K, int IExch, int J) {
                const I4 IBuff    = (K * NTotList + IExch) * NJ + J;
                const I4 IArr     = LocIndex(IExch);
                Array(K, IArr, J) = reinterpret_cast<ValType &>(LocBuff(IBuff));
             });
      } else {
         OMEGA_SCOPE(LocIndexH, LocList.IndexH);
         OMEGA_SCOPE(LocBuffH, Neighbors[CurNeighbor].RecvBufferH);
         for (int K = 0; K < NK; ++K) {
            for (int IExch = 0; IExch < LocList.NTot; ++IExch) {
               for (int J = 0; J < NJ; ++J) {
                  const I4 IBuff = (K * LocList.NTot + IExch) * NJ + J;
                  const I4 IArr  = LocIndexH(IExch);
                  Array(K, IArr, J) =
                      reinterpret_cast<ValType &>(LocBuffH(IBuff));
               }
            }
         }
      }
   }

   template <typename T>
   std::enable_if_t<ArrayRank<T>::Is4D>
   unpackBuffer(const T &Array,      // 4D Kokkos array of any type
                const I4 CurNeighbor // current neighbor
   ) {

      using ValType = typename T::non_const_value_type;

      OMEGA_SCOPE(LocList, Neighbors[CurNeighbor].RecvLists[CurElem]);
      OMEGA_SCOPE(LocNeighbor, Neighbors[CurNeighbor]);

      const I4 NJ = Array.extent(3);
      const I4 NK = Array.extent(1);
      const I4 NL = Array.extent(0);

      const I4 NTotList = LocList.NTot;

      if (devBufferPUP(Array)) {
         OMEGA_SCOPE(LocIndex, LocList.Index);
         OMEGA_SCOPE(LocBuff, Neighbors[CurNeighbor].RecvBuffer);

         parallelFor(
             {NL, NK, NTotList, NJ},
             KOKKOS_LAMBDA(int L, int K, int IExch, int J) {
                const I4 IBuff = ((L * NK + K) * NTotList + IExch) * NJ + J;
                const I4 IArr  = LocIndex(IExch);
                Array(L, K, IArr, J) =
                    reinterpret_cast<ValType &>(LocBuff(IBuff));
             });
      } else {
         OMEGA_SCOPE(LocIndexH, LocList.IndexH);
         OMEGA_SCOPE(LocBuffH, Neighbors[CurNeighbor].RecvBufferH);
         for (int L = 0; L < NL; ++L) {
            for (int K = 0; K < NK; ++K) {
               for (int IExch = 0; IExch < LocList.NTot; ++IExch) {
                  for (int J = 0; J < NJ; ++J) {
                     const I4 IBuff =
                         ((L * NK + K) * NTotList + IExch) * NJ + J;
                     const I4 IArr = LocIndexH(IExch);
                     Array(L, K, IArr, J) =
                         reinterpret_cast<ValType &>(LocBuffH(IBuff));
                  }
               }
            }
         }
      }
   }

   template <typename T>
   std::enable_if_t<ArrayRank<T>::Is5D>
   unpackBuffer(const T &Array,      // 5D Kokkos array of any type
                const I4 CurNeighbor // current neighbor
   ) {

      using ValType = typename T::non_const_value_type;

      OMEGA_SCOPE(LocList, Neighbors[CurNeighbor].RecvLists[CurElem]);
      OMEGA_SCOPE(LocNeighbor, Neighbors[CurNeighbor]);

      const I4 NJ = Array.extent(4);
      const I4 NK = Array.extent(2);
      const I4 NL = Array.extent(1);
      const I4 NM = Array.extent(0);

      const I4 NTotList = LocList.NTot;

      if (devBufferPUP(Array)) {
         OMEGA_SCOPE(LocIndex, LocList.Index);
         OMEGA_SCOPE(LocBuff, Neighbors[CurNeighbor].RecvBuffer);

         parallelFor(
             {NM, NL, NK, NTotList, NJ},
             KOKKOS_LAMBDA(int M, int L, int K, int IExch, int J) {
                const I4 IBuff =
                    (((M * NL + L) * NK + K) * NTotList + IExch) * NJ + J;
                const I4 IArr = LocIndex(IExch);
                Array(M, L, K, IArr, J) =
                    reinterpret_cast<ValType &>(LocBuff(IBuff));
             });
      } else {
         OMEGA_SCOPE(LocIndexH, LocList.IndexH);
         OMEGA_SCOPE(LocBuffH, Neighbors[CurNeighbor].RecvBufferH);
         for (int M = 0; M < NM; ++M) {
            for (int L = 0; L < NL; ++L) {
               for (int K = 0; K < NK; ++K) {
                  for (int IExch = 0; IExch < LocList.NTot; ++IExch) {
                     for (int J = 0; J < NJ; ++J) {
                        const I4 IBuff =
                            (((M * NL + L) * NK + K) * NTotList + IExch) * NJ +
                            J;
                        const I4 IArr = LocIndexH(IExch);
                        Array(M, L, K, IArr, J) =
                            reinterpret_cast<ValType &>(LocBuffH(IBuff));
                     }
                  }
               }
            }
         }
      }
   }

   //---------------------------------------------------------------------------
   // Function template to perform a full halo exchange on the input Kokkos
   // array of any supported type defined on the input index space ThisElem
   template <typename T>
   int
   exchangeFullArrayHalo(T &Array,            // Kokkos array of any type
                         MeshElement ThisElem // index space Array is defined on
   ) {

      I4 IErr{0}; // error code

      // Logical flag to track if all messages have been received
      bool AllReceived{false};

      // Save the index space the input array is defined on
      CurElem = ThisElem;

      // For cell-based quantities, the number of halo layers equals HaloWidth,
      // edge- and vertex-based quantities have an extra layer.
      if (CurElem == OnCell) {
         NumLayers = HaloWidth;
      } else {
         NumLayers = HaloWidth + 1;
      }

      // Determine the number of array elements per cell, edge, or vertex
      // in the input array
      I4 NDims = Array.Rank;
      if (NDims == 1) {
         TotSize = 1;
      } else if (NDims == 2) {
         TotSize = Array.extent(1);
      } else {
         TotSize = 1;
         for (int I = 0; I < NDims - 2; ++I) {
            TotSize *= Array.extent(I);
         }
         TotSize *= Array.extent(NDims - 1);
      }

      // If the Array is in device memory space, the buffer pack and unpack
      // functions will use device buffers
      bool UseDevBuffer = devBufferPUP(Array);

      // Allocate the receive buffers and Call MPI_Irecv for each Neighbor
      // so the local task is ready to accept messages from each
      // neighboring task
      startReceives(UseDevBuffer);

      // Loop through each Neighbor, resetting communication flags and packing
      // buffers if there are elements to be sent to the neighboring task
      for (int INghbr = 0; INghbr < NNghbr; ++INghbr) {
         Neighbors[INghbr].Received = false;
         Neighbors[INghbr].Unpacked = false;
         if (SendFlags[CurElem][INghbr]) {
            packBuffer(Array, INghbr);
         }
      }

      // Call MPI_Isend for each Neighbor to send the packed buffers
      startSends(UseDevBuffer);

      // Wait for all sends to complete before proceeding
      for (int INghbr = 0; INghbr < NNghbr; ++INghbr) {
         if (SendFlags[CurElem][INghbr]) {
            MPI_Wait(&Neighbors[INghbr].SReq, MPI_STATUS_IGNORE);
         }
      }

      I4 MaxIter = 1000000000; // Large integer to prevent infinite loop
      I4 IPass   = 0;          // Number of passes through while loop
      I4 NRcvd   = 0;          // Integer to track number of messages received

      // Total number of messages the local task will receive
      I4 NMessages = std::accumulate(RecvFlags[CurElem].begin(),
                                     RecvFlags[CurElem].end(), 0);

      // Until all messages from neighboring tasks are received, loop
      // through Neighbor objects and use MPI_Test to check if the message
      // has been received. Unpack buffers upon receipt of each message
      while (!AllReceived) {
         for (int INghbr = 0; INghbr < NNghbr; ++INghbr) {
            if (RecvFlags[CurElem][INghbr]) {
               if (!Neighbors[INghbr].Received) {
                  MPI_Test(&Neighbors[INghbr].RReq, &Neighbors[INghbr].Received,
                           MPI_STATUS_IGNORE);
                  if (Neighbors[INghbr].Received) {
                     ++NRcvd;
                  }
               }
               if (Neighbors[INghbr].Received && !Neighbors[INghbr].Unpacked) {
                  // If the device buffer will be used in unpackBuffer, but the
                  // exchange was done with the host buffer, a deep copy from
                  // host to device is needed.
                  if (UseDevBuffer && !ExchOnDev) {
                     Kokkos::resize(Neighbors[INghbr].RecvBuffer,
                                    Neighbors[INghbr].RecvBufferH.size());
                     deepCopy(Neighbors[INghbr].RecvBuffer,
                              Neighbors[INghbr].RecvBufferH);
                  }
                  unpackBuffer(Array, INghbr);
                  Neighbors[INghbr].Unpacked = true;
               }
            }
         }

         if (NRcvd == NMessages) {
            AllReceived = true;
         }
         ++IPass;
         if (IPass == MaxIter) {
            LOG_ERROR("Halo: Maximum iterations reached during halo exchange");
            IErr = -1;
            break;
         }
      }

      return IErr;
   } // end exchangeFullArrayHalo

}; // end class Halo

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_HALO_H
