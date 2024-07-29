#ifndef OMEGA_MACH_ENV_H
#define OMEGA_MACH_ENV_H
//===-- base/MachEnv.h - machine environment definitions --------*- C++ -*-===//
//
/// \file
/// \brief Defines aspects of the parallel and node machine environment
///
/// The machine environment defines a number of parameters associated with
/// the message-passing and threading environments. It also can describe
/// the node-level hardware environment, including any useful accelerator
/// or processor-level constants. Multiple machine environments can be defined,
/// mostly associated with subsets of the processor decomposition, but
/// a default machine env must be created very early in model
/// initialization. All environments can be retrieved by name, though a
/// specific retrieval for the most common default environment is provided
/// for efficiency.
//
//===----------------------------------------------------------------------===//

#include "mpi.h"

#include "Logging.h"
#include <map>
#include <memory>

namespace OMEGA {

// For CPUs, a compile-time vector length is useful for
// blocking the inner loops. This should be set to an appropriate
// length (typically 32, 64, 128) for CPU-only builds, but set to one for
// GPU builds to maximize parallelism instead.
#ifdef OMEGA_VECTOR_LENGTH
constexpr int VecLength = OMEGA_VECTOR_LENGTH;
#else
constexpr int VecLength = 1;
#endif

/// The MachEnv class is a container that holds information on
/// the message passing, threading and node environment.
class MachEnv {

 private:
   MPI_Comm Comm;       ///< MPI communicator for this environment
   int MyTask;          ///< task ID for local MPI task (rank)
   int NumTasks;        ///< total number of MPI tasks (ranks)
   int MasterTask;      ///< task ID for master task
   bool MasterTaskFlag; ///< true if this task is the master task
   bool MemberFlag;     ///< true if task is in communicator group

   // Add threading variables here
   int NumThreads; ///< number of OpenMP threads per task

   // Add any other useful machine parameters here
   // It may be useful at some point to track the number
   // of various devices per node (CPUs, GPUs), the number
   // of tasks allocated per node, etc.

   /// The default environment describes the environment for OMEGA
   /// defined for most of the model. Because it is used most often,
   /// we store the extra pointer here for easier retrieval.
   static MachEnv *DefaultEnv;

   /// All environments are tracked/stored within the class as a
   /// map paired with a name for later retrieval.
   static std::map<std::string, std::unique_ptr<MachEnv>> AllEnvs;

   /// CONSTRUCTORS
   /// All constructors are declared private to prevent accidental creation
   /// of new environments, MachEnv::create is the only way to create a new
   /// environment

   /// Constructor for environment based on input MPI communicator
   MachEnv(const std::string Name, ///< [in] name of the environment
           const MPI_Comm inComm   ///< [in] MPI communicator to use
   );

   /// Constructs a new environment with a given name from a contiguous
   /// subset of tasks in an existing environment starting at task 0
   /// of the parent environment. The master task can optionally be
   /// set but will default to zero if not provided.
   MachEnv(const std::string Name,    ///< [in] name of environment
           const MachEnv *InEnv,      ///< [in] existing parent MachEnv
           const int NewSize,         ///< [in] use first newSize tasks
           const int InMasterTask = 0 ///< [in] optional task to use for master
   );

   /// Constructs a new environment with a given name from a strided
   /// subset of of tasks in an existing environment. The master task
   /// can optionally be set but will default to zero if not provided.
   MachEnv(const std::string Name,    ///< [in] name of env
           const MachEnv *InEnv,      ///< [in] existing parent MachEnv
           const int NewSize,         ///< [in] num tasks in new env
           const int Begin,           ///< [in] starting parent task
           const int Stride,          ///< [in] stride for tasks to incl
           const int InMasterTask = 0 ///< [in] optional task to use for master
   );

   /// Constructs a new environment with a given name from a custom subset
   /// of tasks in an existing environment. The tasks are defined by a
   /// list of parent tasks to include. The master task can optionally be
   /// set but will default to zero if not provided.
   MachEnv(const std::string Name,    ///< [in] name of environment
           const MachEnv *InEnv,      ///< [in] existing parent MachEnv
           const int NewSize,         ///< [in] num tasks in new env
           const int Tasks[],         ///< [in] vector of parent tasks to incl
           const int InMasterTask = 0 ///< [in] optional task to use for master
   );

   // forbid copy and move construction
   MachEnv(const MachEnv &) = delete;
   MachEnv(MachEnv &&)      = delete;

 public:
   // Methods

   // Creates a new environment by calling the constructor with the
   // supplied arguments and stores it in the map of all environments.
   // This a variadic template because MachEnv constructors have different
   // numbers of arguments and we want to simply forward the supplied
   // arguments to the constructor
   template <class... ArgTypes>
   static MachEnv *create(const std::string &Name, ArgTypes &&...Args) {
      // Check to see if an environment of the same name already exists and
      // if so, exit with an error
      if (AllEnvs.find(Name) != AllEnvs.end()) {
         LOG_ERROR("Attempted to create a MachEnv with name {} but an Env of "
                   "that name already exists",
                   Name);
         return nullptr;
      }

      // create a new environment on the heap and put it in a map of
      // unique_ptrs, which will manage its lifetime
      auto *NewEnv = new MachEnv(Name, std::forward<ArgTypes>(Args)...);
      AllEnvs.emplace(Name, NewEnv);

      return get(Name);
   }

   /// Initializes the Machine Environment and creates the default
   /// machine environment based on an input MPI communicator. In
   /// standalone mode, this will typically be MPI_COMM_WORLD, but
   /// in coupled mode, this is the communicator assigned to the
   /// Omega component.
   static void init(const MPI_Comm InComm ///< [in] MPI communicator to use
   );

   /// Removes a MachEnv
   static void
   removeEnv(const std::string Name ///< [in] name of environment to remove
   );

   /// Removes all MachEnvs to clean up
   static void removeAll();

   // Retrieval functions

   /// Retrieve the default environment
   static MachEnv *getDefault();

   /// Retrieve any other environment by name
   static MachEnv *get(const std::string Name ///< [in] name of environment
   );

   /// Get communicator for an environment
   MPI_Comm getComm() const; ///< returns MPI communicator for this env

   /// Get local task/rank ID
   int getMyTask() const;

   /// Get total number of MPI tasks/ranks
   int getNumTasks() const;

   /// Get task ID for the master task (typically 0)
   int getMasterTask() const;

   /// Determine whether local task is the master
   bool isMasterTask() const;

   /// Determine whether local task is a member of this environment.
   /// This is primarily to prevent retrievals of non-existent
   /// values when a given environment uses only a subset of the
   /// tasks.
   bool isMember() const;

   // Only one variable can be set

   /// Set master task ID. By default, the master task is task 0 but
   /// can be set here to a different task if the master task has
   /// become over-burdened by work or memory (eg in coupled mode when
   /// all components are using the same master).
   int setMasterTask(const int TaskID ///< [in] new task to use as master
   );

   /// Prints all members of a MachEnv (typically for debugging)
   void print() const;

}; // end class MachEnv

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_MACH_ENV_H
