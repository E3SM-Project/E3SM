#ifndef OMEGA_ANALYSIS_H
#define OMEGA_ANALYSIS_H

//===-- analysis/Analysis.h - OMEGA Analysis ----------*- C++ -*-===//
//
/// \file
/// \brief Defines core Analysis framework for in-situ computation
///
/// The Analysis module provides in-situ computation of analysis fields from
/// the ocean model state during simulation runtime. Analysis fields are
/// computed on-the-fly and written to output streams at user-specified
/// intervals, providing an alternative to extensive offline
/// post-processing.
///
/// The framework is built on a composable operator architecture where
/// operators can be chained together to produce analysis outputs. Each
/// operator performs a single, well-defined transformation (e.g., spatial
/// reduction, temporal averaging, binary operations). Operators declare
/// their input field dependencies at construction, and the Analysis
/// orchestrator resolves dependencies to form a directed acyclic graph
/// (DAG).
///
/// The orchestrator uses an alarm-based scheduling model to trigger
/// operator computation. Each operator node contains pointers to one or
/// more alarms; when any alarm rings, the operator computes its output.
/// Upstream dependencies are computed recursively on-demand, with
/// timestamp-based caching to prevent redundant computation when multiple
/// downstream operators share an intermediate result.
///
/// For detailed design documentation including algorithmic formulation,
/// dependency resolution, and operator registration, see
/// components/omega/doc/design/Analysis.md
///
//===----------------------------------------------------------------------===//

#include "AnalysisOpFactory.h"
#include "AnalysisOperator.h"
#include "Config.h"
#include "DataTypes.h"
#include "Dimension.h"
#include "Error.h"
#include "Field.h"
#include "HorzMesh.h"
#include "Logging.h"
#include "MachEnv.h"
#include "Pacer.h"
#include "TimeMgr.h"
#include "VertCoord.h"
#include "operators/Ops.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace OMEGA {

/// Parses a frequency string into numeric and unit components
/// Frequency strings are of the form "1day", "6hour", "1month", etc.
/// Returns a vector with two elements: [0] = numeric part, [1] = units
/// If units do not end with 's', appends 's'.
std::vector<std::string>
parseFreqStr(const std::string &FreqStr ///< [in] frequency string to parse
);

/// Internal representation of a node in the Analysis operator dependency
/// graph. Each node contains an operator instance, pointers to its upstream
/// dependencies, the names of output streams that consume its output, and
/// non-owning pointers to alarms that trigger its computation.
struct OperatorNode {
   std::unique_ptr<AnalysisOperator> Op;  ///< Operator instance (owned)
   std::vector<OperatorNode *> Upstreams; ///< Upstream deps (non-owning)
   std::vector<std::string> StreamNames;  ///< Output streams; empty if
                                          ///< intermediate operator
   std::vector<Alarm *> ComputeAlarms;    ///< Compute alarms (non-owning)
};

/// The Analysis class is the top-level orchestrator for the in-situ
/// analysis framework. It is responsible for:
/// - Reading analysis configuration and constructing AnalysisGroup instances
/// - Parsing operator chain strings and instantiating operators via factory
/// - Resolving operator dependencies to form a directed acyclic graph (DAG)
/// - Managing the alarm-based scheduling system for operator computation
/// - Owning accumulation alarms for temporal reduction operators
/// - Triggering recursive operator computation on each timestep
///
/// The class maintains a vector of OperatorNode instances representing all
/// registered operators. Each node stores pointers to its upstream
/// dependencies and the alarms that trigger its computation. The Analysis
/// object owns accumulation alarms for temporal reduction operators, while
/// output alarms are borrowed from IOStream instances.
class Analysis {
 public:
   /// Initializes the Analysis module by registering all base operators
   /// and creating the default Analysis instance. This function must be
   /// called after HorzMesh, VertCoord, and TimeStepper are initialized,
   /// as it retrieves pointers to the default instances of these objects.
   static void init();

   /// Creates a new named Analysis instance and registers it in the
   /// AllAnalysisObjects map. Each Analysis instance maintains its own
   /// set of operators and dependency graph.
   static Analysis *
   create(const std::string &Name, ///< [in] name for new Analysis instance
          const MachEnv *Env,      ///< [in] machine environment
          const HorzMesh *Mesh,    ///< [in] horizontal mesh
          const VertCoord *VCoord, ///< [in] vertical coordinate
          Clock *ModelClock,       ///< [in] pointer to model clock
          Config *Options          ///< [in] configuration options
   );

   /// Computes all analysis fields whose alarms are ringing at the current
   /// timestep. This function is called once per timestep from the main
   /// driver. It iterates over all operator nodes and triggers recursive
   /// computation for any node with a ringing alarm.
   void computeAll();

   /// Parses an underscore-delimited operator chain string and instantiates
   /// all operators in the chain that do not yet exist as Fields. For
   /// example, "Temperature_SpatialMean_TimeMean1day" parses into three
   /// operators. If intermediate operators already exist (shared by other
   /// chains), they are reused rather than duplicated.
   void parseChainAndBuildOps(
       const std::string &OpChainStr ///< [in] underscore-delimited chain string
   );

   /// Instantiates a single operator via the factory and appends it as
   /// an OperatorNode. The factory selects the correct templated
   /// specialization based on the primary upstream Field's metadata
   /// (scalar type, rank, memory location).
   void
   registerAnalysisOp(const std::string &OpName, ///< [in] operator type name
                      const std::vector<std::string>
                          &UpstreamNames, ///< [in] upstream field names
                      Config Options      ///< [in] operator configuration
   );

   /// Returns a reference to the model clock pointer. This is used by
   /// AnalysisGroup instances during stream creation to access the clock
   /// for alarm initialization.
   Clock *&getModelClock();

   /// Returns a vector of non-owning pointers to all registered operator
   /// nodes. Used primarily for testing and validation of the dependency graph.
   const std::vector<OperatorNode *> getOpNodes();

   /// Checks whether an operator node with the given full instance name
   /// already exists in the OpNodes vector. Returns true if found, false
   /// otherwise. Used during chain parsing to avoid creating duplicate
   /// operators.
   bool OpNodeExists(const std::string &FullOpName ///< [in] full op name
   );

   /// Retrieves the default Analysis instance. The preference is to pass
   /// the Analysis pointer as an argument, but retrieval is necessary for
   /// sharing info between initialization and run phases.
   static Analysis *getDefault();

   /// Destructor - deallocates all memory and deletes the Analysis instance
   ~Analysis();

   /// Removes all defined Analysis instances and cleans up static resources
   static void finalize();

 private:
   /// Accumulation alarms owned by Analysis for temporal reduction
   /// operators. These alarms control how frequently samples are added to
   /// running sums (e.g., every timestep or at a coarser interval). Output
   /// alarms are borrowed from IOStream instances and stored in
   /// OperatorNode::ComputeAlarms.
   std::vector<std::unique_ptr<Alarm>> AccumulationAlarms;

   /// Pointer to the default Analysis instance for easy retrieval
   static Analysis *DefAnalysis;

   /// Map of all Analysis instances, paired with names for retrieval
   static std::map<std::string, std::unique_ptr<Analysis>> AllAnalysisObjects;

   /// Private constructor for creating a new Analysis instance.
   /// Called by the create() factory method. Reads configuration,
   /// constructs AnalysisGroup instances, and builds the operator
   /// dependency graph.
   Analysis(const std::string &Name, ///< [in] name for new instance
            const MachEnv *Env,      ///< [in] machine environment
            const HorzMesh *Mesh,    ///< [in] horizontal mesh
            const VertCoord *VCoord, ///< [in] vertical coordinate
            Clock *ModelClock,       ///< [in] pointer to model clock
            Config *Options          ///< [in] configuration options
   );

   std::string Name;        ///< Name of this Analysis instance
   const MachEnv *Env;      ///< Machine environment for MPI operations
   const HorzMesh *Mesh;    ///< Horizontal mesh for spatial operations
   const VertCoord *VCoord; ///< Vertical coordinate for vertical ops
   Clock *ModelClock;       ///< Pointer to model clock for time mgmt

   /// All registered operator nodes forming the dependency graph
   std::vector<std::unique_ptr<OperatorNode>> OpNodes;

   /// Registers all built-in operator types with the AnalysisOpFactory.
   /// Called once during init() before any operators are instantiated.
   /// Each operator template is registered with all supported array type
   /// variants (scalar types, ranks, memory locations). Defined in Ops.cpp.
   static void registerAllBaseAnalysisOperators();

   /// Post-hoc dependency resolution: iterates over all operator nodes and
   /// matches input field names against other nodes' output field names to
   /// populate the Upstreams vectors. This forms the edges of the
   /// dependency graph. In future versions, this will be replaced by
   /// signature-based deduplication during graph construction.
   void buildOperatorDependencies();

   /// Sets ComputeAlarms on terminal nodes by borrowing alarm pointers from
   /// associated IOStream instances. For temporal reduction operators, also
   /// creates accumulation alarms and adds them to ComputeAlarms. Then
   /// calls propagateAlarmsUpstream() to propagate alarms to upstream
   /// dependencies.
   void setComputeAlarms();

   /// Calls initialize() on all operators after the dependency graph is
   /// complete and all Fields exist. This allows operators to store
   /// pointers to mesh, environment, and other resources needed during
   /// compute().
   void initializeAllOps();

   /// Iteratively propagates alarm pointers from downstream operators to
   /// upstream operators. An upstream operator must be computed whenever
   /// any of its downstream consumers needs data, so each downstream alarm
   /// is added to the upstream's ComputeAlarms vector. Iteration continues
   /// until no further changes occur (fixed point).
   void propagateAlarmsUpstream();

   /// Recursively computes an operator node and all its upstream
   /// dependencies. Uses timestamp-based cache validation to prevent
   /// redundant computation when multiple downstream operators share an
   /// intermediate result. If the node's LastComputed timestamp matches the
   /// current TimeStamp, returns immediately (cache hit). Otherwise,
   /// recursively computes all upstream dependencies first, then calls the
   /// node's compute() method.
   void
   computeRecursive(OperatorNode *Node,          ///< [in] node to compute
                    const TimeInstant &TimeStamp ///< [in] current timestamp
   );

   // Forbid copy and move construction
   Analysis(const Analysis &) = delete;
   Analysis(Analysis &&)      = delete;

}; // end class Analysis

} // namespace OMEGA
//===----------------------------------------------------------------------===//
#endif // OMEGA_ANALYSIS_H
