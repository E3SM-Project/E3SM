//===-- analysis/Analysis.cpp - analysis module -----------------*- C++ -*-===//
//
// The Analysis module provides in-situ computation of analysis fields from
// the ocean model state during simulation runtime. This implementation file
// contains the core orchestration logic including operator chain parsing,
// dependency graph construction, alarm-based scheduling, and recursive
// computation with caching.
//
//===----------------------------------------------------------------------===//

#include "Analysis.h"
#include "IOStream.h"
#include "TimeStepper.h"
#include "analysisGroups/Groups.h"

namespace OMEGA {

// Create the static class members
Analysis *Analysis::DefAnalysis = nullptr;
std::map<std::string, std::unique_ptr<Analysis>> Analysis::AllAnalysisObjects;

//------------------------------------------------------------------------------
// Parses a frequency string like "1day" or "6hour" into numeric and unit
// components. Returns a vector with [0] = numeric part, [1] = units part.
// Appends 's' to units if not already present for consistency with
// TimeInterval.
std::vector<std::string> parseFreqStr(const std::string &FreqStr) {

   std::string DigitStr;
   std::string UnitsStr;

   // Find where digits end and units begin
   size_t Pos = FreqStr.find_first_not_of("0123456789");
   if (Pos != std::string::npos) {
      DigitStr = FreqStr.substr(0, Pos);
      UnitsStr = FreqStr.substr(Pos);
   }

   // Validate that we found both components
   if (DigitStr.empty() || UnitsStr.empty()) {
      ABORT_ERROR("Analysis: Invalid frequency string found in Config: {}",
                  FreqStr);
   }

   // Ensure units end with 's' for consistency with TimeInterval
   if (UnitsStr.back() != 's') {
      UnitsStr += 's';
   }

   return {DigitStr, UnitsStr};

} // end parseFreqStr

//------------------------------------------------------------------------------
// Initializes the Analysis module by registering all base operator types
// with the factory, then creating the default Analysis instance. This function
// must be called after HorzMesh, VertCoord, and TimeStepper have been
// initialized, as it retrieves pointers to their default instances.
void Analysis::init() {

   Pacer::start("Analysis:init", 0);

   // Register all built-in operators with the factory
   registerAllBaseAnalysisOperators();

   // Retrieve default instances of required components
   auto DefEnv         = MachEnv::getDefault();
   auto Mesh           = HorzMesh::getDefault();
   auto VCoord         = VertCoord::getDefault();
   auto DefTimeStepper = TimeStepper::getDefault();
   Clock *OmegaClock   = DefTimeStepper->getClock();
   Config *OmegaConfig = Config::getOmegaConfig();

   // Create the default Analysis instance
   Analysis::DefAnalysis =
       create("Default", DefEnv, Mesh, VCoord, OmegaClock, OmegaConfig);

   Pacer::stop("Analysis:init", 0);

} // end init

//------------------------------------------------------------------------------
// Creates a new named Analysis instance using the constructor and registers
// it in the AllAnalysisObjects map. Returns a pointer to the new instance,
// or nullptr if an Analysis with this name already exists.
Analysis *Analysis::create(const std::string &Name, const MachEnv *Env,
                           const HorzMesh *Mesh, const VertCoord *VCoord,
                           Clock *ModelClock, Config *Options) {

   // Check for duplicate names
   if (AllAnalysisObjects.find(Name) != AllAnalysisObjects.end()) {
      LOG_ERROR("Attempted to create an Analysis with name {} "
                "but one with that name already exists",
                Name);
      return nullptr;
   }

   // Create new Analysis instance via constructor
   auto NewAnalysis =
       new Analysis(Name, Env, Mesh, VCoord, ModelClock, Options);

   // Store in map and return pointer
   AllAnalysisObjects.emplace(Name, NewAnalysis);

   return NewAnalysis;

} // end create

//------------------------------------------------------------------------------
// Constructs a new Analysis instance by reading the Analysis configuration
// node, instantiating enabled AnalysisGroup objects, building the operator
// dependency graph, setting up alarms, and initializing all operators.
// Pre-defined AnalysisGroup types (e.g., GlobalStats) are dispatched by name;
// user-defined custom groups will be supported in future versions.
Analysis::Analysis(const std::string &InName, const MachEnv *InEnv,
                   const HorzMesh *InMesh, const VertCoord *InVCoord,
                   Clock *InModelClock, Config *Options) {

   // Store member variables
   Name       = InName;
   Env        = InEnv;
   Mesh       = InMesh;
   VCoord     = InVCoord;
   ModelClock = InModelClock;

   // Set up name prefix for operator and stream names
   std::string NamePrefix = Name + "_";
   if (Name == "Default") {
      NamePrefix = "";
   }

   // Read the Analysis configuration node
   Error Err;
   Config AnalysisCfg("Analysis");
   Err += Options->get(AnalysisCfg);
   CHECK_ERROR_ABORT(Err, "Analysis: Analysis group not in Config");

   // Loop over all AnalysisGroup nodes in the config
   for (Config::Iter It = AnalysisCfg.begin(); It != AnalysisCfg.end(); ++It) {
      std::string GroupName;
      AnalysisCfg.getName(It, GroupName);
      Config GroupCfg(GroupName);
      Err += AnalysisCfg.get(GroupCfg);

      // Check if this group is enabled
      bool GroupEnabled = false;
      GroupCfg.get("Enable", GroupEnabled);

      if (GroupEnabled) {
         // Dispatch to pre-defined AnalysisGroup constructors by name
         if (GroupName == "GlobalStats") {
            GlobalStats GlobalStatsGroup(NamePrefix + GroupName, GroupCfg,
                                         this);
            continue;
         }

         // User-defined custom groups not yet supported
         ABORT_ERROR("Analysis: custom analysis group enabled in config, but "
                     "composable operators are not yet supported.");
      }
   }

   // Build the operator dependency graph
   buildOperatorDependencies();

   // Set up alarm-based scheduling
   setComputeAlarms();

   // Initialize all operators with mesh/env pointers
   initializeAllOps();

} // end constructor

//------------------------------------------------------------------------------
// Parses an underscore-delimited operator chain string and instantiates
// operators for each node in the chain. For example,
// "Temperature_SpatialMean_TimeMean1day" parses into three operators.
// If an intermediate operator already exists (shared by another chain),
// it is reused rather than duplicated. This natural sharing mechanism
// avoids redundant computation of common intermediate results.
void Analysis::parseChainAndBuildOps(const std::string &OpChainStr) {

   // Split the chain string on underscore delimiters
   std::vector<std::string> ChainVec;
   std::stringstream OpChainSS(OpChainStr);
   std::string Part;
   while (std::getline(OpChainSS, Part, '_')) {
      ChainVec.push_back(Part);
   }

   // Build the chain left-to-right, checking at each step if the
   // operator already exists (via Field registry). If not, create it.
   std::string CurChainStr;
   for (const auto &ChainNode : ChainVec) {
      std::string Upstream = CurChainStr;

      // Build up the current chain string progressively
      if (CurChainStr.empty()) {
         CurChainStr = ChainNode;
      } else {
         CurChainStr += ("_" + ChainNode);
      }

      // Check if this operator already exists by looking for its output Field
      if (!Field::exists(CurChainStr)) {

         // Spatial operators (SpatialMean, SpatialMax, etc.)
         if (ChainNode.find("Spatial") != std::string::npos) {
            registerAnalysisOp(ChainNode, {Upstream}, makeOpConfig());
            continue;
         }

         // Temporal operators (TimeMean, etc.) with period embedded in name
         if (ChainNode.find("Time") != std::string::npos) {
            // Extract operator type and period string (e.g., "TimeMean1day")
            std::size_t Pos = ChainNode.find_first_of("0123456789");
            if (Pos == std::string::npos) {
               ABORT_ERROR("Analysis: Improper temporal window string, {}",
                           ChainNode);
            }
            std::string TimeOp  = ChainNode.substr(0, Pos);
            std::string FreqStr = ChainNode.substr(Pos);
            registerAnalysisOp(TimeOp, {Upstream},
                               makeOpConfig(opParam("Period", FreqStr)));
            continue;
         }
         ABORT_ERROR("Analysis: Error trying to parse {}. No Field or "
                     "Operator named {}",
                     OpChainStr, ChainNode);
      }
   }

} // end parseChainAndBuildOps

//------------------------------------------------------------------------------
// Instantiates a single operator via the AnalysisOpFactory and wraps it in
// an OperatorNode. The factory inspects the primary upstream Field's metadata
// (scalar type, rank, memory location) to select the correct templated
// specialization of the operator. The new node is appended to the OpNodes
// vector.
void Analysis::registerAnalysisOp(const std::string &OpName,
                                  const std::vector<std::string> &UpstreamNames,
                                  Config Options) {

   // Create the operator via factory (type-dispatched based on Field metadata)
   auto NewOp = AnalysisOpFactory::createOp(OpName, UpstreamNames, Options);

   if (NewOp) {
      // Wrap in an OperatorNode and add to the graph
      auto Node = std::make_unique<OperatorNode>();
      Node->Op  = std::move(NewOp);
      OpNodes.push_back(std::move(Node));
   }

} // end registerAnalysisOp

//------------------------------------------------------------------------------
// Performs post-hoc dependency resolution by iterating over all operator nodes
// and matching each node's input field names against other nodes' output field
// names. When a match is found, adds a pointer to the Upstreams vector,
// forming the edges of the dependency graph. Avoids adding duplicate upstream
// pointers. In future versions, this will be replaced by signature-based
// deduplication during graph construction.
void Analysis::buildOperatorDependencies() {
   // For each operator node, find its upstream dependencies
   for (auto &Node : OpNodes) {
      auto InputNames = Node->Op->getInputFieldNames();

      // For each required input field
      for (const auto &InputName : InputNames) {
         // Search all other nodes for matching output
         for (auto &PotentialUpstream : OpNodes) {
            auto UpstreamOutputs = PotentialUpstream->Op->getOutputFieldNames();

            for (const auto &UpstreamOutput : UpstreamOutputs) {
               if (UpstreamOutput == InputName) {
                  // Found a match - check if already in Upstreams vector
                  if (std::find(Node->Upstreams.begin(), Node->Upstreams.end(),
                                PotentialUpstream.get()) ==
                      Node->Upstreams.end()) {
                     // Not found, so add the upstream dependency
                     Node->Upstreams.push_back(PotentialUpstream.get());
                  }
                  break;
               }
            }
         }
      }
   }

} // end buildOperatorDependencies

//------------------------------------------------------------------------------
// Sets ComputeAlarms on terminal operator nodes and propagates them upstream.
// Terminal nodes that write to output streams receive alarm pointers from
// those streams. Temporal reduction operators receive two alarms: an
// accumulation alarm (owned by Analysis) that triggers sample collection,
// and an output alarm (from the stream) that triggers finalization and write.
// Discrete-sampling operators receive only the stream's output alarm.
void Analysis::setComputeAlarms() {

   // Get model timestep and current time for alarm initialization
   TimeInterval Timestep   = ModelClock->getTimeStep();
   TimeInstant CurrentTime = ModelClock->getCurrentTime();

   // Set alarms on terminal operators (those that write to output streams)
   for (auto &Node : OpNodes) {
      std::string OpType   = Node->Op->getOperatorType();
      bool IsTimeReduction = (OpType.find("Time") != std::string::npos);

      if (!Node->StreamNames.empty()) {
         // This is a terminal operator that writes to one or more streams

         if (IsTimeReduction) {
            // Temporal reduction operators need two alarms:
            // 1. Accumulation alarm (owned by Analysis) - triggers sample
            // collection
            auto AccumulationAlarm = std::make_unique<Alarm>(
                "Compute_" + OpType, Timestep, CurrentTime);

            // Attach alarm to clock and store pointer
            ModelClock->attachAlarm(AccumulationAlarm.get());
            Node->ComputeAlarms.push_back(AccumulationAlarm.get());
            AccumulationAlarms.push_back(std::move(AccumulationAlarm));

            // 2. Output alarm (from stream) - triggers finalization and write
            // Time reduction operators have exactly one associated stream
            auto StreamAlarm = IOStream::getAlarm(Node->StreamNames[0]);
            Node->ComputeAlarms.push_back(StreamAlarm);

            // Tell the operator which alarm triggers period finalization
            Node->Op->setPeriodAlarm(StreamAlarm);

         } else {
            // Discrete sampling operators: borrow alarm pointers from streams
            for (const auto &StreamName : Node->StreamNames) {
               auto StreamAlarm = IOStream::getAlarm(StreamName);

               // Avoid adding duplicate alarm pointers
               if (std::find(Node->ComputeAlarms.begin(),
                             Node->ComputeAlarms.end(),
                             StreamAlarm) == Node->ComputeAlarms.end()) {
                  Node->ComputeAlarms.push_back(StreamAlarm);
               }
            }
         }
      }
   }

   // Propagate alarms from terminal nodes to all upstream dependencies
   propagateAlarmsUpstream();

} // end setComputeAlarms

//------------------------------------------------------------------------------
// Iteratively propagates alarm pointers from downstream operators to upstream
// operators. An upstream operator must be computed whenever any of its
// downstream consumers needs fresh data, so each downstream alarm is added
// to the upstream's ComputeAlarms vector. Iteration continues until no
// further changes occur (fixed point algorithm).
void Analysis::propagateAlarmsUpstream() {

   // Iterate until no changes occur (fixed point)
   bool Changed = true;
   while (Changed) {
      Changed = false;

      // For each node, find all downstream operators
      for (auto &Node : OpNodes) {
         for (const auto &OtherNode : OpNodes) {
            // Check if Node is in OtherNode's upstream list
            for (const auto *Upstream : OtherNode->Upstreams) {
               if (Upstream == Node.get()) {
                  // Node is upstream of OtherNode, so add OtherNode's alarms
                  for (auto *DownstreamAlarm : OtherNode->ComputeAlarms) {
                     // Only add if not already present
                     if (std::find(Node->ComputeAlarms.begin(),
                                   Node->ComputeAlarms.end(),
                                   DownstreamAlarm) ==
                         Node->ComputeAlarms.end()) {
                        Node->ComputeAlarms.push_back(DownstreamAlarm);
                        Changed = true;
                     }
                  }
               }
            }
         }
      }
   }

} // end propagateAlarmsUpstream

//------------------------------------------------------------------------------
// Main computational loop called once per timestep. Iterates over all
// operator nodes and checks if any of their alarms are ringing. If so,
// calls computeRecursive to compute the operator and all its upstream
// dependencies, using cache validation to avoid redundant work.
void Analysis::computeAll() {

   Pacer::start("Analysis:computeAll", 0);

   TimeInstant CurrentTime = ModelClock->getCurrentTime();

   // Check each operator node for ringing alarms
   for (auto &Node : OpNodes) {
      bool ShouldCompute = false;

      // Check if any of this node's alarms are ringing
      for (auto *Alarm : Node->ComputeAlarms) {
         if (Alarm->isRinging()) {
            ShouldCompute = true;
            break;
         }
      }

      // Trigger recursive computation if any alarm is ringing
      if (ShouldCompute) {
         computeRecursive(Node.get(), CurrentTime);
      }
   }

   Pacer::stop("Analysis:computeAll", 0);

} // end computeAll

//------------------------------------------------------------------------------
// Recursively computes an operator node and all its upstream dependencies,
// using timestamp-based cache validation to prevent redundant computation.
// If the node's output is already valid for this timestamp (cache hit),
// returns immediately. Otherwise, recursively computes all upstream
// dependencies first, then calls the node's compute() method.
void Analysis::computeRecursive(OperatorNode *Node,
                                const TimeInstant &TimeStamp) {

   // Check cache: if already computed for this timestamp, return immediately
   if (Node->Op->isCacheValid(TimeStamp)) {
      return;
   }

   // Recursively compute all upstream dependencies first (depth-first)
   for (auto *Upstream : Node->Upstreams) {
      computeRecursive(Upstream, TimeStamp);
   }

   // Now all inputs are fresh - compute this operator
   std::string InstanceName = Node->Op->getName();
   Pacer::start("Analysis:" + InstanceName + ":compute", 1); 
   Node->Op->compute(TimeStamp);
   Pacer::stop("Analysis:" + InstanceName + ":compute", 1); 

} // end computeRecursive

//------------------------------------------------------------------------------
// Calls initialize() on all operators after the dependency graph is complete
// and all Fields exist. This allows operators to store pointers to mesh,
// environment, and other resources they need during compute() calls.
void Analysis::initializeAllOps() {

   // Initialize each operator with context pointers
   for (auto &OpNode : OpNodes) {
      OpNode->Op->initialize(Env, Mesh, VCoord, makeOpConfig());
   }

} // end initializeAllOps

//------------------------------------------------------------------------------
// Returns a reference to the ModelClock pointer. This is used by
// AnalysisGroup instances during stream creation to access the clock
// for alarm initialization.
Clock *&Analysis::getModelClock() { return ModelClock; } // end getModelClock

//------------------------------------------------------------------------------
// Returns a vector of non-owning pointers to all registered operator nodes.
// Used primarily for testing and validation of the dependency graph.
const std::vector<OperatorNode *> Analysis::getOpNodes() {
   std::vector<OperatorNode *> OpPtrs;
   for (auto &Node : OpNodes) {
      OpPtrs.push_back(Node.get());
   }
   return OpPtrs;
} // end getOpNodes

//------------------------------------------------------------------------------
// Checks whether an operator node with the given full instance name already
// exists in the OpNodes vector. Returns true if found, false otherwise.
// Used during chain parsing to avoid creating duplicate operators.
bool Analysis::OpNodeExists(const std::string &FullOpName) {

   for (const auto &Node : OpNodes) {
      if (FullOpName == Node->Op->getName()) {
         return true;
      }
   }

   return false;

} // end OpNodeExists

//------------------------------------------------------------------------------
// Retrieves the default Analysis instance
Analysis *Analysis::getDefault() { return DefAnalysis; } // end getDefault

// Removes all defined Analysis instances and cleans up static resources
void Analysis::finalize() {
   DefAnalysis = nullptr;
   AllAnalysisObjects.clear();
} // end finalize

//------------------------------------------------------------------------------
// Destructor - deallocates all memory
Analysis::~Analysis() {}

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
