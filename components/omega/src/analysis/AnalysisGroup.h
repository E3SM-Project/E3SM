#ifndef OMEGA_ANALYSISGROUP_H
#define OMEGA_ANALYSISGROUP_H

//===-- analysis/AnalysisGroup.h - AnalysisGroup base class -----*- C++ -*-===//
//
/// \file
/// \brief Defines the AnalysisGroup base class for bundled analysis groups
///
/// AnalysisGroup is the abstract base class for bundled analysis groups that
/// encapsulate common analysis patterns. In the initial implementation (v1),
/// concrete derived classes (e.g., GlobalStats) handle config parsing, operator
/// construction, and stream creation for pre-defined analysis group types.
/// Future versions will support user-defined custom groups specified entirely
/// in config using composable operator chains.
///
/// Each AnalysisGroup reads its configuration, constructs operator chains
/// by calling Analysis::parseChainAndBuildOps(), and creates associated
/// IOStream objects for output. The base class provides helper structures
/// (OpChainInfo, StreamParams) and methods (createAnalysisGroupStreams) that
/// facilitate grouping operator chains by output frequency and type, validating
/// temporal reduction periods, and creating appropriately configured streams.
///
//===----------------------------------------------------------------------===//

#include "Analysis.h"
#include "Config.h"
#include "IOStream.h"
#include <sstream>
#include <string>

namespace OMEGA {

/// AnalysisGroup is the abstract base class for bundled analysis groups.
/// Concrete derived classes (e.g., GlobalStats) encapsulate the configuration
/// parsing, operator construction, and stream creation logic for named analysis
/// group types. The base class provides utilities for grouping operator chains
/// by their output characteristics and creating IOStream objects with
/// appropriate configurations. Future versions will support user-defined custom
/// groups where users specify composable operator chains directly in the config
/// file.
class AnalysisGroup {

 public:
   /// Virtual destructor allows polymorphic deletion of derived classes
   virtual ~AnalysisGroup() = default;

   /// Returns the name of this AnalysisGroup instance
   std::string getName();

   /// Groups operator chains by output period and type (time reduction vs.
   /// instantaneous output intervals), validates temporal reduction periods
   /// against the restart interval, and creates IOStream objects for the
   /// group's output. Called by derived class constructors after all operator
   /// chains have been registered with the Analysis orchestrator.
   void createAnalysisGroupStreams(
       const std::string &GroupName, ///< [in] name of this group
       Config &AnalysisGroupOptions, ///< [in] group configuration options
       Analysis *AnalysisManager     ///< [in] Analysis orchestrator instance
   );

 protected:
   /// Metadata about a single operator chain within this AnalysisGroup.
   /// Stores the chain string (operator instance name), the output frequency,
   /// and whether the chain outputs temporal reductions or instantaneous
   /// snapshots.
   struct OpChainInfo {
      std::string ChainStr; ///< Operator instance name (output Field name)
      std::string FreqStr;  ///< Period/frequency string (e.g., "1Day", "6Hour")
      bool IsTimeReduction; ///< true for temporal reduction; false for
                            ///< instantaneous output
   };

   /// Aggregated metadata for a single output stream within this AnalysisGroup.
   /// Groups multiple operator chains that share the same output frequency and
   /// type (temporal reduction vs. instantaneous) into a single IOStream.
   /// Created by grouping OpChainInfo objects in createAnalysisGroupStreams().
   struct StreamInfo {
      std::vector<std::string> OpNames; ///< Operator instance names to include
      std::string FreqStr;  ///< Period/frequency string (e.g., "1Day", "6Hour")
      bool IsTimeReduction; ///< true for temporal reduction; false for
                            ///< instantaneous output
   };

   /// Sets the IOName metadata on an output Field so IOStream writes it under
   /// a legible name rather than the full operator-chain string.
   /// No-op if the field does not exist.
   static void setOutputIOName(const std::string &InternalFieldName,
                               const std::string &IOName);

   /// Template for constructing IOStream configurations for this group's
   /// output. Provides default values for all IOStream creation parameters.
   /// Derived classes can override defaults using group-specific config options
   /// via the apply() method. The toConfig() method converts the parameters to
   /// a Config object suitable for IOStream::create().
   struct StreamParams {
      /// Constructor initializes all IOStream parameters with default values.
      /// Empty string values indicate parameters that must be set by derived
      /// classes or will be omitted from the final stream configuration.
      StreamParams()
          : Params{
                {"UsePointerFile", "false"},
                {"PointerFilename", ""},
                {"Filename", ""},
                {"Mode", "write"},
                {"IfExists", "append"},
                {"Precision", "double"},
                {"Freq", ""},
                {"FreqUnits", ""},
                {"FileFreq", "9999"},
                {"FileFreqUnits", "Years"},
                {"UseStartEnd", "false"},
                {"StartTime", ""},
                {"EndTime", ""},
            } {}

      /// Applies group-specific configuration overrides to the default
      /// parameter values. Only known parameters can be overridden; unknown
      /// keys trigger an error.
      void apply(const std::map<std::string, std::string>
                     &Overrides ///< [in] parameter overrides
      ) {
         for (const auto &[Key, Value] : Overrides) {
            auto It = Params.find(Key);

            // Validate that the key exists in our parameter map
            if (It == Params.end()) {
               ABORT_ERROR("Analysis: Unknown Stream config parameter, {}",
                           Key);
            }

            It->second = Value;
         }
      }

      /// Converts the parameter map to a Config object suitable for passing
      /// to IOStream::create(). Only parameters with non-empty values are
      /// included in the Config. Contents are left empty and must be added
      /// after stream creation.
      Config toConfig() const {
         Config Cfg;

         // Add only parameters with non-empty values
         for (const auto &[Key, Value] : Params) {
            if (!Value.empty()) {
               Cfg.add(Key, Value);
            }
         }

         // Contents are populated after stream creation via IOStream::addField
         std::vector<std::string> EmptyStrVec{""};
         Cfg.add("Contents", EmptyStrVec);

         return Cfg;
      }

      /// Map of all IOStream configuration options for this group's output
      std::map<std::string, std::string> Params;
   };

   /// Reads ReductionPeriod and SnapshotPeriod from config, builds temporal
   /// operator chains by appending time operators to the provided stems,
   /// calls parseChainAndBuildOps for each chain, and populates OpChainInfos
   /// with metadata for stream creation.
   ///
   /// If ChainConfigs is provided and non-empty, it must have the same size
   /// as ChainStems, with each Config containing operator-specific parameters
   /// for the corresponding chain. If ChainConfigs is empty, chains are built
   /// without additional configuration (default behavior).
   void buildTemporalChains(
       const std::vector<std::string>
           &ChainStems,              ///< [in] operator chain stems
       Config &AnalysisGroupOptions, ///< [in] group configuration
       Analysis *AnalysisManager,    ///< [in] analysis manager
       const std::vector<Config> &ChainConfigs = {}
       ///< [in] optional per-chain configs
   );

   /// Reads ReductionPeriod and SnapshotPeriod from config and validates
   /// that at least one is present
   void parseTemporalPeriods(Config &AnalysisGroupOptions);

   std::string GroupName; ///< Name of this AnalysisGroup instance

   std::vector<std::string> ReductionPeriodList; ///< Temporal averaging periods
   std::vector<std::string>
       SnapshotPeriodList; ///< Instantaneous sample periods

   /// Metadata for all operator chains in this group (output field name,
   /// frequency, and type). Populated by derived class constructors.
   std::vector<OpChainInfo> OpChainInfos;

   /// Full operator chain strings for all chains in this group. This vector
   /// stores the complete underscore-delimited chain strings for reference.
   std::vector<std::string> OpChainStrings;

}; // end class AnalysisGroup

} // end namespace OMEGA

#endif
