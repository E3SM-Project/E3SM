//===- analysis/AnalysisGroup.cpp - AnalysisGroup implementation -*- C++-*-===//
//
// Implementation of AnalysisGroup base class methods. The key method is
// createAnalysisGroupStreams(), which groups operator chains by their output
// characteristics (period and type), validates temporal reduction periods
// against the restart interval, and creates IOStream objects with appropriate
// configurations for writing analysis output.
//
//===----------------------------------------------------------------------===//

#include "AnalysisGroup.h"

namespace OMEGA {

//------------------------------------------------------------------------------
// Returns the name of this AnalysisGroup instance
std::string AnalysisGroup::getName() { return GroupName; } // end getName

//------------------------------------------------------------------------------
void AnalysisGroup::parseTemporalPeriods(Config &AnalysisGroupOptions) {
   Error Err1;
   Error Err2;

   // Read optional temporal reduction periods (e.g., "1day", "1month")
   Err1 = AnalysisGroupOptions.get("ReductionPeriod", ReductionPeriodList);

   // Read optional instantaneous output intervals (e.g., "6hour")
   Err2 = AnalysisGroupOptions.get("SnapshotPeriod", SnapshotPeriodList);

   // Validate that at least one temporal specification is provided
   if (Err1.isFail() && Err2.isFail()) {
      ABORT_ERROR("{}: Error reading both ReductionPeriod and "
                  "SnapshotPeriod from Config, at least one must be present",
                  GroupName);
   }
} // end parseTemporalPeriods

//------------------------------------------------------------------------------
void AnalysisGroup::buildTemporalChains(
    const std::vector<std::string> &ChainStems, Config &AnalysisGroupOptions,
    Analysis *AnalysisManager) {

   // Parse temporal periods from config
   parseTemporalPeriods(AnalysisGroupOptions);

   // Build temporal reduction chains for each stem
   for (const auto &StemStr : ChainStems) {

      // Create temporal reduction chains: Stem -> TimeMean<Period>
      for (const auto &ReductionPeriod : ReductionPeriodList) {
         // Build chain string (e.g., "Temperature_SpatialMean_TimeMean1day")
         std::string ChainStr = StemStr + "_TimeMean" + ReductionPeriod;

         // Store metadata for stream creation
         OpChainInfos.push_back(OpChainInfo{ChainStr, ReductionPeriod, true});

         // Parse chain and instantiate operators
         AnalysisManager->parseChainAndBuildOps(ChainStr);
      }

      // Create instantaneous snapshot chains (if requested)
      if (!SnapshotPeriodList.empty()) {
         // Parse stem chain if not already built (operators may exist from
         // reduction chains above, parseChainAndBuildOps handles duplicates)
         AnalysisManager->parseChainAndBuildOps(StemStr);

         // Store metadata for each snapshot frequency
         for (const auto &SnapshotPeriod : SnapshotPeriodList) {
            OpChainInfos.push_back(OpChainInfo{StemStr, SnapshotPeriod, false});
         }
      }
   }
} // end buildTemporalChains

//------------------------------------------------------------------------------
// Groups operator chains by output period and type (temporal reduction vs.
// discrete sampling), validates temporal reduction periods against the restart
// interval, and creates IOStream objects for the group's output. Each stream
// is associated with the appropriate operator nodes, and the operators' output
// fields are added to the stream contents.
void AnalysisGroup::createAnalysisGroupStreams(const std::string &GroupName,
                                               Config &AnalysisGroupOptions,
                                               Analysis *AnalysisManager) {

   Error Err1;
   Error Err2;

   // Read optional Stream config node from AnalysisGroup configuration
   Config AnalysisStreamOptions("Stream");
   Err1 = AnalysisGroupOptions.get(AnalysisStreamOptions);
   std::map<std::string, std::string> ParamOverrides;

   // If Stream config exists, extract parameter overrides
   if (Err1.isSuccess()) {
      for (Config::Iter It = AnalysisStreamOptions.begin();
           It != AnalysisStreamOptions.end(); ++It) {
         std::string Key = It->first.as<std::string>();
         std::string Value;
         AnalysisStreamOptions.get(Key, Value);
         ParamOverrides[Key] = Value;
      }
   }

   // Create StreamParams with defaults and apply any overrides
   StreamParams StreamCfg;
   StreamCfg.apply(ParamOverrides);

   // Parse filename into prefix and timestamp template
   // e.g., "analysis.$Y.$M" -> prefix="analysis", template=".$Y.$M"
   std::string FilenameStr;
   Err1 = AnalysisGroupOptions.get("Filename", FilenameStr);
   CHECK_ERROR_ABORT(Err1,
                     "AnalysisGroup: Filename not found in Config for group {}",
                     GroupName);

   std::string FilenamePrefix;
   std::string FilenameTemplate;
   size_t Pos = FilenameStr.find("$");
   if (Pos != std::string::npos) {
      if (Pos == 0) {
         FilenamePrefix.clear();
         FilenameTemplate = FilenameStr;
      } else {
         // Split at '$'. If a separator (e.g., '.' or '_') immediately precedes
         // the template, keep it in the template for backwards compatibility.
         FilenamePrefix = FilenameStr.substr(0, Pos);
         if (!FilenamePrefix.empty() &&
             (FilenamePrefix.back() == '.' || FilenamePrefix.back() == '_')) {
            FilenamePrefix.pop_back();
            FilenameTemplate = FilenameStr.substr(Pos - 1);
         } else {
            FilenameTemplate = FilenameStr.substr(Pos);
         }
      }
   } else {
      // No timestamp template - use entire string as prefix
      FilenamePrefix = FilenameStr;
   }

   // Group operator chains by their output stream characteristics
   // Chains with the same frequency and type go to the same stream
   // Map: StreamName -> (OpNames, IsTimeReduction, FreqStr)
   std::map<std::string, StreamInfo> StreamInfos;

   for (const auto &OpInfo : OpChainInfos) {
      std::string StreamName;

      // Construct stream name based on frequency and type
      if (OpInfo.IsTimeReduction) {
         StreamName = GroupName + "_" + OpInfo.FreqStr + "TimeStats";
      } else {
         StreamName = GroupName + "_" + OpInfo.FreqStr + "Instants";
      }

      // Add this operator chain to the appropriate stream
      // If this is the first operator for this stream, initialize metadata
      if (StreamInfos.find(StreamName) == StreamInfos.end()) {
         StreamInfos[StreamName] = {
             {OpInfo.ChainStr},     // OpNames vector with first entry
             OpInfo.FreqStr,        // Store the frequency string
             OpInfo.IsTimeReduction // Store the flag
         };
      } else {
         // Stream already exists, just add this operator
         StreamInfos[StreamName].OpNames.push_back(OpInfo.ChainStr);
      }
   }

   // Create IOStream objects and associate operator nodes
   for (const auto &[StreamName, Metadata] : StreamInfos) {

      // Extract metadata (no need to parse stream name!)
      bool IsTimeReduction       = Metadata.IsTimeReduction;
      const std::string &FreqStr = Metadata.FreqStr;
      const auto &OpNames        = Metadata.OpNames;

      // Parse frequency string into numeric frequency and time units
      std::vector<std::string> ParsedStr = parseFreqStr(FreqStr);
      I4 Freq                            = std::stoi(ParsedStr[0]);
      TimeUnits FreqUnits                = TimeUnitsFromString(ParsedStr[1]);
      TimeInterval PeriodInterval(Freq, FreqUnits);

      // Validate temporal reduction periods against restart interval
      // This ensures proper checkpoint/restart behavior for time-averaged
      // fields
      if (IsTimeReduction) {
         auto RestartAlarm = IOStream::getAlarm("RestartWrite");
         bool IsDivisible =
             RestartAlarm->getInterval()->isDivisibleBy(PeriodInterval);
         if (!IsDivisible) {
            ABORT_ERROR(
                "Analysis: The RestartWrite interval is not divisible by "
                "the averaging period, {} {}. Currently, temporal averaging "
                "is only available over intervals where "
                "RestartPeriod % PeriodInterval == 0",
                ParsedStr[0], ParsedStr[1]);
         }
      }

      // Build filename for this stream (includes period and type in name)
      StreamCfg.Params["Filename"] =
          FilenamePrefix + "_" + FreqStr +
          (IsTimeReduction ? "TimeStats" : "Instants") + FilenameTemplate;
      StreamCfg.Params["Freq"]      = ParsedStr[0];
      StreamCfg.Params["FreqUnits"] = ParsedStr[1];

      // Create the IOStream with configured parameters
      auto NewStreamCfg = StreamCfg.toConfig();
      auto RefClock     = AnalysisManager->getModelClock();
      IOStream::create(StreamName, NewStreamCfg, RefClock);

      // Retrieve the newly created stream
      auto Stream = IOStream::get(StreamName);

      // Associate operator nodes with this stream and populate stream contents
      auto OpNodes = AnalysisManager->getOpNodes();
      for (auto *Node : OpNodes) {
         std::string OpInstanceName = Node->Op->getName();

         // Check if this operator belongs to this stream
         if (std::find(OpNames.begin(), OpNames.end(), OpInstanceName) !=
             OpNames.end()) {
            // Add stream name to operator's list (used for alarm association)
            Node->StreamNames.push_back(StreamName);

            // Add all of the operator's output fields to the stream
            for (const auto &FieldName : Node->Op->getOutputFieldNames()) {
               Stream->addField(FieldName);
            }
         }
      }
   }

} // end createAnalysisGroupStreams

} // end namespace OMEGA
