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
   CHECK_ERROR_ABORT(Err1, "Analysis: Filename not found in Config");

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
   // Map: StreamName -> list of operator instance names
   std::map<std::string, std::vector<std::string>> StreamToOpNames;

   for (const auto &Info : OpChainInfos) {
      std::string StreamName;

      // Construct stream name based on frequency and type
      if (Info.IsTimeReduction) {
         StreamName = GroupName + "_" + Info.FreqStr + "TimeStats";
      } else {
         StreamName = GroupName + "_" + Info.FreqStr + "Samples";
      }

      // Add this operator chain to the appropriate stream
      StreamToOpNames[StreamName].push_back(Info.ChainStr);
   }

   // Create IOStream objects and associate operator nodes
   for (const auto &[StreamName, OpNames] : StreamToOpNames) {

      // Determine stream type from name suffix
      bool IsTimeReduction =
          (StreamName.find("TimeStats") != std::string::npos);

      // Extract period string from stream name
      // e.g., "GlobalStats_1dayTimeStats" -> "1day"
      size_t UnderscorePos         = StreamName.find_last_of("_");
      std::string PeriodWithSuffix = StreamName.substr(UnderscorePos + 1);
      std::string PeriodStr;
      if (IsTimeReduction) {
         PeriodStr =
             PeriodWithSuffix.substr(0, PeriodWithSuffix.find("TimeStats"));
      } else {
         PeriodStr =
             PeriodWithSuffix.substr(0, PeriodWithSuffix.find("Samples"));
      }

      // Parse period string into numeric frequency and time units
      std::vector<std::string> ParsedStr = parseFreqStr(PeriodStr);
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
          FilenamePrefix + "_" + PeriodStr +
          (IsTimeReduction ? "TimeStats" : "Samples") + FilenameTemplate;
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
