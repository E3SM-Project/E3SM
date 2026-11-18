#ifndef OMEGA_TIMEMEANOP_H
#define OMEGA_TIMEMEANOP_H

//===-- analysis/operators/TimeMeanOp.h - TimeMeanOp ------------*- C++ -*-===//
//
/// \file
/// \brief Defines the TimeMeanOp operator for computing time-averaged mean
///
/// TimeMeanOp computes the time-averaged mean of a field over a specified
/// period (e.g., "1day", "1month"). Unlike spatial operators that reduce fields
/// to scalars, TimeMeanOp preserves the full spatial structure of the input
/// field while averaging temporally. The operator accumulates input values at
/// each timestep and divides by the number of accumulations when the period
/// alarm rings.
///
/// The operator is templated on the Kokkos array type (ArrayT) of the input
/// field, supporting any rank and memory layout. The output Field has the same
/// dimensions and layout as the input Field. Accumulation occurs in-place in
/// the output array using element-wise addition.
///
/// The temporal averaging period is specified via the "Period" configuration
/// option and is managed by an Alarm that signals when to finalize the average.
/// When the alarm rings, the accumulated sum is divided by the count, the
/// result is written to the output, and accumulation restarts for the next
/// period. The operator maintains state (accumulation count and period flag)
/// across timesteps.
///
//===----------------------------------------------------------------------===//

#include "AnalysisOperator.h"
#include "TimeStepper.h"

#include <algorithm>
#include <cctype>
#include <limits>
#include <memory>
#include <string>

namespace OMEGA {

/// TimeMeanOp computes the time-averaged mean of a field over a specified
/// period. The operator accumulates input values at each timestep, and when
/// the period alarm rings, divides by the accumulation count to produce the
/// mean. Output Field has the same dimensions and layout as input but always
/// uses Real type for accuracy. Maintains state across timesteps for
/// accumulation.
template <typename ArrayT> class TimeMeanOp : public AnalysisOperator {
 public:
   /// Scalar type extracted from the input array type
   using ScalarT = typename ArrayT::non_const_value_type;

   /// Output array type - same rank as input but always Real precision
   using OutputArrayT = typename std::conditional<
       ArrayT::rank == 1, Array1D_t<Real>,
       typename std::conditional<
           ArrayT::rank == 2, Array2D_t<Real>,
           typename std::conditional<ArrayT::rank == 3, Array3D_t<Real>,
                                     Array4D_t<Real>>::type>::type>::type;

   /// Constructs a TimeMeanOp operator. Reads averaging period from config,
   /// creates output Field matching input dimensions and layout, allocates
   /// output data array for accumulation, and registers the output Field in
   /// the Field registry. Initializes accumulation state. The output Field
   /// name is constructed as InputName + "_TimeMean" + Period (e.g.,
   /// "Temperature_TimeMean1day").
   TimeMeanOp(const std::vector<std::string>
                  &UpstreamNames, ///< [in] input field names
              Config Options      ///< [in] operator config
              )
       : AnalysisOperator("TimeMean") {

      // Store input field names
      InputNames = UpstreamNames;

      // Read averaging period from configuration (e.g., "1day", "1month")
      std::string AvgPeriod;
      Options.get("Period", AvgPeriod);
      PeriodStr = AvgPeriod;

      // Read the optional sampling interval. Empty => sample every timestep
      // (the default, bit-for-bit with prior behavior). The string is parsed
      // and normalized to a TimeInterval later, in setPeriodAlarm(), once the
      // reduction period and model timestep are both available.
      SamplingIntervalStr.clear();
      Options.get("SamplingInterval", SamplingIntervalStr);

      // Construct output field name and set instance name
      std::string OutputFieldName = InputNames[0] + "_TimeMean" + AvgPeriod;
      OutputNames                 = {OutputFieldName};
      InstanceName                = OutputFieldName;

      // Retrieve input Field to extract metadata
      auto InputField = Field::get(InputNames[0]);
      auto InputData  = InputField->template getDataArray<ArrayT>();

      // Get dimension info (number of dims and dimension names)
      auto NDims = InputField->getNumDims();
      std::vector<std::string> DimNames;
      InputField->getDimNames(DimNames);

      // Fetch the input description and valid range
      std::string InputDescr;
      ScalarT InputValidMin, InputValidMax;
      InputField->getMetadata("long_name", InputDescr);
      InputField->getMetadata("valid_min", InputValidMin);
      InputField->getMetadata("valid_max", InputValidMax);

      // Register output Field with same dimensions as input but Real type.
      // A time mean has the units and standard name of the field it
      // averages; its cell_methods record the reduction after any the input
      // already carries (e.g. "area: mean time: mean").
      auto Meta = inheritMetadata(InputNames[0], "time: mean");
      auto OutputField =
          createOutputField(OutputNames[0],
                            "Time average of " + InputDescr,  // Description
                            Meta,                             // CF metadata
                            static_cast<Real>(InputValidMin), // Min valid
                            static_cast<Real>(InputValidMax), // Max valid
                            NDims,                            // Rank
                            DimNames                          // Dim names
          );

      // Store array size for parallel iteration
      ArraySize = static_cast<I4>(InputData.size());

      // Allocate output data array matching input layout but with Real type
      OutputData = OutputArrayT(OutputNames[0] + "_out", InputData.layout());

      // Attach output data array to Field
      OutputField->template attachData<OutputArrayT>(OutputData);

      // Initialize accumulation state
      NumAccum    = 0;
      PeriodAlarm = nullptr;
      IsNewPeriod = true;

      // Default to the flat (non-columnar) path until initialize() detects a
      // recognized horizontal mesh dimension.
      IsMeshDimension = false;
      NHorizOwned     = 0;

      // Create the restart checkpoint Fields now that the output Field and
      // OutputData array exist. This registers them in the "Restart"
      // FieldGroup during construction, before IOStream::validateAll() and
      // the restart read, which both occur later in ocnInit.
      createRestartFields();

   } // end constructor

   /// Detaches the operator-owned sample alarm from the model clock, which
   /// would otherwise keep a pointer to it after the operator is destroyed.
   /// The clock must outlive the operator, as it does when Analysis is
   /// finalized before the time stepper is cleared.
   ~TimeMeanOp() override { detachSampleAlarm(); }

   /// Initializes the operator after all Fields exist. Determines the index
   /// space (cells, edges, or vertices) from the input field's horizontal
   /// dimension name and, for rank >= 2 mesh fields, stores the appropriate
   /// MinLayer/MaxLayer arrays from VertCoord. Fields defined on vertical
   /// interfaces add 1 to KMax to include the interface on the bottommost
   /// active layer. TimeMeanOp does NOT abort on an unrecognized index space
   /// or final dimension: it simply falls back to the flat parallel dispatch
   /// path.
   void initialize(const MachEnv *Env, const HorzMesh *InMesh,
                   const VertCoord *InVCoord, Config Options) override {

      AnalysisOperator::initialize(Env, InMesh, InVCoord, Options);

      constexpr I4 InputRank = ArrayT::rank;

      auto InputField = Field::get(InputNames[0]);
      std::vector<std::string> DimNames;
      InputField->getDimNames(DimNames);

      // For rank 1 the sole dimension is the horizontal index space; for rank
      // >= 2 the horizontal index space is the 2nd-to-last dimension.
      std::string IndexSpaceName;
      if constexpr (InputRank == 1) {
         IndexSpaceName = DimNames[0];
      } else {
         IndexSpaceName                     = DimNames[InputRank - 2];
         const std::string &VerticalDimName = DimNames[InputRank - 1];
         LayerBoundOffset = (VerticalDimName == "NVertLayersP1") ? 1 : 0;
      }

      IsMeshDimension =
          (IndexSpaceName == "NCells" || IndexSpaceName == "NEdges" ||
           IndexSpaceName == "NVertices");

      // Non-mesh output: leave IsMeshDimension false and use the flat path.
      if (!IsMeshDimension) {
         return;
      }

      // Store owned horizontal count and, for rank >= 2, the active-layer
      // bounds appropriate to the horizontal location.
      if (IndexSpaceName == "NCells") {
         NHorizOwned = Mesh->NCellsOwned;
         if constexpr (InputRank > 1) {
            MinLayer = VCoord->MinLayerCell;
            MaxLayer = VCoord->MaxLayerCell;
         }
      } else if (IndexSpaceName == "NEdges") {
         NHorizOwned = Mesh->NEdgesOwned;
         if constexpr (InputRank > 1) {
            MinLayer = VCoord->MinLayerEdgeBot;
            MaxLayer = VCoord->MaxLayerEdgeTop;
         }
      } else { // NVertices
         NHorizOwned = Mesh->NVerticesOwned;
         if constexpr (InputRank > 1) {
            MinLayer = VCoord->MinLayerVertexBot;
            MaxLayer = VCoord->MaxLayerVertexTop;
         }
      }

   } // end initialize

   /// Sets the period alarm that signals when to finalize the time average.
   /// Called by Analysis during initialization to provide the alarm associated
   /// with the averaging period specified in the constructor config. If a
   /// SamplingInterval was configured, this is also where the three accepted
   /// forms are normalized to a single TimeInterval, the period-divisibility
   /// constraint is enforced, and the operator-owned recurring SampleAlarm is
   /// built and attached to the model clock.
   void setPeriodAlarm(Alarm *AlarmIn ///< [in] alarm for averaging period
                       ) override {
      PeriodAlarm = AlarmIn;

      // Make a case-insensitive working copy of the config string to check
      // for variations of "dt" and "timestep". If match, sampling cadence is
      // every timestep and nothing more to do.
      std::string SamplingIntervalLower = SamplingIntervalStr;
      std::transform(SamplingIntervalLower.begin(), SamplingIntervalLower.end(),
                     SamplingIntervalLower.begin(),
                     [](unsigned char C) { return std::tolower(C); });
      if (SamplingIntervalLower == "dt" || SamplingIntervalLower == "timestep")
         return;

      // Empty string is equivalent to sampling every timestep.
      if (SamplingIntervalStr.empty())
         return;

      // The reduction period comes from the period alarm; the model timestep
      // and the absolute simulation start come from the default clock.
      const TimeInterval PeriodInterval = *(PeriodAlarm->getInterval());

      auto *Stepper = TimeStepper::getDefault();
      OMEGA_REQUIRE(Stepper,
                    "TimeMeanOp: null default TimeStepper in setPeriodAlarm");
      Clock *ModelClock = Stepper->getClock();
      OMEGA_REQUIRE(ModelClock,
                    "TimeMeanOp: null model clock in setPeriodAlarm");
      const TimeInterval ModelTimeStep = ModelClock->getTimeStep();

      // Anchor the sample schedule to the fixed absolute clock start rather
      // than to the period alarm's runtime-mutated previous ring time. The
      // clock start (Clock::getStartTime) is identical in straight-through and
      // restart runs, so the sample-boundary sequence it generates is a
      // deterministic walk from a fixed origin.
      const TimeInstant ClockStart = ModelClock->getStartTime();

      // Normalize the configured string to a single TimeInterval. The suffix
      // selects the form:
      //   "<n>perPeriod" => PeriodInterval / n (n samples per period)
      //   "<n>steps"     => n * ModelTimeStep  (every n timesteps)
      //   otherwise      => parsed directly as a TimeInterval string
      const std::string PerPeriodSuffix = "perPeriod";
      const std::string StepsSuffix     = "steps";

      auto endsWith = [](const std::string &Str, const std::string &Suffix) {
         return Str.size() >= Suffix.size() &&
                Str.compare(Str.size() - Suffix.size(), Suffix.size(),
                            Suffix) == 0;
      };

      if (endsWith(SamplingIntervalStr, PerPeriodSuffix)) {
         std::string NumStr = SamplingIntervalStr.substr(
             0, SamplingIntervalStr.size() - PerPeriodSuffix.size());
         I4 NumSamples = std::stoi(NumStr);
         if (NumSamples <= 0)
            ABORT_ERROR("TimeMeanOp: SamplingInterval '{}' must request a "
                        "positive number of samples per period",
                        SamplingIntervalStr);
         SamplingInterval = dividePeriod(PeriodInterval, NumSamples);
      } else if (endsWith(SamplingIntervalStr, StepsSuffix)) {
         std::string NumStr = SamplingIntervalStr.substr(
             0, SamplingIntervalStr.size() - StepsSuffix.size());
         I4 NumSteps = std::stoi(NumStr);
         if (NumSteps <= 0)
            ABORT_ERROR("TimeMeanOp: SamplingInterval '{}' must request a "
                        "positive number of timesteps",
                        SamplingIntervalStr);
         SamplingInterval = NumSteps * ModelTimeStep;
      } else {
         // Raw frequency-string form (e.g. "600seconds", "6hour"). Parse it
         // the same way AnalysisGroup does, via the numeric (value, TimeUnits)
         // TimeInterval constructor, so "<n>unit" strings are accepted.
         std::vector<std::string> Parsed = parseFreqStr(SamplingIntervalStr);
         I4 Freq                         = std::stoi(Parsed[0]);
         TimeUnits FreqUnits             = TimeUnitsFromString(Parsed[1]);
         SamplingInterval                = TimeInterval(Freq, FreqUnits);
      }

      // Enforce the single period-divisibility constraint: the sampling
      // interval must divide the reduction period evenly so the unweighted
      // count-mean stays an unbiased estimate of the time-average.
      if (!PeriodInterval.isDivisibleBy(SamplingInterval)) {
         ABORT_ERROR("TimeMeanOp: SamplingInterval '{}' must divide the "
                     "averaging period evenly (PeriodInterval % "
                     "SamplingInterval == 0).",
                     SamplingIntervalStr);
      }

      HasSamplingInterval = true;

      // Build and own a recurring sample alarm anchored to the absolute clock
      // start, so its boundary sequence (ClockStart + n*SamplingInterval) is a
      // deterministic walk reproduced identically on fresh starts and restarts.
      // An alarm from an earlier call is detached before it is replaced.
      detachSampleAlarm();
      OwnedSampleAlarm = std::make_unique<Alarm>(InstanceName + "_SampleAlarm",
                                                 SamplingInterval, ClockStart);
      ModelClock->attachAlarm(OwnedSampleAlarm.get());
      SampleClock = ModelClock;
      SampleAlarm = OwnedSampleAlarm.get();

      // Seed the freshly constructed alarm to the current clock time. Per
      // Alarm::reset, a periodic alarm walks its boundaries forward to the
      // first interval after the input time. On a restart, Analysis is built
      // before the restart time is read, and it is Clock::setCurrentTime that
      // then resets this alarm, like every other periodic alarm attached to
      // the clock, to the boundary a straight-through run would hold at the
      // restart time -- without introducing any new restart state.
      OwnedSampleAlarm->reset(ModelClock->getCurrentTime());
   }

   /// TimeMeanOp maintains mid-period accumulation state (the partial sum in
   /// OutputData and the sample count NumAccum) that must survive restarts.
   bool needsRestart() const override { return true; }

   /// Creates the two Fields that checkpoint the operator's mid-period state
   /// and registers them in the "Restart" FieldGroup so they are written to
   /// and read from the restart file by the existing restart streams:
   ///   1. <OutputName>_AccumSum shares the operator's OutputData array, so
   ///      the running partial sum is persisted with no extra copy. (Between
   ///      finalize() and the next accumulation OutputData already holds the
   ///      finished mean; that is harmless because NumAccum is normalized to
   ///      0 on write and the array is overwritten on the next new period.)
   ///   2. <OutputName>_NumAccum is a one-element I4 host array holding the
   ///      sample count.
   /// Both are marked OptionalRead so restart files written before this
   /// feature existed (which lack these variables) read cleanly, leaving
   /// NumAccum at 0 (a fresh period, matching the previous behavior).
   void createRestartFields() {

      // Ensure the Restart group exists (other modules may create it first)
      if (!FieldGroup::exists("Restart")) {
         FieldGroup::create("Restart");
      }

      const std::string AccumSumName = OutputNames[0] + "_AccumSum";
      const std::string NumAccumName = OutputNames[0] + "_NumAccum";

      // --- Partial sum Field: shares the existing OutputData array ---
      auto OutputField = Field::get(OutputNames[0]);
      auto NDims       = OutputField->getNumDims();
      std::vector<std::string> DimNames;
      OutputField->getDimNames(DimNames);

      auto AccumSumField = Field::create(
          AccumSumName, "Restart checkpoint of time-mean partial sum", "", "",
          std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::max(),
          NDims, DimNames);
      AccumSumField->setOptionalRead(true);
      FieldGroup::addFieldToGroup(AccumSumName, "Restart");
      AccumSumField->template attachData<OutputArrayT>(OutputData);

      // --- Sample-count Field: a one-element non-distributed I4 array ---
      const std::string CountDimName = "NTimeMeanCount";
      if (!Dimension::exists(CountDimName)) {
         Dimension::create(CountDimName, 1);
      }

      NumAccumData = HostArray1DI4(NumAccumName + "_data", 1);

      auto NumAccumField =
          Field::create(NumAccumName, "Restart checkpoint of time-mean count",
                        "", "", 0, std::numeric_limits<I4>::max(), 1,
                        std::vector<std::string>{CountDimName});
      NumAccumField->setOptionalRead(true);
      FieldGroup::addFieldToGroup(NumAccumName, "Restart");
      NumAccumField->template attachData<HostArray1DI4>(NumAccumData);

   } // end createRestartFields

   /// Restores live state from the restart-read count. IsNewPeriod is derived
   /// from the count: a zero count means no partial period is in progress and
   /// accumulation starts fresh (also the default for restart files that lack
   /// the checkpoint variables, since OptionalRead leaves the array at its
   /// fill value of 0 for I4).
   void applyRestartState() override {
      NumAccum    = NumAccumData(0);
      IsNewPeriod = (NumAccum == 0);
   } // end applyRestartState

   /// Computes the time-averaged mean by accumulating input values. On the
   /// first call after alarm rings (IsNewPeriod), copies input to output and
   /// sets count to 1. On subsequent calls, adds input to accumulated sum and
   /// increments count. When period alarm rings, divides accumulated sum by
   /// count to finalize mean and resets state for next period. Updates
   /// timestamp and computed flag.
   void compute(const TimeInstant &TimeStamp ///< [in] current timestamp
                ) override {

      // Retrieve input Field and extract data array
      auto InputField = Field::get(InputNames[0]);
      auto InputData  = InputField->template getDataArray<ArrayT>();

      // Sampling gate: when a SamplingInterval is configured, only accumulate
      // on sample steps (SampleAlarm ringing); otherwise (default) the gate is
      // always open and every timestep is sampled, bit-for-bit with prior
      // behavior. Accumulation is gated, but the period-finalize logic below
      // runs every call so a period boundary is always honored -- even when it
      // lands on a non-sample step.
      if (shouldSample(TimeStamp)) {
         // Accumulate input values into output array
         if (IsNewPeriod) {
            // Start of new averaging period: initialize with first input
            NumAccum = 1;
            accumInit(InputData);
            IsNewPeriod = false;
         } else {
            // Continue accumulation: add input to running sum
            accumAdd(InputData);
            ++NumAccum;
         }
      }

      // Period finalization runs outside the sampling gate so a boundary that
      // falls on a non-sample step is still honored. Only act when a partial
      // period is in progress (NumAccum > 0 / not IsNewPeriod); a boundary that
      // arrives before the first sample of a period has nothing to finalize.
      if (!IsNewPeriod && PeriodAlarm != nullptr && PeriodAlarm->isRinging()) {
         // Finalize: divide accumulated sum by count to get mean
         finalize();

         // Reset state for next averaging period
         IsNewPeriod = true;
      }

      // Update cache validity markers
      LastComputed  = TimeStamp;
      FieldComputed = true;

      // Refresh the persisted sample count now that NumAccum/IsNewPeriod have
      // settled. compute() is the only place NumAccum is mutated, so this keeps
      // the restart checkpoint consistent regardless of restart-write timing.
      // After a finalize (IsNewPeriod true) NumAccum still holds the finished
      // period's count, so it is normalized to 0 to mark no partial period.
      NumAccumData(0) = IsNewPeriod ? 0 : NumAccum;

   } // end compute

   /// Initializes the output accumulator with the first input sample of a new
   /// averaging period. For column-structured mesh fields the write is
   /// restricted to active layers [MinLayer, MaxLayer] using hierarchical
   /// parallelism. Otherwise a flat parallel loop over the entire array is
   /// used.
   void accumInit(const ArrayT &InputData) {
      OMEGA_SCOPE(LocOutputData, OutputData);
      constexpr I4 InputRank = ArrayT::rank;

      if (IsMeshDimension && InputRank >= 2) {
         OMEGA_SCOPE(LocMinLayer, MinLayer);
         OMEGA_SCOPE(LocMaxLayer, MaxLayer);
         OMEGA_SCOPE(LocInputData, InputData);
         const I4 LocLayerBoundOffset = LayerBoundOffset;

         if constexpr (InputRank == 2) {
            parallelForOuter(
                "TimeMeanInit2D", LaunchConfig({NHorizOwned}),
                KOKKOS_LAMBDA(int IHoriz, const TeamMember &Team) {
                   const I4 KMin   = LocMinLayer(IHoriz);
                   const I4 KMax   = LocMaxLayer(IHoriz) + LocLayerBoundOffset;
                   const I4 KRange = vertRange(KMin, KMax);
                   parallelForInner(
                       Team, KRange, INNER_LAMBDA(int KIdx) {
                          const I4 K = KMin + KIdx;
                          LocOutputData(IHoriz, K) =
                              static_cast<Real>(LocInputData(IHoriz, K));
                       });
                });
         } else if constexpr (InputRank == 3) {
            const I4 Dim0 = InputData.extent(0);
            parallelForOuter(
                "TimeMeanInit3D", LaunchConfig({Dim0, NHorizOwned}),
                KOKKOS_LAMBDA(int I0, int IHoriz, const TeamMember &Team) {
                   const I4 KMin   = LocMinLayer(IHoriz);
                   const I4 KMax   = LocMaxLayer(IHoriz) + LocLayerBoundOffset;
                   const I4 KRange = vertRange(KMin, KMax);
                   parallelForInner(
                       Team, KRange, INNER_LAMBDA(int KIdx) {
                          const I4 K = KMin + KIdx;
                          LocOutputData(I0, IHoriz, K) =
                              static_cast<Real>(LocInputData(I0, IHoriz, K));
                       });
                });
         }
      } else {
         // Flat fallback for scalars, rank-1, and non-mesh outputs
         parallelFor(
             {ArraySize}, KOKKOS_LAMBDA(const int FlatIdx) {
                LocOutputData.data()[FlatIdx] =
                    static_cast<Real>(InputData.data()[FlatIdx]);
             });
      }
   } // end accumInit

   /// Adds an input sample to the running accumulator. Uses the same active
   /// layer bounds and dispatch logic as accumInit so that only the layers
   /// initialized at the start of the period are accumulated.
   void accumAdd(const ArrayT &InputData) {
      OMEGA_SCOPE(LocOutputData, OutputData);
      constexpr I4 InputRank = ArrayT::rank;

      if (IsMeshDimension && InputRank >= 2) {
         OMEGA_SCOPE(LocMinLayer, MinLayer);
         OMEGA_SCOPE(LocMaxLayer, MaxLayer);
         OMEGA_SCOPE(LocInputData, InputData);
         const I4 LocLayerBoundOffset = LayerBoundOffset;

         if constexpr (InputRank == 2) {
            parallelForOuter(
                "TimeMeanAdd2D", LaunchConfig({NHorizOwned}),
                KOKKOS_LAMBDA(int IHoriz, const TeamMember &Team) {
                   const I4 KMin   = LocMinLayer(IHoriz);
                   const I4 KMax   = LocMaxLayer(IHoriz) + LocLayerBoundOffset;
                   const I4 KRange = vertRange(KMin, KMax);
                   parallelForInner(
                       Team, KRange, INNER_LAMBDA(int KIdx) {
                          const I4 K = KMin + KIdx;
                          LocOutputData(IHoriz, K) +=
                              static_cast<Real>(LocInputData(IHoriz, K));
                       });
                });
         } else if constexpr (InputRank == 3) {
            const I4 Dim0 = InputData.extent(0);
            parallelForOuter(
                "TimeMeanAdd3D", LaunchConfig({Dim0, NHorizOwned}),
                KOKKOS_LAMBDA(int I0, int IHoriz, const TeamMember &Team) {
                   const I4 KMin   = LocMinLayer(IHoriz);
                   const I4 KMax   = LocMaxLayer(IHoriz) + LocLayerBoundOffset;
                   const I4 KRange = vertRange(KMin, KMax);
                   parallelForInner(
                       Team, KRange, INNER_LAMBDA(int KIdx) {
                          const I4 K = KMin + KIdx;
                          LocOutputData(I0, IHoriz, K) +=
                              static_cast<Real>(LocInputData(I0, IHoriz, K));
                       });
                });
         }
      } else {
         // Flat fallback for scalars, rank-1, and non-mesh outputs
         parallelFor(
             {ArraySize}, KOKKOS_LAMBDA(const int FlatIdx) {
                LocOutputData.data()[FlatIdx] +=
                    static_cast<Real>(InputData.data()[FlatIdx]);
             });
      }
   } // end accumAdd

   /// Divides the accumulated sum by the sample count to form the mean. Uses
   /// the same active layer bounds and dispatch logic as accumInit/accumAdd so
   /// fill values in inactive layers are preserved through finalization.
   void finalize() {
      OMEGA_SCOPE(LocOutputData, OutputData);
      const Real InvNumAccum = 1.0 / static_cast<Real>(NumAccum);
      constexpr I4 InputRank = ArrayT::rank;

      if (IsMeshDimension && InputRank >= 2) {
         OMEGA_SCOPE(LocMinLayer, MinLayer);
         OMEGA_SCOPE(LocMaxLayer, MaxLayer);
         const I4 LocLayerBoundOffset = LayerBoundOffset;

         if constexpr (InputRank == 2) {
            parallelForOuter(
                "TimeMeanFinal2D", LaunchConfig({NHorizOwned}),
                KOKKOS_LAMBDA(int IHoriz, const TeamMember &Team) {
                   const I4 KMin   = LocMinLayer(IHoriz);
                   const I4 KMax   = LocMaxLayer(IHoriz) + LocLayerBoundOffset;
                   const I4 KRange = vertRange(KMin, KMax);
                   parallelForInner(
                       Team, KRange, INNER_LAMBDA(int KIdx) {
                          const I4 K = KMin + KIdx;
                          LocOutputData(IHoriz, K) *= InvNumAccum;
                       });
                });
         } else if constexpr (InputRank == 3) {
            const I4 Dim0 = LocOutputData.extent(0);
            parallelForOuter(
                "TimeMeanFinal3D", LaunchConfig({Dim0, NHorizOwned}),
                KOKKOS_LAMBDA(int I0, int IHoriz, const TeamMember &Team) {
                   const I4 KMin   = LocMinLayer(IHoriz);
                   const I4 KMax   = LocMaxLayer(IHoriz) + LocLayerBoundOffset;
                   const I4 KRange = vertRange(KMin, KMax);
                   parallelForInner(
                       Team, KRange, INNER_LAMBDA(int KIdx) {
                          const I4 K = KMin + KIdx;
                          LocOutputData(I0, IHoriz, K) *= InvNumAccum;
                       });
                });
         }
      } else {
         // Flat fallback for scalars, rank-1, and non-mesh outputs
         parallelFor(
             {ArraySize}, KOKKOS_LAMBDA(const int FlatIdx) {
                LocOutputData.data()[FlatIdx] *= InvNumAccum;
             });
      }
   } // end finalize

 private:
   /// Detaches the sample alarm, if any, from the clock it is attached to
   void detachSampleAlarm() {
      if (SampleClock != nullptr && OwnedSampleAlarm) {
         SampleClock->detachAlarm(OwnedSampleAlarm.get());
      }
      SampleClock = nullptr;
   } // end detachSampleAlarm

   /// Divides the averaging period into NumSamples equal sampling intervals.
   /// TimeInterval's integer division truncates a calendar interval, so a
   /// 1Day period divided by 4 would give zero days. Days have a fixed length
   /// and are divided in seconds instead. Months vary in length, so a period
   /// in months or years is divided in whole months and must divide evenly.
   TimeInterval
   dividePeriod(const TimeInterval &PeriodInterval, ///< [in] averaging period
                I4 NumSamples ///< [in] number of samples per period
   ) const {
      std::vector<std::string> Parsed = parseFreqStr(PeriodStr);
      const I8 Length                 = std::stoll(Parsed[0]);
      const TimeUnits Units           = TimeUnitsFromString(Parsed[1]);

      if (Units == TimeUnits::Days) {
         const I8 PeriodSeconds = Length * Calendar::getSecondsPerDay();
         return TimeInterval(PeriodSeconds, TimeUnits::Seconds) / NumSamples;
      }

      if (Units == TimeUnits::Months || Units == TimeUnits::Years) {
         const I8 PeriodMonths =
             (Units == TimeUnits::Years) ? Length * MONTHS_PER_YEAR : Length;
         if (PeriodMonths % NumSamples != 0) {
            ABORT_ERROR("TimeMeanOp: SamplingInterval '{}' must divide the "
                        "{} averaging period into whole months, since months "
                        "vary in length. Use a fixed interval such as 1day "
                        "instead.",
                        SamplingIntervalStr, PeriodStr);
         }
         return TimeInterval(PeriodMonths / NumSamples, TimeUnits::Months);
      }

      return PeriodInterval / NumSamples;
   } // end dividePeriod

   /// Output data array matching input layout but with Real type, used to
   /// accumulate sum and store the final time-averaged mean
   OutputArrayT OutputData;

   /// Number of values accumulated in current averaging period
   I4 NumAccum;

   /// Total size of input/output arrays for the flat parallel path
   I4 ArraySize;

   /// Number of owned points in the horizontal dimension (for MPI-correct
   /// hierarchical loops over mesh-distributed fields)
   I4 NHorizOwned;

   /// Whether the horizontal dimension is a recognized mesh dimension
   /// (NCells/NEdges/NVertices). When false, the flat fallback path is used.
   bool IsMeshDimension;

   /// Additional upper bound offset for interface fields on NVertLayersP1.
   I4 LayerBoundOffset = 0;

   /// Min active layer index for each horizontal point (rank >= 2 mesh fields)
   Array1DI4 MinLayer;

   /// Max active layer index for each horizontal point (rank >= 2 mesh fields)
   Array1DI4 MaxLayer;

   /// Alarm that signals when averaging period is complete (e.g., end of
   /// day/month)
   Alarm *PeriodAlarm;

   /// Flag indicating whether the next compute() should start a new averaging
   /// period
   bool IsNewPeriod;

   /// One-element host array backing the persisted NumAccum counter Field.
   /// Shares no storage with the member NumAccum; compute() copies the
   /// (normalized) count into it after each accumulation step and
   /// applyRestartState() copies it back after a restart read.
   HostArray1DI4 NumAccumData;

   /// Operator-owned recurring sample alarm, non-null only when a
   /// SamplingInterval is configured. The base-class SampleAlarm pointer
   /// borrows this instance; keeping it here ties its lifetime to the
   /// operator and keeps it attached to the model clock.
   std::unique_ptr<Alarm> OwnedSampleAlarm;

   /// The clock OwnedSampleAlarm is attached to, so that the destructor can
   /// detach it; null when no sample alarm is attached
   Clock *SampleClock = nullptr;

   /// Raw, un-normalized SamplingInterval string from config (empty => sample
   /// every timestep). Normalized to the base-class SamplingInterval in
   /// setPeriodAlarm() once the period and timestep are known.
   std::string SamplingIntervalStr;

   /// The averaging period as given in config (e.g. "1Month"), kept so that
   /// a "<n>perPeriod" sampling interval can be built in the period's units
   std::string PeriodStr;

}; // end class TimeMeanOp

} // end namespace OMEGA

#endif
