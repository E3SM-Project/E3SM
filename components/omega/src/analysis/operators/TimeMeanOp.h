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

      // Register output Field with same dimensions as input but Real type
      auto OutputField =
          Field::create(OutputNames[0],
                        "Time average of " + InputNames[0], // Description
                        "",                                 // Units
                        "",                                 // Standard name
                        -std::numeric_limits<Real>::max(),  // Min valid value
                        std::numeric_limits<Real>::max(),   // Max valid value
                        -std::numeric_limits<Real>::max(),  // Fill value
                        NDims,                              // Rank
                        DimNames                            // Dimension names
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

   } // end constructor

   /// Sets the period alarm that signals when to finalize the time average.
   /// Called by Analysis during initialization to provide the alarm associated
   /// with the averaging period specified in the constructor config.
   void setPeriodAlarm(Alarm *Alarm ///< [in] alarm for averaging period
                       ) override {
      PeriodAlarm = Alarm;
   }

   /// Computes the time-averaged mean by accumulating input values. On the
   /// first call after alarm rings (IsNewPeriod), copies input to output and
   /// sets count to 1. On subsequent calls, adds input to accumulated sum and
   /// increments count. When period alarm rings, divides accumulated sum by
   /// count to finalize mean and resets state for next period. Updates
   /// timestamp and computed flag.
   void compute(const TimeInstant &TimeStamp ///< [in] current timestamp
                ) override {
      // Create local scope reference to output array for kernel capture
      OMEGA_SCOPE(LocOutputData, OutputData);

      // Retrieve input Field and extract data array
      auto InputField = Field::get(InputNames[0]);
      auto InputData  = InputField->template getDataArray<ArrayT>();

      // Accumulate input values into output array
      if (IsNewPeriod) {
         // Start of new averaging period: initialize with first input
         // Cast to Real for accumulation
         NumAccum = 1;
         parallelFor(
             {ArraySize}, KOKKOS_LAMBDA(const int FlatIdx) {
                LocOutputData.data()[FlatIdx] =
                    static_cast<Real>(InputData.data()[FlatIdx]);
             });
         IsNewPeriod = false;

         // If the period alarm rings on the first sample of a period, finalize
         // immediately (mean == current value) and start a new period next call.
         if (PeriodAlarm != nullptr && PeriodAlarm->isRinging()) {
            IsNewPeriod = true;
         }
      } else {
         // Continue accumulation: add input to running sum
         parallelFor(
             {ArraySize}, KOKKOS_LAMBDA(const int FlatIdx) {
                LocOutputData.data()[FlatIdx] +=
                    static_cast<Real>(InputData.data()[FlatIdx]);
             });
         ++NumAccum;

         // Check if period alarm is ringing (time to finalize average)
         bool ShouldFinalize =
             (PeriodAlarm != nullptr && PeriodAlarm->isRinging());

         if (ShouldFinalize) {
            // Finalize: divide accumulated sum by count to get mean
            Real InvNumAccum = 1.0 / static_cast<Real>(NumAccum);
            parallelFor(
                {ArraySize}, KOKKOS_LAMBDA(const int FlatIdx) {
                   LocOutputData.data()[FlatIdx] *= InvNumAccum;
                });

            // Reset state for next averaging period
            IsNewPeriod = true;
         }
      }

      // Update cache validity markers
      LastComputed  = TimeStamp;
      FieldComputed = true;

   } // end compute

 private:
   /// Output data array matching input layout but with Real type, used to
   /// accumulate sum and store the final time-averaged mean
   OutputArrayT OutputData;

   /// Number of values accumulated in current averaging period
   I4 NumAccum;

   /// Total size of input/output arrays for parallel iteration
   I4 ArraySize;

   /// Alarm that signals when averaging period is complete (e.g., end of
   /// day/month)
   Alarm *PeriodAlarm;

   /// Flag indicating whether the next compute() should start a new averaging
   /// period
   bool IsNewPeriod;

}; // end class TimeMeanOp

} // end namespace OMEGA

#endif
