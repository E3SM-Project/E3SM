//===-- ocn/TracerDefs.inc --------------------*- C++ -*-===//
//
/// \file
/// \brief OMEGA Tracer Definition File
///
/// This file defines all OMEGA tracers. It is imported using
/// the C++ #include preprocessor directive in Tracers.cpp. Only the tracer
/// definitions selected in the Tracers section of the OMEGA YAML configuration
/// file will be allocated into memory.
//===----------------------------------------------------------------------===//

// Index defined for each tracer. The value of these tracer indices are still
// set in Tracers::init()  and tracers that have not been selected would have
// standard invalid value: Tracers::IndxInvalid = -1
inline static I4 IndxTemp   = Tracers::IndxInvalid;
inline static I4 IndxSalt   = Tracers::IndxInvalid;
inline static I4 IndxDebug1 = Tracers::IndxInvalid;
inline static I4 IndxDebug2 = Tracers::IndxInvalid;
inline static I4 IndxDebug3 = Tracers::IndxInvalid;

// Tracer definitions packaged in a defineAllTracers function
static void defineAllTracers() {

   define("Temperature",                        ///< [in] Name of tracer
          "Conservative Temperature",           ///< [in] Long name, description
          "degree_C",                           ///< [in] Units
          "sea_water_conservative_temperature", ///< [in] CF standard Name
          -273.15,                              ///< [in] min valid field value
          100.0,                                ///< [in] max valid field value
          1.e33,                                ///< [in] value of undef entries
          IndxTemp);                            ///< [out] (optional) static index

   define("Salinity", "Absolute Salinity", "g kg-1",
          "sea_water_absolute_salinity", 0.0, 50.0, 1.e33, IndxSalt);

   define("Debug1", "Debug Tracer 1", "none", "none", 0.0, 100.0, 1.e33,
          IndxDebug1);
   define("Debug2", "Debug Tracer 2", "none", "none", 0.0, 100.0, 1.e33,
          IndxDebug2);
   define("Debug3", "Debug Tracer 3", "none", "none", 0.0, 100.0, 1.e33,
          IndxDebug3);
}
