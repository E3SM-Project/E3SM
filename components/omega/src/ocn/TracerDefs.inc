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

define("Temp",                            ///< Name of tracer
       "Potential Temperature",           ///< Long name or description
       "degree_C",                        ///< Units
       "sea_water_potential_temperature", ///< CF standard Name
       -273.15,                           ///< min valid field value
       100.0,                             ///< max valid field value
       1.e33                              ///< value for undef entries
);

define("Salt", "Salinity", "psu", "sea_water_salinity", 0.0, 50.0, 1.e33);

define("Debug1", "Debug Tracer 1", "none", "none", 0.0, 100.0, 1.e33);
define("Debug2", "Debug Tracer 2", "none", "none", 0.0, 100.0, 1.e33);
define("Debug3", "Debug Tracer 3", "none", "none", 0.0, 100.0, 1.e33);
