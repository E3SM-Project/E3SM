/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_ENUMS_HPP
#define HOMMEXX_ENUMS_HPP

#include "Kokkos_Core.hpp"

#include <cstdint>

namespace Homme
{

// Convert strong typed enum to the underlying int value
// TODO: perhaps move this to Utility.hpp
template<typename E>
constexpr
KOKKOS_FORCEINLINE_FUNCTION
typename std::underlying_type<E>::type etoi(E e) {
  return static_cast<typename std::underlying_type<E>::type>(e);
}

// ============== Run options check utility enum ================== //

namespace Errors {

enum class ComparisonOp {
  EQ = 0,   // EQUAL
  NE,       // NOT EQUAL
  GT,       // GREATHER THAN
  LT,       // LESS THAN
  GE,       // GREATHER THAN OR EQUAL
  LE        // LESS THAN OR EQUAL
};

} // namespace Errors

// =================== Run parameters enums ====================== //

enum class ForcingAlg : int {
  FORCING_OFF =-1,
  FORCING_0   = 0, 
  FORCING_1   = 1, // Unsupported
  FORCING_2   = 2, // TODO: Rename FORCING_1 and FORCING_2 to something more descriptive
};

enum class MoistDry {
  MOIST,
  DRY
};

enum class AdvectionForm {
  Conservative,
  NonConservative
};

enum class RemapAlg {
  PPM_MIRRORED = 1,
  PPM_LIMITED_EXTRAP = 10
};

inline std::string remapAlg2str (const RemapAlg alg) {
  switch (alg) {
    case RemapAlg::PPM_MIRRORED:
      return "PPM Mirrored";
    case RemapAlg::PPM_LIMITED_EXTRAP:
      return "PPM Limited Extrapolation";
  }

  return "UNKNOWN";
}

// Currently unused. Test case forcings are done in the F90 code.
enum class TestCase {
  UNUSED, // be clear that we don't use this information
  ASP_BAROCLINIC,
  ASP_GRAVITY_WAVE,
  ASP_MOUNTAIN,
  ASP_ROSSBY,
  ASP_TRACER,
  BAROCLINIC,
  DCMIP2012_TEST1_1,
  DCMIP2012_TEST1_2,
  DCMIP2012_TEST1_3,
  DCMIP2012_TEST2_0,
  DCMIP2012_TEST2_1,
  DCMIP2012_TEST2_2,
  DCMIP2012_TEST3,
  HELD_SUAREZ0,
  JW_BAROCLINIC,
  DCMIP2016_TEST1,
  DCMIP2016_TEST2,
  DCMIP2016_TEST3
};

enum class UpdateType {
  LEAPFROG,
  FORWARD
};

enum class TimeStepType {
  // Explicit
  ttype5          =  5,

  // Implicit-Explicit
  ttype7_imex     =  7,
  ttype9_imex     =  9,
  ttype10_imex    = 10

};

// ======= How to combine output/input during calculations ========== //

enum class CombineMode {
  Replace  = 0,   // out = in
  Scale,          // out = alpha*in
  Update,         // out = beta*out + in
  ScaleUpdate,    // out = beta*out + alpha*in
  ScaleAdd,       // out = out + alpha*in (special case of ScaleUpdate with beta=1)
  Add,            // out = out + in (special case of ScaleAdd/Update with alpha=1, beta=1.0)
  Multiply,       // out = out*in
  Divide          // out = out/in
};

template<CombineMode CM>
inline std::string cm2str () {
  switch (CM) {
    case CombineMode::Replace:      return "Replace";
    case CombineMode::Scale:        return "Scale";
    case CombineMode::Update:       return "Update";
    case CombineMode::ScaleUpdate:  return "ScaleUpdate";
    case CombineMode::ScaleAdd:     return "ScaleAdd";
    case CombineMode::Add:          return "Add";
    case CombineMode::Multiply:     return "Multiply";
    case CombineMode::Divide:       return "Divide";
  }
  return "UNKNOWN";
}

// =================== Euler Step DSS Option ====================== //

enum class DSSOption {
  ETA,
  OMEGA,
  DIV_VDP_AVE
};

// =================== Mesh connectivity enums ====================== //

// The kind of connection: edge, corner or missing (one of the corner connections on one of the 8 cube vertices)
enum class ConnectionKind : std::uint8_t {
  EDGE    = 0,
  CORNER  = 1,
  MISSING = 2,  // Used to detect missing connections
  ANY     = 3   // Used when the kind of connection is not needed
};

// The locality of connection: local, shared or missing
enum class ConnectionSharing : std::uint8_t {
  LOCAL   = 0,
  SHARED  = 1,
  MISSING = 2,  // Used to detect missing connections
  ANY     = 3   // Used when the kind of connection is not needed
};

enum class ConnectionName : std::uint8_t {
  // Edges
  SOUTH = 0,
  NORTH = 1,
  WEST  = 2,
  EAST  = 3,

  // Corners
  SWEST = 4,
  SEAST = 5,
  NWEST = 6,
  NEAST = 7
};

// Direction (useful only for an edge)
constexpr int NUM_DIRECTIONS = 3;
enum class Direction : std::uint8_t {
  FORWARD  = 0,
  BACKWARD = 1,
  INVALID  = 2
};


} // namespace Homme

#endif // HOMMEXX_ENUMS_HPP
