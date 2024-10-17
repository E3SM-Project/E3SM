#ifndef FRACTIONAL_LANDUSE_HPP
#define FRACTIONAL_LANDUSE_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/grid/remap/abstract_remapper.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/iop/intensive_observation_period.hpp"
#include "share/scream_types.hpp"
#include "share/util/scream_time_stamp.hpp"

namespace scream {
namespace frac_landuse {

template <typename ScalarType, typename DeviceType>
struct fracLandUseFunctions {
  using Device = DeviceType;

  using KT         = KokkosTypes<Device>;
  using MemberType = typename KT::MemberType;
  using view_2d    = typename KT::template view_2d<Real>;

  // -------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------

  struct FracLandUseData {
    FracLandUseData() = default;
    FracLandUseData(const int ncol_, const int nclass_) {
      init(ncol_, nclass_, true);
    }
    void init(const int ncol_, const int nclass_, const bool allocate) {
      ncols  = ncol_;
      nclass = nclass_;
      if(allocate) FRAC_LAND_USE = view_2d("FRAC_LAND_USE", nclass, ncols);
    }
    // Basic spatial dimensions of the data
    int ncols;
    int nclass;

    view_2d FRAC_LAND_USE;  // Fractional land use (unitless) dimensions =
                            // (nclass, ncols)
  };                        // FracLandUseData

  // -------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------

  struct FracLandUseInput {
    FracLandUseInput() = default;
    FracLandUseInput(const int ncols_, const int nclass_) {
      init(ncols_, nclass_);
    }
    void init(const int ncols_, const int nclass_) {
      data.init(ncols_, nclass_, true);
    }
    FracLandUseData data;  // All fractional land use fields
  };                       // FracLandUseInput

  // -------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------

  // The output is really just SPAData, but for clarity it might
  // help to see a SPAOutput along a SPAInput in functions signatures
  using FracLandUseOutput = FracLandUseData;

  // -------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------

  // Fractional land use routines
  static std::shared_ptr<AbstractRemapper> create_horiz_remapper(
      const std::shared_ptr<const AbstractGrid> &model_grid,
      const std::string &fracLandUse_data_file, const std::string &map_file);

  // -------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------

  static std::shared_ptr<AtmosphereInput> create_data_reader(
      const std::shared_ptr<AbstractRemapper> &horiz_remapper,
      const std::string &data_file);
#if 0
  // -------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------

  static void fracLandUse_main(const fracLandUseTimeState &time_state,
                               const fracLandUseInput &data_beg,
                               const fracLandUseInput &data_end,
                               const fracLandUseOutput &data_out);

  // -------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------

  static void update_fracLandUse_data_from_file(
      std::shared_ptr<AtmosphereInput> &scorpio_reader,
      const util::TimeStamp &ts,
      const int time_index,  // zero-based
      AbstractRemapper &fracLandUse_horiz_interp,
      fracLandUseInput &fracLandUse_input);

  // -------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------

  static void update_fracLandUse_timestate(
      std::shared_ptr<AtmosphereInput> &scorpio_reader,
      const util::TimeStamp &ts, AbstractRemapper &fracLandUse_horiz_interp,
      fracLandUseTimeState &time_state, fracLandUseInput &fracLandUse_beg,
      fracLandUseInput &fracLandUse_end);

  // -------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------

  // The following three are called during fracLandUse_main
  static void perform_time_interpolation(const fracLandUseTimeState &time_state,
                                         const fracLandUseInput &data_beg,
                                         const fracLandUseInput &data_end,
                                         const fracLandUseOutput &data_out);

  // Performs convex interpolation of x0 and x1 at point t
  template <typename ScalarX, typename ScalarT>
  KOKKOS_INLINE_FUNCTION static ScalarX linear_interp(const ScalarX &x0,
                                                      const ScalarX &x1,
                                                      const ScalarT &t);
  // -------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------
#endif

  static void init_frac_landuse_file_read(
      const int ncol, const int nclass,
      const std::shared_ptr<const AbstractGrid> &grid,
      const std::string &data_file, const std::vector<std::string> &sectors,
      const std::string &mapping_file,
      // output
      std::shared_ptr<AbstractRemapper> &FracLandUseHorizInterp,
      FracLandUseInput &FracLandUseData_start,
      FracLandUseInput &FracLandUseData_end,
      FracLandUseOutput &FracLandUseData_out,
      std::shared_ptr<AtmosphereInput> &FracLandUseDataReader);

};  // struct fracLandUseFunctions

}  // namespace frac_landuse
}  // namespace scream
#endif  // FRACTIONAL_LANDUSE_HPP

#include "fractional_land_use_impl.hpp"