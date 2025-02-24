#ifndef MAM_AEROSOL_OPTICS_READ_TABLES_HPP
#define MAM_AEROSOL_OPTICS_READ_TABLES_HPP

#include "ekat/ekat_parameter_list.hpp"
#include "mam_coupling.hpp"
#include "share/field/field_manager.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/grid/grids_manager.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
#include "share/io/scorpio_input.hpp"

// later to mam_coupling.hpp
namespace scream::mam_coupling {

using view_1d_host    = typename KT::view_1d<Real>::HostMirror;
using view_2d_host    = typename KT::view_2d<Real>::HostMirror;
using view_5d_host    = typename KT::view_ND<Real, 5>::HostMirror;
using complex_view_1d = typename KT::view_1d<Kokkos::complex<Real>>;

constexpr int nlwbands = mam4::modal_aer_opt::nlwbands;
constexpr int nswbands = mam4::modal_aer_opt::nswbands;

struct AerosolOpticsHostData {
  // host views
  view_2d_host refindex_real_sw_host;
  view_2d_host refindex_im_sw_host;
  view_2d_host refindex_real_lw_host;
  view_2d_host refindex_im_lw_host;

  view_5d_host absplw_host;
  view_5d_host abspsw_host;
  view_5d_host asmpsw_host;
  view_5d_host extpsw_host;
};

using AerosolOpticsDeviceData = mam4::modal_aer_opt::AerosolOpticsDeviceData;

inline void set_parameters_table(
    AerosolOpticsHostData &aerosol_optics_host_data,
    ekat::ParameterList &rrtmg_params,
    std::map<std::string, FieldLayout> &layouts,
    std::map<std::string, view_1d_host> &host_views) {
  // Set up input structure to read data from file.
  using strvec_t = std::vector<std::string>;
  using namespace ShortFieldTagsNames;

  constexpr int refindex_real = mam4::modal_aer_opt::refindex_real;
  constexpr int refindex_im   = mam4::modal_aer_opt::refindex_im;
  constexpr int coef_number   = mam4::modal_aer_opt::coef_number;

  auto make_layout = [](const std::vector<int> &extents,
                        const std::vector<std::string> &names) {
    std::vector<FieldTag> tags(extents.size(), CMP);
    return FieldLayout(tags, extents, names);
  };

  auto refindex_real_lw_host =
      view_2d_host("refrtablw_real_host", nlwbands, refindex_real);
  auto refindex_im_lw_host =
      view_2d_host("refrtablw_im_host", nlwbands, refindex_im);

  auto refindex_real_sw_host =
      view_2d_host("refrtabsw_real_host", nswbands, refindex_real);
  auto refindex_im_sw_host =
      view_2d_host("refrtabsw_im_host", nswbands, refindex_im);

  // absplw(lw_band, mode, refindex_im, refindex_real, coef_number)
  auto absplw_host = view_5d_host("absplw_host", nlwbands, 1, refindex_im,
                                  refindex_real, coef_number);

  auto asmpsw_host = view_5d_host("asmpsw_host", nswbands, 1, refindex_im,
                                  refindex_real, coef_number);
  auto extpsw_host = view_5d_host("extpsw_host", nswbands, 1, refindex_im,
                                  refindex_real, coef_number);
  auto abspsw_host = view_5d_host("abspsw_host", nswbands, 1, refindex_im,
                                  refindex_real, coef_number);

  aerosol_optics_host_data.refindex_real_lw_host = refindex_real_lw_host;
  aerosol_optics_host_data.refindex_im_lw_host   = refindex_im_lw_host;
  aerosol_optics_host_data.refindex_real_sw_host = refindex_real_sw_host;
  aerosol_optics_host_data.refindex_im_sw_host   = refindex_im_sw_host;
  aerosol_optics_host_data.absplw_host           = absplw_host;
  aerosol_optics_host_data.asmpsw_host           = asmpsw_host;
  aerosol_optics_host_data.extpsw_host           = extpsw_host;
  aerosol_optics_host_data.abspsw_host           = abspsw_host;

  auto refindex_real_lw_layout =
      make_layout({nlwbands, refindex_real}, {"lwband", "refindex_real"});
  auto refindex_im_lw_layout =
      make_layout({nlwbands, refindex_im}, {"lwband", "refindex_im"});
  auto refindex_real_sw_layout =
      make_layout({nswbands, refindex_real}, {"swband", "refindex_real"});
  auto refindex_im_sw_layout =
      make_layout({nswbands, refindex_im}, {"swband", "refindex_im"});
  auto absplw_layout = make_layout(
      {nlwbands, 1, refindex_im, refindex_real, coef_number},
      {"lwband", "mode", "refindex_im", "refindex_real", "coef_number"});
  // use also for extpsw, abspsw
  auto asmpsw_layout = make_layout(
      {nswbands, 1, refindex_im, refindex_real, coef_number},
      {"swband", "mode", "refindex_im", "refindex_real", "coef_number"});

  rrtmg_params.set<strvec_t>(
      "Field Names",
      {"asmpsw", "extpsw", "abspsw", "absplw", "refindex_real_sw",
       "refindex_im_sw", "refindex_real_lw", "refindex_im_lw"});

  rrtmg_params.set("Skip_Grid_Checks", true);

  host_views["refindex_real_sw"] =
      view_1d_host(refindex_real_sw_host.data(), refindex_real_sw_host.size());

  host_views["refindex_im_sw"] =
      view_1d_host(refindex_im_sw_host.data(), refindex_im_sw_host.size());

  host_views["refindex_real_lw"] =
      view_1d_host(refindex_real_lw_host.data(), refindex_real_lw_host.size());

  host_views["refindex_im_lw"] =
      view_1d_host(refindex_im_lw_host.data(), refindex_im_lw_host.size());

  host_views["absplw"] = view_1d_host(absplw_host.data(), absplw_host.size());

  host_views["asmpsw"] = view_1d_host(asmpsw_host.data(), asmpsw_host.size());

  host_views["extpsw"] = view_1d_host(extpsw_host.data(), extpsw_host.size());

  host_views["abspsw"] = view_1d_host(abspsw_host.data(), abspsw_host.size());

  layouts.emplace("refindex_real_lw", refindex_real_lw_layout);
  layouts.emplace("refindex_im_lw", refindex_im_lw_layout);
  layouts.emplace("refindex_real_sw", refindex_real_sw_layout);
  layouts.emplace("refindex_im_sw", refindex_im_sw_layout);
  layouts.emplace("absplw", absplw_layout);
  layouts.emplace("asmpsw", asmpsw_layout);
  layouts.emplace("extpsw", asmpsw_layout);
  layouts.emplace("abspsw", asmpsw_layout);
}
// KOKKOS_INLINE_FUNCTION
inline void read_rrtmg_table(
    const std::string &table_filename, const int imode,
    ekat::ParameterList &params,
    const std::shared_ptr<const AbstractGrid> &grid,
    const std::map<std::string, view_1d_host> &host_views_1d,
    const std::map<std::string, FieldLayout> &layouts,
    const AerosolOpticsHostData &aerosol_optics_host_data,
    const AerosolOpticsDeviceData &aerosol_optics_device_data) {
  constexpr int refindex_real = mam4::modal_aer_opt::refindex_real;
  constexpr int refindex_im   = mam4::modal_aer_opt::refindex_im;
  constexpr int coef_number   = mam4::modal_aer_opt::coef_number;

  using view_3d_host = typename KT::view_3d<Real>::HostMirror;

  // temp views:
  view_3d_host temp_lw_3d_host("temp_absplw_host", coef_number, refindex_real,
                               refindex_im);

  params.set("Filename", table_filename);
  AtmosphereInput rrtmg(params, grid, host_views_1d, layouts);
  rrtmg.read_variables();
  rrtmg.finalize();

  // copy data from host to device for mode 1
  int d1 = imode;
  for(int d3 = 0; d3 < nswbands; ++d3) {
    auto real_host_d3 =
        ekat::subview(aerosol_optics_host_data.refindex_real_sw_host, d3);
    Kokkos::deep_copy(aerosol_optics_device_data.refrtabsw[d1][d3],
                      real_host_d3);
    auto im_host_d3 =
        ekat::subview(aerosol_optics_host_data.refindex_im_sw_host, d3);
    Kokkos::deep_copy(aerosol_optics_device_data.refitabsw[d1][d3], im_host_d3);
  }  // d3

  for(int d3 = 0; d3 < nlwbands; ++d3) {
    auto real_host_d3 =
        ekat::subview(aerosol_optics_host_data.refindex_real_lw_host, d3);
    Kokkos::deep_copy(aerosol_optics_device_data.refrtablw[d1][d3],
                      real_host_d3);
    auto im_host_d3 =
        ekat::subview(aerosol_optics_host_data.refindex_im_lw_host, d3);
    Kokkos::deep_copy(aerosol_optics_device_data.refitablw[d1][d3], im_host_d3);
  }  // d3

  // NOTE: we need to reorder dimensions in absplw
  // netcfd : (lw_band, mode, refindex_im, refindex_real, coef_number)
  // mam4xx : (mode, lw_band, coef_number, refindex_real, refindex_im )
  // e3sm : (ntot_amode,coef_number,refindex_real,refindex_im,nlwbands)

  for(int d5 = 0; d5 < nlwbands; ++d5) {
    // reshape data:
    Kokkos::parallel_for(
        "reshaping absplw",
        Kokkos::MDRangePolicy<Kokkos::Rank<3>,
                              Kokkos::DefaultHostExecutionSpace>(
            {0, 0, 0}, {coef_number, refindex_real, refindex_im}),
        [&](const int d2, const int d3, const int d4) {
          temp_lw_3d_host(d2, d3, d4) =
              aerosol_optics_host_data.absplw_host(d5, 0, d4, d3, d2);
        });
    Kokkos::fence();

    // syn data to device
    Kokkos::deep_copy(aerosol_optics_device_data.absplw[d1][d5],
                      temp_lw_3d_host);
  }  // d5
  // asmpsw, abspsw, extpsw
  // netcfd : (sw_band, mode, refindex_im, refindex_real, coef_number)
  // mam4xx : (mode, sw_band, coef_number, refindex_real, refindex_im )

  for(int d5 = 0; d5 < nswbands; ++d5) {
    // reshape data
    Kokkos::parallel_for(
        "reshaping asmpsw",
        Kokkos::MDRangePolicy<Kokkos::Rank<3>,
                              Kokkos::DefaultHostExecutionSpace>(
            {0, 0, 0}, {coef_number, refindex_real, refindex_im}),
        [&](const int d2, const int d3, const int d4) {
          temp_lw_3d_host(d2, d3, d4) =
              aerosol_optics_host_data.asmpsw_host(d5, 0, d4, d3, d2);
        });
    Kokkos::fence();
    // syn data to device
    Kokkos::deep_copy(aerosol_optics_device_data.asmpsw[d1][d5],
                      temp_lw_3d_host);
    // reshape data
    Kokkos::parallel_for(
        "reshaping abspsw",
        Kokkos::MDRangePolicy<Kokkos::Rank<3>,
                              Kokkos::DefaultHostExecutionSpace>(
            {0, 0, 0}, {coef_number, refindex_real, refindex_im}),
        [&](const int d2, const int d3, const int d4) {
          temp_lw_3d_host(d2, d3, d4) =
              aerosol_optics_host_data.abspsw_host(d5, 0, d4, d3, d2);
        });
    Kokkos::fence();
    // syn data to device
    Kokkos::deep_copy(aerosol_optics_device_data.abspsw[d1][d5],
                      temp_lw_3d_host);
    // reshape data
    Kokkos::parallel_for(
        "reshaping extpsw",
        Kokkos::MDRangePolicy<Kokkos::Rank<3>,
                              Kokkos::DefaultHostExecutionSpace>(
            {0, 0, 0}, {coef_number, refindex_real, refindex_im}),
        [&](const int d2, const int d3, const int d4) {
          temp_lw_3d_host(d2, d3, d4) =
              aerosol_optics_host_data.extpsw_host(d5, 0, d4, d3, d2);
        });

    Kokkos::fence();
    // syn data to device
    Kokkos::deep_copy(aerosol_optics_device_data.extpsw[d1][d5],
                      temp_lw_3d_host);

  }  // d5
}

inline void read_water_refindex(const std::string &table_filename,
                                const std::shared_ptr<const AbstractGrid> &grid,
                                const complex_view_1d &crefwlw,
                                const complex_view_1d &crefwsw) {
  // refractive index for water read in read_water_refindex
  // crefwsw(nswbands) ! complex refractive index for water visible
  // crefwlw(nlwbands) ! complex refractive index for water infrared

  using namespace ShortFieldTagsNames;
  using view_1d_host = typename KT::view_1d<Real>::HostMirror;
  // Set up input structure to read data from file.
  using strvec_t = std::vector<std::string>;

  // here a made a list of variables that I want to read from netcdf files
  ekat::ParameterList params;
  params.set("Filename", table_filename);
  params.set("Skip_Grid_Checks", true);

  params.set<strvec_t>("Field Names",
                       {"refindex_im_water_lw", "refindex_im_water_sw",
                        "refindex_real_water_lw", "refindex_real_water_sw"});
  // make a list of host views
  std::map<std::string, view_1d_host> host_views_water;
  // fist allocate host views.
  view_1d_host refindex_im_water_sw_host("refindex_im_water_sw_host", nswbands);
  view_1d_host refindex_real_water_sw_host("refindex_real_water_sw_host",
                                           nswbands);
  view_1d_host refindex_im_water_lw_host("refindex_im_water_lw_host", nlwbands);
  view_1d_host refindex_real_water_lw_host("refindex_real_water_lw_host",
                                           nlwbands);

  host_views_water["refindex_im_water_sw"]   = refindex_im_water_sw_host;
  host_views_water["refindex_real_water_sw"] = refindex_real_water_sw_host;
  host_views_water["refindex_im_water_lw"]   = refindex_im_water_lw_host;
  host_views_water["refindex_real_water_lw"] = refindex_real_water_lw_host;

  // defines layouts
  std::map<std::string, FieldLayout> layouts_water;
  FieldLayout refindex_water_sw_layout{{CMP}, {nswbands}, {"swband"}};
  FieldLayout refindex_water_lw_layout{{CMP}, {nlwbands}, {"lwband"}};

  layouts_water.emplace("refindex_im_water_sw", refindex_water_sw_layout);
  layouts_water.emplace("refindex_real_water_sw", refindex_water_sw_layout);
  layouts_water.emplace("refindex_im_water_lw", refindex_water_lw_layout);
  layouts_water.emplace("refindex_real_water_lw", refindex_water_lw_layout);

  // create a object to read data
  AtmosphereInput refindex_water(params, grid, host_views_water, layouts_water);
  refindex_water.read_variables();
  refindex_water.finalize();

  //  maybe make a 1D vied of Kokkos::complex<Real>
  const auto crefwlw_host = Kokkos::create_mirror_view(crefwlw);
  const auto crefwsw_host = Kokkos::create_mirror_view(crefwsw);
  for(int i = 0; i < nlwbands; ++i) {
    // Kokkos::complex<Real> temp;
    crefwlw_host(i).real() = refindex_real_water_lw_host(i);
    crefwlw_host(i).imag() = haero::abs(refindex_im_water_lw_host(i));
  }
  Kokkos::deep_copy(crefwlw, crefwlw_host);
  // set complex representation of refractive indices as module data
  for(int i = 0; i < nswbands; ++i) {
    // Kokkos::complex<Real> temp;
    crefwsw_host(i).real() = refindex_real_water_sw_host(i);
    crefwsw_host(i).imag() = haero::abs(refindex_im_water_sw_host(i));
  }
  Kokkos::deep_copy(crefwsw, crefwsw_host);
}
// read_refindex_aero

inline void set_refindex_names(std::string surname, ekat::ParameterList &params,
                               std::map<std::string, view_1d_host> &host_views,
                               std::map<std::string, FieldLayout> &layouts) {
  // set variables names
  using view_1d_host = typename KT::view_1d<Real>::HostMirror;
  using strvec_t     = std::vector<std::string>;
  using namespace ShortFieldTagsNames;

  std::string refindex_real_sw = "refindex_real_" + surname + "_sw";
  std::string refindex_im_sw   = "refindex_im_" + surname + "_sw";
  std::string refindex_real_lw = "refindex_real_" + surname + "_lw";
  std::string refindex_im_lw   = "refindex_im_" + surname + "_lw";

  params.set("Skip_Grid_Checks", true);
  params.set<strvec_t>("Field Names", {refindex_real_sw, refindex_im_sw,
                                       refindex_real_lw, refindex_im_lw});
  // allocate host views
  host_views[refindex_real_sw] = view_1d_host(refindex_real_sw, nswbands);
  host_views[refindex_im_sw]   = view_1d_host(refindex_im_sw, nswbands);
  host_views[refindex_real_lw] = view_1d_host(refindex_real_lw, nlwbands);
  host_views[refindex_im_lw]   = view_1d_host(refindex_im_lw, nlwbands);

  FieldLayout scalar_refindex_sw_layout{{CMP}, {nswbands}, {"swband"}};
  FieldLayout scalar_refindex_lw_layout{{CMP}, {nlwbands}, {"lwband"}};

  layouts.emplace(refindex_real_sw, scalar_refindex_sw_layout);
  layouts.emplace(refindex_im_sw, scalar_refindex_sw_layout);
  layouts.emplace(refindex_real_lw, scalar_refindex_lw_layout);
  layouts.emplace(refindex_im_lw, scalar_refindex_lw_layout);

}  // set_refindex_aero

inline void set_refindex_aerosol(
    const int species_id, std::map<std::string, view_1d_host> &host_views,
    mam_coupling::complex_view_2d::HostMirror
        &specrefndxsw_host,  // complex refractive index for water visible
    mam_coupling::complex_view_2d::HostMirror &specrefndxlw_host) {
  std::string sw_real_name = "refindex_real_aer_sw";
  std::string lw_real_name = "refindex_real_aer_lw";
  std::string sw_im_name   = "refindex_im_aer_sw";
  std::string lw_im_name   = "refindex_im_aer_lw";

  for(int i = 0; i < nswbands; i++) {
    specrefndxsw_host(i, species_id).real() = host_views[sw_real_name](i);
    specrefndxsw_host(i, species_id).imag() =
        haero::abs(host_views[sw_im_name](i));
  }
  for(int i = 0; i < nlwbands; i++) {
    specrefndxlw_host(i, species_id).real() = host_views[lw_real_name](i);
    specrefndxlw_host(i, species_id).imag() =
        haero::abs(host_views[lw_im_name](i));
  }

}  // copy_refindex_to_device

}  // namespace scream::mam_coupling

#endif
