#ifndef MAM_AEROSOL_OPTICS_READ_TABLES_HPP
#define MAM_AEROSOL_OPTICS_READ_TABLES_HPP

#include "mam_coupling.hpp"
#include "share/manager/field_manager.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/grid/grids_manager.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
#include "share/io/scorpio_input.hpp"

#include <ekat_parameter_list.hpp>

// later to mam_coupling.hpp
namespace scream::mam_coupling {

inline Field
make_field (const std::string& name,
            const FieldLayout& layout,
            const std::shared_ptr<const AbstractGrid>& grid)
{
  const auto units = ekat::units::Units::nondimensional();
  FieldIdentifier fid(name,layout,units,grid->name());
  Field f(fid);
  f.allocate_view();
  return f;
};

using view_2d_host    = typename KT::view_2d<Real>::HostMirror;
using view_5d_host    = typename KT::view_ND<Real, 5>::HostMirror;
using complex_view_1d = typename KT::view_1d<Kokkos::complex<Real>>;

constexpr int nlwbands = mam4::modal_aer_opt::nlwbands;
constexpr int nswbands = mam4::modal_aer_opt::nswbands;

using AerosolOpticsDeviceData = mam4::modal_aer_opt::AerosolOpticsDeviceData;

inline std::map<std::string,Field>
create_optics_fields(const std::shared_ptr<const AbstractGrid> &grid)
{
  // Set up input structure to read data from file.
  using namespace ShortFieldTagsNames;

  constexpr int refindex_real = mam4::modal_aer_opt::refindex_real;
  constexpr int refindex_im   = mam4::modal_aer_opt::refindex_im;
  constexpr int coef_number   = mam4::modal_aer_opt::coef_number;

  auto make_layout = [](const std::vector<int> &extents,
                        const std::vector<std::string> &names) {
    std::vector<FieldTag> tags(extents.size(), CMP);
    return FieldLayout(tags, extents, names);
  };

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

  std::map<std::string,Field> data;

  data["refindex_real_lw"] = make_field("refindex_real_lw",refindex_real_lw_layout,grid);
  data["refindex_im_lw"]   = make_field("refindex_im_lw",  refindex_im_lw_layout,grid);
  data["refindex_real_sw"] = make_field("refindex_real_sw",refindex_real_sw_layout,grid);
  data["refindex_im_sw"]   = make_field("refindex_im_sw",  refindex_im_sw_layout,grid);
  data["absplw"]           = make_field("absplw",          absplw_layout,grid);
  data["asmpsw"]           = make_field("asmpsw",          asmpsw_layout,grid);
  data["extpsw"]           = make_field("extpsw",          asmpsw_layout,grid);
  data["abspsw"]           = make_field("abspsw",          asmpsw_layout,grid);

  return data;
}

inline void read_rrtmg_table(
    const std::string &table_filename, const int imode,
    const std::shared_ptr<const AbstractGrid> &grid,
    const std::map<std::string,Field> &aerosol_optics_fields,
    const AerosolOpticsDeviceData &aerosol_optics_device_data) {
  constexpr int refindex_real = mam4::modal_aer_opt::refindex_real;
  constexpr int refindex_im   = mam4::modal_aer_opt::refindex_im;
  constexpr int coef_number   = mam4::modal_aer_opt::coef_number;

  using view_3d_host = typename KT::view_3d<Real>::HostMirror;

  // temp views:
  view_3d_host temp_lw_3d_host("temp_absplw_host", coef_number, refindex_real,
                               refindex_im);

  AtmosphereInput rrtmg(table_filename, grid, aerosol_optics_fields, true);
  rrtmg.read_variables();
  rrtmg.finalize();

  // copy data from host to device for mode 1
  // TODO: why can't we copy device data directly?
  int d1 = imode;
  auto refindex_real_sw_host = aerosol_optics_fields.at("refindex_real_sw").get_view<const Real**,Host>();
  auto refindex_im_sw_host   = aerosol_optics_fields.at("refindex_im_sw").get_view<const Real**,Host>();
  auto refindex_real_lw_host = aerosol_optics_fields.at("refindex_real_lw").get_view<const Real**,Host>();
  auto refindex_im_lw_host   = aerosol_optics_fields.at("refindex_im_lw").get_view<const Real**,Host>();
  for(int d3 = 0; d3 < nswbands; ++d3) {
    auto real_host_d3 = ekat::subview(refindex_real_sw_host, d3);
    Kokkos::deep_copy(aerosol_optics_device_data.refrtabsw[d1][d3],
                      real_host_d3);
    auto im_host_d3 = ekat::subview(refindex_im_sw_host, d3);
    Kokkos::deep_copy(aerosol_optics_device_data.refitabsw[d1][d3], im_host_d3);
  }  // d3

  for(int d3 = 0; d3 < nlwbands; ++d3) {
    auto real_host_d3 = ekat::subview(refindex_real_lw_host, d3);
    Kokkos::deep_copy(aerosol_optics_device_data.refrtablw[d1][d3],
                      real_host_d3);
    auto im_host_d3 = ekat::subview(refindex_im_lw_host, d3);
    Kokkos::deep_copy(aerosol_optics_device_data.refitablw[d1][d3], im_host_d3);
  }  // d3

  // NOTE: we need to reorder dimensions in absplw
  // netcfd : (lw_band, mode, refindex_im, refindex_real, coef_number)
  // mam4xx : (mode, lw_band, coef_number, refindex_real, refindex_im )
  // e3sm : (ntot_amode,coef_number,refindex_real,refindex_im,nlwbands)

  auto absplw_host = aerosol_optics_fields.at("absplw").get_view<const Real*****,Host>();
  for(int d5 = 0; d5 < nlwbands; ++d5) {
    // reshape data:
    Kokkos::parallel_for(
        "reshaping absplw",
        Kokkos::MDRangePolicy<Kokkos::Rank<3>,
                              Kokkos::DefaultHostExecutionSpace>(
            {0, 0, 0}, {coef_number, refindex_real, refindex_im}),
        [&](const int d2, const int d3, const int d4) {
          temp_lw_3d_host(d2, d3, d4) = absplw_host(d5, 0, d4, d3, d2);
        });
    Kokkos::fence();

    // syn data to device
    Kokkos::deep_copy(aerosol_optics_device_data.absplw[d1][d5],
                      temp_lw_3d_host);
  }  // d5
  // asmpsw, abspsw, extpsw
  // netcfd : (sw_band, mode, refindex_im, refindex_real, coef_number)
  // mam4xx : (mode, sw_band, coef_number, refindex_real, refindex_im )

  auto asmpsw_host = aerosol_optics_fields.at("asmpsw").get_view<const Real*****,Host>();
  auto abspsw_host = aerosol_optics_fields.at("abspsw").get_view<const Real*****,Host>();
  auto extpsw_host = aerosol_optics_fields.at("extpsw").get_view<const Real*****,Host>();
  for(int d5 = 0; d5 < nswbands; ++d5) {
    // reshape data
    Kokkos::parallel_for(
        "reshaping asmpsw",
        Kokkos::MDRangePolicy<Kokkos::Rank<3>,
                              Kokkos::DefaultHostExecutionSpace>(
            {0, 0, 0}, {coef_number, refindex_real, refindex_im}),
        [&](const int d2, const int d3, const int d4) {
          temp_lw_3d_host(d2, d3, d4) = asmpsw_host(d5, 0, d4, d3, d2);
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
          temp_lw_3d_host(d2, d3, d4) = abspsw_host(d5, 0, d4, d3, d2);
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
          temp_lw_3d_host(d2, d3, d4) = extpsw_host(d5, 0, d4, d3, d2);
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
  // Set up input structure to read data from file.

  // here a made a list of variables that I want to read from netcdf files
  FieldLayout refindex_water_sw_layout{{CMP}, {nswbands}, {"swband"}};
  FieldLayout refindex_water_lw_layout{{CMP}, {nlwbands}, {"lwband"}};

  auto refindex_im_water_sw   = make_field("refindex_im_water_sw", refindex_water_sw_layout,grid);
  auto refindex_real_water_sw = make_field("refindex_real_water_sw", refindex_water_sw_layout,grid);
  auto refindex_im_water_lw   = make_field("refindex_im_water_lw", refindex_water_lw_layout,grid);
  auto refindex_real_water_lw = make_field("refindex_real_water_lw", refindex_water_lw_layout,grid);

  std::vector<Field> fields = {
    refindex_im_water_sw,
    refindex_real_water_sw,
    refindex_im_water_lw,
    refindex_real_water_lw
  };

  // create a object to read data
  AtmosphereInput refindex_water(table_filename, grid, fields, true);
  refindex_water.read_variables();
  refindex_water.finalize();

  //  maybe make a 1D vied of Kokkos::complex<Real>
  const auto crefwlw_host = Kokkos::create_mirror_view(crefwlw);
  const auto crefwsw_host = Kokkos::create_mirror_view(crefwsw);
  const auto refindex_im_water_sw_host   = refindex_im_water_sw.get_view<const Real*,Host>();
  const auto refindex_real_water_sw_host = refindex_real_water_sw.get_view<const Real*,Host>();
  const auto refindex_im_water_lw_host   = refindex_im_water_lw.get_view<const Real*,Host>();
  const auto refindex_real_water_lw_host = refindex_real_water_lw.get_view<const Real*,Host>();
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

inline std::map<std::string,Field>
create_refindex_fields (const std::string& surname,
                        const std::shared_ptr<const AbstractGrid>& grid)
{
  using namespace ShortFieldTagsNames;

  std::string refindex_real_sw = "refindex_real_" + surname + "_sw";
  std::string refindex_im_sw   = "refindex_im_" + surname + "_sw";
  std::string refindex_real_lw = "refindex_real_" + surname + "_lw";
  std::string refindex_im_lw   = "refindex_im_" + surname + "_lw";

  FieldLayout refindex_sw_layout{{CMP}, {nswbands}, {"swband"}};
  FieldLayout refindex_lw_layout{{CMP}, {nlwbands}, {"lwband"}};

  std::map<std::string,Field> fields = {
    {refindex_real_sw, make_field(refindex_real_sw,refindex_sw_layout,grid)},
    {refindex_im_sw,   make_field(refindex_im_sw,  refindex_sw_layout,grid)},
    {refindex_real_lw, make_field(refindex_real_lw,refindex_lw_layout,grid)},
    {refindex_im_lw,   make_field(refindex_im_lw,  refindex_lw_layout,grid)}
  };
  return fields;
}

inline void set_refindex_aerosol(
    const int species_id,
    std::map<std::string, Field> fields,
    mam_coupling::complex_view_2d::HostMirror
        &specrefndxsw_host,  // complex refractive index for water visible
    mam_coupling::complex_view_2d::HostMirror &specrefndxlw_host)
{
  std::string sw_real_name = "refindex_real_aer_sw";
  std::string lw_real_name = "refindex_real_aer_lw";
  std::string sw_im_name   = "refindex_im_aer_sw";
  std::string lw_im_name   = "refindex_im_aer_lw";

  for(int i = 0; i < nswbands; i++) {
    auto sw_real_h = fields[sw_real_name].get_view<const Real*,Host>();
    auto sw_im_h   = fields[sw_im_name].get_view<const Real*,Host>();
    specrefndxsw_host(i, species_id).real() = sw_real_h(i);
    specrefndxsw_host(i, species_id).imag() = haero::abs(sw_im_h(i));
  }
  for(int i = 0; i < nlwbands; i++) {
    auto lw_real_h = fields[lw_real_name].get_view<const Real*,Host>();
    auto lw_im_h   = fields[lw_im_name].get_view<const Real*,Host>();
    specrefndxlw_host(i, species_id).real() = lw_real_h(i);
    specrefndxlw_host(i, species_id).imag() = haero::abs(lw_im_h(i));
  }

}  // copy_refindex_to_device

}  // namespace scream::mam_coupling

#endif
