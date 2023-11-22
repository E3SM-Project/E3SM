#ifndef MAM_AEROSOL_OPTICS_READ_TABLES_HPP
#define MAM_AEROSOL_OPTICS_READ_TABLES_HPP

#include "share/io/scream_scorpio_interface.hpp"
#include "share/field/field_manager.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/grid/grids_manager.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "share/io/scorpio_input.hpp"

// NOTE: I will add functions for aerosol_optics here, I will move this code later to mam_coupling.hpp
namespace scream::mam_coupling {

using view_1d_host = typename KT::view_1d<Real>::HostMirror;
using view_2d_host = typename KT::view_2d<Real>::HostMirror;
using view_5d_host = Kokkos::View<Real *****>::HostMirror;

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

struct AerosolOpticsDeviceData {
// devices views
// FIXME: move this code to mam4x and simplify number of inputs for aerosol_optics
// FIXME: add description of these tables.
view_1d refitabsw[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nswbands];
view_1d refrtabsw[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nswbands];
view_1d refrtablw[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nlwbands];
view_1d refitablw[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nlwbands];

view_3d abspsw[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nswbands];
view_3d absplw[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nlwbands];
view_3d asmpsw[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nswbands];
view_3d extpsw[mam4::AeroConfig::num_modes()][mam4::modal_aer_opt::nswbands];
};



inline
void set_parameters_table(AerosolOpticsHostData &aerosol_optics_host_data,
                          ekat::ParameterList& rrtmg_params,
                          std::map<std::string,FieldLayout>&  layouts,
                          std::map<std::string,view_1d_host>& host_views
                          )
{

    // Set up input structure to read data from file.
  using strvec_t = std::vector<std::string>;
  using namespace ShortFieldTagsNames;

  constexpr int refindex_real = mam4::modal_aer_opt::refindex_real;
  constexpr int refindex_im = mam4::modal_aer_opt::refindex_im;
  constexpr int coef_number = mam4::modal_aer_opt::coef_number;

  auto refindex_real_lw_host = view_2d_host("refrtablw_real_host", nlwbands, refindex_real);
  auto refindex_im_lw_host = view_2d_host("refrtablw_im_host", nlwbands, refindex_im);

  auto refindex_real_sw_host = view_2d_host("refrtabsw_real_host", nswbands, refindex_real);
  auto refindex_im_sw_host = view_2d_host("refrtabsw_im_host", nswbands, refindex_im);

  // absplw(lw_band, mode, refindex_im, refindex_real, coef_number)
  auto absplw_host = view_5d_host("absplw_host", nlwbands, 1,
   refindex_im, refindex_real,  coef_number);

  auto asmpsw_host = view_5d_host("asmpsw_host", nswbands, 1, refindex_im, refindex_real, coef_number);
  auto extpsw_host = view_5d_host("extpsw_host", nswbands, 1, refindex_im, refindex_real, coef_number);
  auto abspsw_host = view_5d_host("abspsw_host", nswbands, 1, refindex_im, refindex_real, coef_number);

  aerosol_optics_host_data.refindex_real_lw_host = refindex_real_lw_host;
  aerosol_optics_host_data.refindex_im_lw_host =refindex_im_lw_host;
  aerosol_optics_host_data.refindex_real_sw_host =refindex_real_sw_host;
  aerosol_optics_host_data.refindex_im_sw_host = refindex_im_sw_host;
  aerosol_optics_host_data.absplw_host = absplw_host;
  aerosol_optics_host_data.asmpsw_host = asmpsw_host;
  aerosol_optics_host_data.extpsw_host = extpsw_host;
  aerosol_optics_host_data.abspsw_host = abspsw_host;

  FieldLayout scalar_refindex_real_lw_layout { {LWBAND, NREFINDEX_REAL},
                                     {nlwbands,
                                      refindex_real} };

  FieldLayout scalar_refindex_im_lw_layout { {LWBAND, NREFINDEX_IM},
                                     {nlwbands,
                                      refindex_im} };
  FieldLayout scalar_refindex_real_sw_layout { {SWBAND, NREFINDEX_REAL},
                                     {nswbands,
                                      refindex_real} };

  FieldLayout scalar_refindex_im_sw_layout { {SWBAND, NREFINDEX_IM},
                                     {nswbands,
                                      refindex_im} };

  FieldLayout scalar_absplw_layout { {LWBAND, MODE, NREFINDEX_IM, NREFINDEX_REAL, NCOEF_NUMBER},
                                     {nlwbands, 1, refindex_im, refindex_real,  coef_number} };
  // use also for extpsw, abspsw
  FieldLayout scalar_asmpsw_layout { {SWBAND, MODE, NREFINDEX_IM, NREFINDEX_REAL, NCOEF_NUMBER},
                                     { nswbands, 1, refindex_im, refindex_real, coef_number} };

  rrtmg_params.set<strvec_t>("Field Names",{
                                      "asmpsw",// need aditinal work
                                      "extpsw", // need aditinal work
                                      "abspsw",// need aditinal work
                                      "absplw", // need aditinal work
                                      "refindex_real_sw", // done copy
                                      "refindex_im_sw", // done copy
                                      "refindex_real_lw", // done copy
                                      "refindex_im_lw" // done copy
                                      });

  rrtmg_params.set("Skip_Grid_Checks",true);

  host_views["refindex_real_sw"]
   = view_1d_host(refindex_real_sw_host.data(),
   refindex_real_sw_host.size());

   host_views["refindex_im_sw"]
   = view_1d_host(refindex_im_sw_host.data(),
   refindex_im_sw_host.size());

  host_views["refindex_real_lw"]
   = view_1d_host(refindex_real_lw_host.data(),
   refindex_real_lw_host.size());

   host_views["refindex_im_lw"]
   = view_1d_host(refindex_im_lw_host.data(),
   refindex_im_lw_host.size());

   host_views["absplw"]
   = view_1d_host(absplw_host.data(),
   absplw_host.size());

    host_views["asmpsw"]
   = view_1d_host(asmpsw_host.data(),
   asmpsw_host.size());

    host_views["extpsw"]
   = view_1d_host(extpsw_host.data(),
   extpsw_host.size());

    host_views["abspsw"]
   = view_1d_host(abspsw_host.data(),
   abspsw_host.size());

  layouts.emplace("refindex_real_lw", scalar_refindex_real_lw_layout);
  layouts.emplace("refindex_im_lw", scalar_refindex_im_lw_layout);
  layouts.emplace("refindex_real_sw", scalar_refindex_real_sw_layout);
  layouts.emplace("refindex_im_sw", scalar_refindex_im_sw_layout);
  layouts.emplace("absplw", scalar_absplw_layout);
  layouts.emplace("asmpsw", scalar_asmpsw_layout);
  layouts.emplace("extpsw", scalar_asmpsw_layout);
  layouts.emplace("abspsw", scalar_asmpsw_layout);

}
// KOKKOS_INLINE_FUNCTION
inline
void read_rrtmg_table(const std::string& table_filename,
                     const int imode,
                     ekat::ParameterList& params,
                     const std::shared_ptr<const AbstractGrid>& grid,
                     const std::map<std::string,view_1d_host>& host_views_1d,
                     const std::map<std::string,FieldLayout>&  layouts,
                     const AerosolOpticsHostData &aerosol_optics_host_data,
                     const AerosolOpticsDeviceData &aerosol_optics_device_data ){


  constexpr int refindex_real = mam4::modal_aer_opt::refindex_real;
  constexpr int refindex_im = mam4::modal_aer_opt::refindex_im;
  constexpr int coef_number = mam4::modal_aer_opt::coef_number;


  params.set("Filename",table_filename);
  AtmosphereInput rrtmg(params, grid, host_views_1d, layouts);
  // // -1000 forces the interface to read a dataset that does not have time as variable.
  rrtmg.read_variables(-1000);
  rrtmg.finalize();

    // copy data from host to device for mode 1
  int d1=imode;
  for (int d3 = 0; d3 < nswbands; ++d3) {
    auto real_host_d3 = Kokkos::subview(aerosol_optics_host_data.refindex_real_sw_host,
     d3, Kokkos::ALL());
    Kokkos::deep_copy(aerosol_optics_device_data.refrtabsw[d1][d3],real_host_d3);
    auto im_host_d3 = Kokkos::subview(aerosol_optics_host_data.refindex_im_sw_host, d3, Kokkos::ALL());
    Kokkos::deep_copy(aerosol_optics_device_data.refitabsw[d1][d3],im_host_d3);
  } // d3

  for (int d3 = 0; d3 < nlwbands; ++d3) {
    auto real_host_d3 = Kokkos::subview(aerosol_optics_host_data.refindex_real_lw_host, d3, Kokkos::ALL());
    Kokkos::deep_copy(aerosol_optics_device_data.refrtablw[d1][d3],real_host_d3);
    auto im_host_d3 = Kokkos::subview(aerosol_optics_host_data.refindex_im_lw_host, d3, Kokkos::ALL());
    Kokkos::deep_copy(aerosol_optics_device_data.refitablw[d1][d3],im_host_d3);
  } // d3

  // NOTE: we need to reorder dimenstions in absplw
  // netcfd : (lw_band, mode, refindex_im, refindex_real, coef_number)
  // mam4xx : (mode, lw_band, coef_number, refindex_real, refindex_im )
  // e3sm : (ntot_amode,coef_number,refindex_real,refindex_im,nlwbands)
  // FIXME: it maybe not work in gpus.
  for (int d5 = 0; d5 < nlwbands; ++d5) {
  for (int d2 = 0; d2 < coef_number; d2++)
  for (int d3 = 0; d3 < refindex_real; d3++)
  for (int d4 = 0; d4 < refindex_im; d4++)
        aerosol_optics_device_data.absplw[d1][d5](d2,d3,d4) = aerosol_optics_host_data.absplw_host(d5,0,d4,d3,d2);
  }// d5

  // asmpsw, abspsw, extpsw
  // netcfd : (sw_band, mode, refindex_im, refindex_real, coef_number)
  // mam4xx : (mode, sw_band, coef_number, refindex_real, refindex_im )

  for (int d5 = 0; d5 < nswbands; ++d5)
  for (int d2 = 0; d2 < coef_number; d2++)
  for (int d3 = 0; d3 < refindex_real; d3++)
  for (int d4 = 0; d4 < refindex_im; d4++)
  {
    aerosol_optics_device_data.asmpsw[d1][d5](d2,d3,d4) = aerosol_optics_host_data.asmpsw_host(d5,0,d4,d3,d2);
    aerosol_optics_device_data.abspsw[d1][d5](d2,d3,d4) = aerosol_optics_host_data.abspsw_host(d5,0,d4,d3,d2);
    aerosol_optics_device_data.extpsw[d1][d5](d2,d3,d4) = aerosol_optics_host_data.extpsw_host(d5,0,d4,d3,d2);
  } // d5

}

inline
void read_water_refindex(const std::string& table_filename,
                     const std::shared_ptr<const AbstractGrid>& grid,
                     Kokkos::complex<Real> crefwlw[mam4::modal_aer_opt::nlwbands],
                     Kokkos::complex<Real> crefwsw[mam4::modal_aer_opt::nswbands])
  {
  // refractive index for water read in read_water_refindex
  // crefwsw(nswbands) ! complex refractive index for water visible
  // crefwlw(nlwbands) ! complex refractive index for water infrared

  using namespace ShortFieldTagsNames;
  using view_1d_host = typename KT::view_1d<Real>::HostMirror;
  // Set up input structure to read data from file.
  using strvec_t = std::vector<std::string>;

  // here a made a list of variables that I want to read from netcdf files
  ekat::ParameterList params;
  params.set("Filename",table_filename);
  params.set("Skip_Grid_Checks",true);

  params.set<strvec_t>("Field Names",{
    "refindex_im_water_lw",
    "refindex_im_water_sw",
    "refindex_real_water_lw",
    "refindex_real_water_sw"});
  // make a list of host views
  std::map<std::string,view_1d_host> host_views_water;
  // fist allocate host views.
  view_1d_host refindex_im_water_sw_host("refindex_im_water_sw_host", nswbands);
  view_1d_host refindex_real_water_sw_host("refindex_real_water_sw_host", nswbands);
  view_1d_host refindex_im_water_lw_host("refindex_im_water_lw_host", nlwbands);
  view_1d_host refindex_real_water_lw_host("refindex_real_water_lw_host", nlwbands);

  host_views_water["refindex_im_water_sw"]
   = refindex_im_water_sw_host;
   host_views_water["refindex_real_water_sw"]
   = refindex_real_water_sw_host;
  host_views_water["refindex_im_water_lw"]
   = refindex_im_water_lw_host;
   host_views_water["refindex_real_water_lw"]
   = refindex_real_water_lw_host;

  // defines layouts
  std::map<std::string,FieldLayout>  layouts_water;
  FieldLayout scalar_refindex_water_sw_layout { {SWBAND},
                                                   {nswbands} };
  FieldLayout scalar_refindex_water_lw_layout { {LWBAND},
                                                   {nlwbands} };

  layouts_water.emplace("refindex_im_water_sw",
  scalar_refindex_water_sw_layout);
  layouts_water.emplace("refindex_real_water_sw",
  scalar_refindex_water_sw_layout);
  layouts_water.emplace("refindex_im_water_lw",
  scalar_refindex_water_lw_layout);
  layouts_water.emplace("refindex_real_water_lw",
  scalar_refindex_water_lw_layout);

  // create a object to read data
  AtmosphereInput refindex_water(params, grid,
  host_views_water, layouts_water);
  refindex_water.read_variables();
  refindex_water.finalize();

  //move data to device
  // FIXME: check this correct or can be done in a better way.
  // maybe make a 1D vied of Kokkos::complex<Real>
  for (int i = 0; i < nlwbands; ++i)
  {
    // Kokkos::complex<Real> temp;
    crefwlw[i].real() = refindex_real_water_lw_host(i);
    crefwlw[i].imag() = refindex_im_water_lw_host(i);
    // Kokkos::deep_copy(crefwlw[i],temp);
  }

  for (int i = 0; i < nswbands; ++i)
  {
    // Kokkos::complex<Real> temp;
    crefwsw[i].real() = refindex_real_water_sw_host(i);
    crefwsw[i].imag() = refindex_im_water_sw_host(i);
    // Kokkos::deep_copy(crefwsw[i],temp);
  }

  }
//read_refindex_aero

inline
void set_refindex(std::string surname,
                  ekat::ParameterList& params,
                  std::map<std::string,view_1d_host>& host_views,
                  std::map<std::string,FieldLayout>&  layouts )
{
  // set variables names
  using view_1d_host = typename KT::view_1d<Real>::HostMirror;
  using strvec_t = std::vector<std::string>;
  using namespace ShortFieldTagsNames;

  std::string refindex_real_sw = "refindex_real_"+surname+"_sw";
  std::string refindex_im_sw = "refindex_im_"+surname+"_sw";
  std::string refindex_real_lw = "refindex_real_"+surname+"_lw";
  std::string refindex_im_lw = "refindex_im_"+surname+"_lw";

  params.set("Skip_Grid_Checks",true);
  params.set<strvec_t>("Field Names",{ refindex_real_sw,
                                       refindex_im_sw,
                                       refindex_real_lw,
                                       refindex_im_lw
                                       });
  // allocate host views
  host_views[refindex_real_sw]   =  view_1d_host(refindex_real_sw, nswbands);
  host_views[refindex_im_sw]     =  view_1d_host(refindex_im_sw, nswbands);
  host_views[refindex_real_lw]   = view_1d_host(refindex_real_lw, nlwbands);
  host_views[refindex_im_lw]     = view_1d_host(refindex_im_lw, nlwbands);

  FieldLayout scalar_refindex_sw_layout { {SWBAND}, {nswbands} };
  FieldLayout scalar_refindex_lw_layout { {LWBAND}, {nlwbands} };

  layouts.emplace(refindex_real_sw,
  scalar_refindex_sw_layout);
  layouts.emplace(refindex_im_sw,
  scalar_refindex_sw_layout);
  layouts.emplace(refindex_real_lw,
  scalar_refindex_lw_layout);
  layouts.emplace(refindex_im_lw,
  scalar_refindex_lw_layout);

}  //set_refindex_aero

inline
void copy_refindex_aerosol_to_device(const int species_id,
std::map<std::string,view_1d_host>& host_views,
mam_coupling::complex_view_2d  &specrefndxsw, // complex refractive index for water visible
  mam_coupling::complex_view_2d &specrefndxlw)
{
  std::string sw_real_name = "refindex_real_aer_sw";
  std::string lw_real_name = "refindex_real_aer_lw";
  std::string sw_im_name = "refindex_im_aer_sw";
  std::string lw_im_name = "refindex_im_aer_lw";

  // using complex_view_2d_host = typename KT::view_1d<Kokkos::complex<Real>>::HostMirror;
  // FIXME: can we do this copy without making another host view?
  // Here I assume that specrefndxsw and specrefndxlw are already allocated.
  auto specrefndxsw_host = Kokkos::create_mirror_view(specrefndxsw);
  auto specrefndxlw_host = Kokkos::create_mirror_view(specrefndxlw);

  for (int i = 0; i < nswbands; i++)
  {
    specrefndxsw_host(i,species_id).real() = host_views[sw_real_name](i);
    specrefndxsw_host(i,species_id).imag() = host_views[sw_im_name](i);
  }
  for (int i = 0; i < nlwbands; i++)
  {
    specrefndxlw_host(i,species_id).real() = host_views[lw_real_name](i);
    specrefndxlw_host(i,species_id).imag() = host_views[lw_im_name](i);
  }
  Kokkos::deep_copy(specrefndxsw, specrefndxsw_host);
  Kokkos::deep_copy(specrefndxlw, specrefndxlw_host);
} // copy_refindex_to_device

} // namespace scream::mam_coupling

#endif
