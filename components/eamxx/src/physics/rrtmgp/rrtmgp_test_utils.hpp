#ifndef RRTMGP_TEST_UTILS_HPP
#define RRTMGP_TEST_UTILS_HPP

#include "cpp/extensions/cloud_optics/mo_cloud_optics.h"
#include "physics/rrtmgp/eamxx_rrtmgp_interface.hpp"
#include "cpp/rrtmgp/mo_gas_concentrations.h"
#include "cpp/rte/mo_fluxes.h"
#include "cpp/extensions/cloud_optics/mo_cloud_optics.h"

namespace rrtmgpTest {

bool file_exists(const char *filename);

#ifdef RRTMGP_ENABLE_YAKL
bool all_close(real2d &arr1, real2d &arr2, double tolerance);

void dummy_clouds(
  CloudOptics &cloud_optics, real2d &p_lay, real2d &t_lay,
  real2d &lwp, real2d &iwp, real2d &rel, real2d &rei, real2d &cld
);

void dummy_atmos(
  std::string inputfile,
  int ncol, real2d &p_lay, real2d &t_lay,
  real1d &sfc_alb_dir_vis, real1d &sfc_alb_dir_nir,
  real1d &sfc_alb_dif_vis, real1d &sfc_alb_dif_nir,
  real1d &mu0,
  real2d &lwp, real2d &iwp, real2d &rel, real2d &rei, real2d &cld
);

void read_fluxes(
  std::string inputfile,
  real2d &sw_flux_up, real2d &sw_flux_dn, real2d &sw_flux_dir,
  real2d &lw_flux_up, real2d &lw_flux_dn
);

void write_fluxes(
  std::string outputfile,
  real2d &sw_flux_up, real2d &sw_flux_dn, real2d &sw_flux_dir,
  real2d &lw_flux_up, real2d &lw_flux_dn
);
#endif

#ifdef RRTMGP_ENABLE_KOKKOS
template <typename RealT=scream::Real, typename LayoutT=Kokkos::LayoutRight, typename DeviceT=DefaultDevice>
struct rrtmgp_test_utils {

using interface_t = scream::rrtmgp::rrtmgp_interface<RealT, LayoutT, DeviceT>;
using real1dk = typename interface_t::view_t<RealT*>;
using real2dk = typename interface_t::view_t<RealT**>;
using real3dk = typename interface_t::view_t<RealT***>;
using MDRP = typename conv::MDRP<LayoutT>;

static bool all_close(real2dk &arr1, real2dk &arr2, double tolerance)
{
  int nx = arr1.extent(0);
  int ny = arr2.extent(1);
  auto arr1_h = Kokkos::create_mirror_view(arr1);
  auto arr2_h = Kokkos::create_mirror_view(arr2);
  Kokkos::deep_copy(arr1_h, arr1);
  Kokkos::deep_copy(arr2_h, arr2);
  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      if (abs(arr1_h(i,j) - arr2_h(i,j)) > tolerance || std::isnan(arr1_h(i,j) - arr2_h(i,j))) {
        printf("arr1 = %f, arr2 = %f at %i,%i\n", arr1_h(i,j), arr2_h(i,j), i, j);
        return false;
      }
    }
  }
  return true;
}

static void dummy_clouds(
  CloudOpticsK<RealT, LayoutT, DeviceT> &cloud_optics, real2dk &p_lay, real2dk &t_lay,
  real2dk &lwp, real2dk &iwp, real2dk &rel, real2dk &rei, real2dk &cloud_mask
)
{
  // Problem sizes
  int ncol = t_lay.extent(0);
  int nlay = t_lay.extent(1);

  // Generate some fake liquid and ice water data. We pick values to be midway between
  // the min and max of the valid lookup table values for effective radii
  real rel_val = 0.5 * (cloud_optics.get_min_radius_liq() + cloud_optics.get_max_radius_liq());
  real rei_val = 0.5 * (cloud_optics.get_min_radius_ice() + cloud_optics.get_max_radius_ice());

  // Restrict clouds to troposphere (> 100 hPa = 100*100 Pa) and not very close to the ground (< 900 hPa), and
  // put them in 2/3 of the columns since that's roughly the total cloudiness of earth.
  // Set sane values for liquid and ice water path.
  // NOTE: these "sane" values are in g/m2!
  Kokkos::parallel_for( MDRP::template get<2>({ncol, nlay}) , KOKKOS_LAMBDA (int icol, int ilay) {
    cloud_mask(icol,ilay) = p_lay(icol,ilay) > 100. * 100. && p_lay(icol,ilay) < 900. * 100. && ((icol+1)%3) != 0;
    // Ice and liquid will overlap in a few layers
    lwp(icol,ilay) = conv::merge(10.,  0., cloud_mask(icol,ilay) && t_lay(icol,ilay) > 263.);
    iwp(icol,ilay) = conv::merge(10.,  0., cloud_mask(icol,ilay) && t_lay(icol,ilay) < 273.);
    rel(icol,ilay) = conv::merge(rel_val, 0., lwp(icol,ilay) > 0.);
    rei(icol,ilay) = conv::merge(rei_val, 0., iwp(icol,ilay) > 0.);
  });
}

static void dummy_atmos(
  std::string inputfile,
  int ncol, real2dk &p_lay, real2dk &t_lay,
  real1dk &sfc_alb_dir_vis, real1dk &sfc_alb_dir_nir,
  real1dk &sfc_alb_dif_vis, real1dk &sfc_alb_dif_nir,
  real1dk &mu0,
  real2dk &lwp, real2dk &iwp, real2dk &rel, real2dk &rei, real2dk &cld
)
{
  // Setup boundary conditions, solar zenith angle, etc
  // NOTE: this stuff would come from the model in a real run

  // Ocean-ish values for surface albedos, just for example
  Kokkos::deep_copy(sfc_alb_dir_vis , 0.06 );
  Kokkos::deep_copy(sfc_alb_dir_nir , 0.06 );
  Kokkos::deep_copy(sfc_alb_dif_vis , 0.06 );
  Kokkos::deep_copy(sfc_alb_dif_nir , 0.06 );

  // Pick a solar zenith angle; this should come from the model
  Kokkos::deep_copy(mu0, 0.86 );

  // Get dummy cloud PHYSICAL properties. Note that this function call
  // needs the CloudOptics object only because it uses the min and max
  // valid values from the lookup tables for liquid and ice water path to
  // create a dummy atmosphere.
  dummy_clouds(*interface_t::cloud_optics_sw_k, p_lay, t_lay, lwp, iwp, rel, rei, cld);
}

static void read_fluxes(
  std::string inputfile,
  real2dk &sw_flux_up, real2dk &sw_flux_dn, real2dk &sw_flux_dir,
  real2dk &lw_flux_up, real2dk &lw_flux_dn
)
{
  // Initialize netcdf reader
  conv::SimpleNetCDF io;
  io.open(inputfile, NC_NOWRITE);

  // Initialize arrays to hold fluxes
  int nlev = io.getDimSize("lev");
  int ncol = io.getDimSize("col_flx");
  sw_flux_up  = real2dk("sw_flux_up" , ncol, nlev);
  sw_flux_dn  = real2dk("sw_flux_dn" , ncol, nlev);
  sw_flux_dir = real2dk("sw_flux_dir", ncol, nlev);
  lw_flux_up  = real2dk("lw_flux_up" , ncol, nlev);
  lw_flux_dn  = real2dk("lw_flux_dn" , ncol, nlev);

  // Read data
  io.read(sw_flux_up,  "sw_flux_up" );
  io.read(sw_flux_dn,  "sw_flux_dn" );
  io.read(sw_flux_dir, "sw_flux_dir");
  io.read(lw_flux_up,  "lw_flux_up" );
  io.read(lw_flux_dn,  "lw_flux_dn" );
}

static void write_fluxes(
  std::string outputfile,
  real2dk &sw_flux_up, real2dk &sw_flux_dn, real2dk &sw_flux_dir,
  real2dk &lw_flux_up, real2dk &lw_flux_dn
)
{
  conv::SimpleNetCDF io;
  io.create(outputfile);
  io.write(sw_flux_up , "sw_flux_up" , {"col_flx","lev"});
  io.write(sw_flux_dn , "sw_flux_dn" , {"col_flx","lev"});
  io.write(sw_flux_dir, "sw_flux_dir", {"col_flx","lev"});
  io.write(lw_flux_up , "lw_flux_up" , {"col_flx","lev"});
  io.write(lw_flux_dn , "lw_flux_dn" , {"col_flx","lev"});
  io.close();
}

};
#endif

}
#endif
