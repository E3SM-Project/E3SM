#ifndef RRTMGP_TEST_UTILS_HPP
#define RRTMGP_TEST_UTILS_HPP
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
bool all_close(real2dk &arr1, real2dk &arr2, double tolerance);

void dummy_clouds(
  CloudOpticsK &cloud_optics, real2dk &p_lay, real2dk &t_lay,
  real2dk &lwp, real2dk &iwp, real2dk &rel, real2dk &rei, real2dk &cld
);

void dummy_atmos(
  std::string inputfile,
  int ncol, real2dk &p_lay, real2dk &t_lay,
  real1dk &sfc_alb_dir_vis, real1dk &sfc_alb_dir_nir,
  real1dk &sfc_alb_dif_vis, real1dk &sfc_alb_dif_nir,
  real1dk &mu0,
  real2dk &lwp, real2dk &iwp, real2dk &rel, real2dk &rei, real2dk &cld
);

void read_fluxes(
  std::string inputfile,
  real2dk &sw_flux_up, real2dk &sw_flux_dn, real2dk &sw_flux_dir,
  real2dk &lw_flux_up, real2dk &lw_flux_dn
);

void write_fluxes(
  std::string outputfile,
  real2dk &sw_flux_up, real2dk &sw_flux_dn, real2dk &sw_flux_dir,
  real2dk &lw_flux_up, real2dk &lw_flux_dn
);
#endif

}
#endif
