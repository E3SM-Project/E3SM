#include "mo_gas_concentrations.h"
#include "mo_fluxes.h"
#include "mo_cloud_optics.h"
namespace rrtmgpTest {
    bool all_equals(real2d &arr1, real2d &arr2);
    void dummy_clouds(
            CloudOptics &cloud_optics, real2d &p_lay, real2d &t_lay, 
            real2d &lwp, real2d &iwp, real2d &rel, real2d &rei
        );
    void dummy_atmos(
            std::string inputfile,
            int ncol, real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev, GasConcs &gas_concs, real2d &col_dry, 
            real2d &sfc_alb_dir, real2d &sfc_alb_dif, real1d &mu0,
            real2d &lwp, real2d &iwp, real2d &rel, real2d &rei
        );
    void read_fluxes(
            std::string inputfile, 
            real2d &sw_flux_up, real2d &sw_flux_dn, real2d &sw_flux_dn_dir,
            real2d &lw_flux_up, real2d &lw_flux_dn
        );
}
