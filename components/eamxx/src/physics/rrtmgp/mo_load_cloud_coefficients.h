#ifndef MO_LOAD_CLOUD_COEFFICIENTS_HPP
#define MO_LOAD_CLOUD_COEFFICIENTS_HPP

#include "cpp/rrtmgp_const.h"
#include "cpp/rte/mo_optical_props.h"
#include "cpp/extensions/cloud_optics/mo_cloud_optics.h"

void load_cld_lutcoeff(CloudOptics &cloud_spec, std::string cld_coeff_file);

void load_cld_padecoeff(CloudOptics &cloud_spec, std::string cld_coeff_file);

#endif


