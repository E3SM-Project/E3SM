#pragma once

#include "const.h"
#include "mo_optical_props.h"
#include "mo_cloud_optics.h"


void load_cld_lutcoeff(CloudOptics &cloud_spec, std::string cld_coeff_file);


void load_cld_padecoeff(CloudOptics &cloud_spec, std::string cld_coeff_file);



