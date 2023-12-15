#include "mo_load_coefficients.h"
#include "simple_netcdf.hpp"

// This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
//
// Contacts: Robert Pincus and Eli Mlawer
// email:  rrtmgp@aer.com
//
// Copyright 2015-2018,  Atmospheric and Environmental Research and
// Regents of the University of Colorado.  All right reserved.
//
// Use and duplication is permitted under the terms of the
//    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
// -------------------------------------------------------------------------------------------------

void load_and_init(GasOpticsRRTMGP &kdist, std::string filename, GasConcs const &available_gases) {
  simple_netcdf::SimpleNetCDF io;
  io.open(filename , NC_NOWRITE);

  // Read the many arrays
  string1d     gas_names;
  string1d     gas_minor;
  string1d     identifier_minor;
  string1d     minor_gases_lower;
  string1d     minor_gases_upper;
  string1d     scaling_gas_lower;
  string1d     scaling_gas_upper;
  intHost3d    key_species;
  realHost2d   band_lims;
  intHost2d    band2gpt;
  real         press_ref_trop;
  real         temp_ref_p;
  real         temp_ref_t;
  realHost1d   press_ref;
  realHost1d   temp_ref;
  realHost3d   vmr_ref;
  realHost4d   kmajor;
  intHost2d    minor_limits_gpt_lower;
  intHost2d    minor_limits_gpt_upper;
  boolHost1d   minor_scales_with_density_lower;
  boolHost1d   minor_scales_with_density_upper;
  boolHost1d   scale_by_complement_lower;
  boolHost1d   scale_by_complement_upper;
  intHost1d    kminor_start_lower;
  intHost1d    kminor_start_upper;
  realHost3d   kminor_lower;
  realHost3d   kminor_upper;
  realHost3d   rayl_lower;
  realHost3d   rayl_upper;

  // Read in strings
  charHost2d tmp;
  tmp = charHost2d();  io.read( tmp , "gas_names"         );  gas_names         = char2d_to_string1d(tmp);
  tmp = charHost2d();  io.read( tmp , "gas_minor"         );  gas_minor         = char2d_to_string1d(tmp);
  tmp = charHost2d();  io.read( tmp , "identifier_minor"  );  identifier_minor  = char2d_to_string1d(tmp);
  tmp = charHost2d();  io.read( tmp , "minor_gases_lower" );  minor_gases_lower = char2d_to_string1d(tmp);
  tmp = charHost2d();  io.read( tmp , "minor_gases_upper" );  minor_gases_upper = char2d_to_string1d(tmp);
  tmp = charHost2d();  io.read( tmp , "scaling_gas_lower" );  scaling_gas_lower = char2d_to_string1d(tmp);
  tmp = charHost2d();  io.read( tmp , "scaling_gas_upper" );  scaling_gas_upper = char2d_to_string1d(tmp);

  io.read( key_species                     , "key_species" );
  io.read( band_lims                       , "bnd_limits_wavenumber" );
  io.read( band2gpt                        , "bnd_limits_gpt" );
  io.read( press_ref                       , "press_ref" );
  io.read( temp_ref                        , "temp_ref" );
  io.read( temp_ref_p                      , "absorption_coefficient_ref_P" );
  io.read( temp_ref_t                      , "absorption_coefficient_ref_T" );
  io.read( press_ref_trop                  , "press_ref_trop" );
  io.read( kminor_lower                    , "kminor_lower" );
  io.read( kminor_upper                    , "kminor_upper" );
  io.read( minor_limits_gpt_lower          , "minor_limits_gpt_lower" );
  io.read( minor_limits_gpt_upper          , "minor_limits_gpt_upper" );
  io.read( minor_scales_with_density_lower , "minor_scales_with_density_lower" );
  io.read( minor_scales_with_density_upper , "minor_scales_with_density_upper" );
  io.read( scale_by_complement_lower       , "scale_by_complement_lower" );
  io.read( scale_by_complement_upper       , "scale_by_complement_upper" );
  io.read( kminor_start_lower              , "kminor_start_lower" );
  io.read( kminor_start_upper              , "kminor_start_upper" );
  io.read( vmr_ref                         , "vmr_ref" );
  io.read( kmajor                          , "kmajor" );

  if (io.varExists("rayl_lower")) {
    io.read( rayl_lower , "rayl_lower" );
    io.read( rayl_upper , "rayl_upper" );
  }

  // Initialize the gas optics class with data. The calls look slightly different depending
  //   on whether the radiation sources are internal to the atmosphere (longwave) or external (shortwave)
  // gas_optics%load() returns a string; a non-empty string indicates an error.
  if (io.varExists("totplnk")) {
    // If there's a totplnk variable in the file, then it's a longwave (internal sources) type
    realHost2d totplnk;
    realHost4d planck_frac;
    io.read( totplnk     , "totplnk"        );
    io.read( planck_frac , "plank_fraction" );
    kdist.load(available_gases, gas_names, key_species, band2gpt, band_lims, press_ref, press_ref_trop, 
               temp_ref, temp_ref_p, temp_ref_t, vmr_ref, kmajor, kminor_lower, kminor_upper, 
               gas_minor, identifier_minor, minor_gases_lower, minor_gases_upper, 
               minor_limits_gpt_lower, minor_limits_gpt_upper, minor_scales_with_density_lower, 
               minor_scales_with_density_upper, scaling_gas_lower, scaling_gas_upper, 
               scale_by_complement_lower, scale_by_complement_upper, kminor_start_lower, 
               kminor_start_upper, totplnk, planck_frac, rayl_lower, rayl_upper);
  } else {
    // Otherwise, it's a shortwave type
    realHost1d solar_src;
    try {
      if (io.varExists("solar_source")) {
        io.read( solar_src , "solar_source" );
      } else if (io.varExists("solar_source_quiet")) {
        // Newer RRTMGP input data files include components of solar source function sufficient to
        // compute solar variability; ignore this for now, and simply read in the baseline solar source.
        // This seems to get us closer to the original average insolation from
        // the full spectral resolution data, but there is disagreement between
        // this and EAM RRTMG for some reason.
        // TODO: Ultimately, we need to be able to use solar source data consistent with EAM (i.e., the
        // input4MIPS data); this will require implementing a separate function
        // here to read that data in and integrate the solar source over each
        // wavelength/gpoint band
        io.read( solar_src            , "solar_source_quiet" );
        if (false) {
          realHost1d solar_source_facular;
          realHost1d solar_source_sunspot;
          real mg_index;
          real sb_index;
          io.read( solar_source_facular , "solar_source_facular" );
          io.read( solar_source_sunspot , "solar_source_sunspot" );
          io.read( mg_index, "mg_default");
          io.read( sb_index, "sb_default");
          const real a_offset = 0.1495954;
          const real b_offset = 0.00066696;
          auto ngpt = solar_src.totElems();
          yakl::fortran::parallel_for(yakl::fortran::SimpleBounds<1>(ngpt), YAKL_LAMBDA(int igpt) {
              solar_src(igpt) = solar_src(igpt)
                  + (mg_index - a_offset) * solar_source_facular(igpt)
                  + (sb_index - b_offset) * solar_source_sunspot(igpt);
          });
        }
      } else {
        throw 1;
      }
    } catch (int err) {
      std::cout << "ERROR: no solar_source variable found in input data.\n";
    }
    kdist.load(available_gases, gas_names, key_species, band2gpt, band_lims, press_ref, press_ref_trop, 
               temp_ref, temp_ref_p, temp_ref_t, vmr_ref, kmajor, kminor_lower, kminor_upper, 
               gas_minor, identifier_minor, minor_gases_lower, minor_gases_upper, 
               minor_limits_gpt_lower, minor_limits_gpt_upper, minor_scales_with_density_lower, 
               minor_scales_with_density_upper, scaling_gas_lower, scaling_gas_upper, 
               scale_by_complement_lower, scale_by_complement_upper, kminor_start_lower, 
               kminor_start_upper, solar_src, rayl_lower, rayl_upper);
  }
  io.close();
}


