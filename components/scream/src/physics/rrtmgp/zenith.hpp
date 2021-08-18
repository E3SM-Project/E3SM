#ifndef ZENITH_HPP
#define ZENITH_HPP

/*
 * Okay, things we do in EAM:
 * get_cosine_solar_zenith -> zenith -> shr_orb_decl -> shr_orb_avg_cosz
 * But we need orbital parameters from shr_orb_params in shr_orb_mod.
 * It looks to me like we could start by compiling the fortran source
 * shr_orb_mod, and then providing C++ bridges to shr_orb_decl,
 * shr_orb_avg_cosz, and shr_orb_params. We could then have a zenith
 * function that would call, in order:
 *
 *     shr_orb_params()
 *     shr_orb_decl()
 *     shr_orb_avg_cosz()
 *
 * That should at least get us a baseline, and then we can work on porting
 * shr_orb_mod.F90 to C++.
 */

extern void shr_orb_params(double calday, double eccen, double mvelpp, double lambm0,
                           double obliqr, double delta, double eccf);

extern void shr_orb_decl(double calday, double eccen, double mvelpp, double lambm0, 
                         double obliqr, double delta, double eccf);

extern double shr_orb_avg_cosz(double jday, double lat, double lon, double declin, double dt_avg);

// TODO: port things to C++
// Function to compute solar declination given:
//   calday (calendar day)
//   eccen  ()
//   mvelpp
//   lambm0
//   obliqr
//   delta
//   eccf 
//template<class T> T shr_orb_decl(T calday, T eccen, T mvelpp, T lambm0, T obliqr, T delta, T eccf) {
//}

// Function to compute cosine of solar zenith angle given julian day (jday),
// latitude (lat), longitude (lon), solar declination angle (declin), and
// timestep for average cosz computation (dt_avg). This is a port of the
// like-named function in e3sm/share/util/shr_orb_mod.F90. According to the
// comments therein, that function has the following credits:
//
//!=======================================================================
//! A New Algorithm for Calculation of Cosine Solar Zenith Angle
//! Author: Linjiong Zhou
//! E-mail: linjiongzhou@hotmail.com
//! Date  : 2015.02.22
//! Ref.  : Zhou et al., GRL, 2015
//!=======================================================================
//template<class T> T shr_orb_avg_cosz(T jday, T lat, T lon, T declin, T dt_avg) {
//    T cosz = 0;
//    return cosz;
//};
#endif
