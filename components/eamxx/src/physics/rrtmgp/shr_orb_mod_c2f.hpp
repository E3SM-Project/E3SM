#ifndef SHR_ORB_MOD_C2F_HPP
#define SHR_ORB_MOD_C2F_HPP
using shr_orb_decl_callback = bool (*)(int, double, double*, double*);
using shr_orb_geometry_callback = bool (*)(int, double, double*, double*, double*);
extern "C" void shr_orb_set_geometry_callback_c2f(shr_orb_geometry_callback callback, bool enabled);
extern "C" bool shr_orb_has_geometry_callback_c2f();
extern "C" bool shr_orb_use_equation_of_time_c2f();
extern "C" bool shr_orb_coupled_geometry_c2f(int ymd, double tod, double* delta, double* eccf, double* eot_seconds);
extern "C" void shr_orb_set_decl_callback_c2f(shr_orb_decl_callback callback);
extern "C" bool shr_orb_coupled_decl_c2f(int ymd, double tod, double* delta, double* eccf);
extern "C" int shr_orb_undef_int_c2f;
extern "C" void shr_orb_params_c2f(
        int *orb_year , double *eccen, double *mvelpp, double *lambm0,
        double *obliqr, double *delta, double *eccf
        );
extern "C" void shr_orb_params_continuous_c2f(
        int orb_year, double year_fraction, double *eccen, double *obliq,
        double *mvelp, double *obliqr, double *lambm0, double *mvelpp
        );
extern "C" void shr_orb_decl_c2f(
        double calday, double eccen, double mvelpp, double lambm0,
        double obliqr, double *delta, double *eccf
        );
extern "C" void shr_orb_decl_gregorian_c2f(
        double days_since_20000101, double eccen, double mvelpp, double lambm0,
        double obliqr, double *delta, double *eccf
        );
extern "C" void shr_orb_decl_fixed_year_c2f(
        double calday, double days_per_year, double eccen, double mvelpp,
        double lambm0, double obliqr, double* delta, double* eccf,
        double* eot_seconds, bool use_eot);
extern "C" double shr_orb_cosz_c2f(
        double jday, double lat, double lon, double declin, double dt_avg
        );
extern "C" void shr_orb_decl_eot_c2f(
        double calday, double eccen, double mvelpp, double lambm0,
        double obliqr, double* delta, double* eccf, double* eot_seconds);
extern "C" void shr_orb_decl_gregorian_eot_c2f(
        double elapsed_days, double eccen, double mvelpp, double lambm0,
        double obliqr, double* delta, double* eccf, double* eot_seconds);
extern "C" double shr_orb_cosz_eot_c2f(
        double jday, double lat, double lon, double declin, double dt_avg, double eot_seconds);
#endif
