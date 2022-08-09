#ifndef SHR_ORB_MOD_C2F_HPP
#define SHR_ORB_MOD_C2F_HPP
extern "C" int shr_orb_undef_int_c2f;
extern "C" void shr_orb_params_c2f(
        int *orb_year , double *eccen, double *mvelpp, double *lambm0,
        double *obliqr, double *delta, double *eccf
        );
extern "C" void shr_orb_decl_c2f(
        double calday, double eccen, double mvelpp, double lambm0,
        double obliqr, double *delta, double *eccf
        );
extern "C" double shr_orb_cosz_c2f(
        double jday, double lat, double lon, double declin, double dt_avg
        );
#endif
