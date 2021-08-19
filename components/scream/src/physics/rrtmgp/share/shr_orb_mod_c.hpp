#ifndef SHR_ORB_MOD_C_HPP
#define SHR_ORB_MOD_C_HPP
extern "C" void shr_orb_params_c(
        int *orb_year , double *eccen, double *mvelpp, double *lambm0,
        double *obliqr, double *delta, double *eccf
        );
extern "C" void shr_orb_decl_c(
        double calday, double eccen, double mvelpp, double lambm0, 
        double obliqr, double delta, double eccf
        );
extern "C" double shr_orb_cosz_c(
        double jday, double lat, double lon, double declin, double dt_avg
        );
#endif
