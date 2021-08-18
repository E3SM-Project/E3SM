#ifndef SHR_ORB_MOD_HPP
#define SHR_ORB_MOD_HPP
extern void shr_orb_params(
        double calday, double eccen, double mvelpp, double lambm0,
        double obliqr, double delta, double eccf
        );
extern void shr_orb_decl(
        double calday, double eccen, double mvelpp, double lambm0, 
        double obliqr, double delta, double eccf
        );
extern double shr_orb_avg_cosz(
        double jday, double lat, double lon, double declin, double dt_avg
        );
#endif
