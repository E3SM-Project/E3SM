#include "../../ccsmcpl/include/CCSM_GCAM_interface.h"

CCSM_GCAM_interface *p_obj;

extern "C" {

  void initccsminterface_() {
    p_obj = new CCSM_GCAM_interface();
  }
  void deleteccsminterface_() {
    delete p_obj;
  }
  void initcgcam_() {
    p_obj->initGCAM();
  }
  //  void runcgcam_(int *ymd, int *tod, double x2g[19][2], double g2x[19][3]) {
  //    p_obj->runGCAM(ymd, tod, x2g, g2x);
  //  }
  // new call passing size of double arrays
  void setdensitycgcam_(int *ymd, int *tod, double *gcami, int *gcami_fdim_1, int *gcami_fdim_2) {
    p_obj->setDensityGCAM(ymd, tod, gcami,gcami_fdim_1,gcami_fdim_2);
  }
  // new call passing size of double arrays
  void runcgcam_(int *ymd, int *tod, double *gcami, int *gcami_fdim_1, int *gcami_fdim_2, double *gcamo, int *gcamo_fdim_1, int *gcamo_fdim_2, double *gcamoemis, int *gcamoemis_fdim_1, int *gcamoemis_fdim_2,int *gcamoyr1, int* gcamoyr2,int* sneakermode,int* write_rest) {
    p_obj->runGCAM(ymd, tod, gcami,gcami_fdim_1,gcami_fdim_2,gcamo, gcamo_fdim_1, gcamo_fdim_2,gcamoemis, gcamoemis_fdim_1, gcamoemis_fdim_2,gcamoyr1, gcamoyr2,sneakermode,write_rest);
  }
  void finalizecgcam_() {
    p_obj->finalizeGCAM();
  }
}
