void initglm_() {
    initialize();
}
void runglm_() {
    run_model();
}
void stepglm_(int *curryear,double *glmi,int* glmi_fdim1, int* glmi_fdim2, double *glmo, int* glmo_fdim1, int* glmo_fdim2) {
  stepglm_ccsm(curryear,glmi,glmi_fdim1,glmi_fdim2,glmo,glmo_fdim1,glmo_fdim2);
}
void finalizeglm_() {
    finalize();
}
