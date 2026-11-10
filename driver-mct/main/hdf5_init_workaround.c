extern int H5open(void);

__attribute__((constructor))
static void init_hdf5_early(void) {
    (void)H5open();
}
