#ifndef __PIO_INTERNAL__
#define __PIO_INTERNAL__
#include <pio.h>

file_desc_t *pio_get_file_from_id(int ncid);
int pio_delete_file_from_list(int ncid);
void pio_add_to_file_list(file_desc_t *file);
iosystem_desc_t *pio_get_iosystem_from_id(int iosysid);
int pio_add_to_iosystem_list(iosystem_desc_t *ios);
int check_netcdf(file_desc_t *file,const int status, const char *fname, const int line);
int iotype_error(const int iotype, const char *fname, const int line);
#endif
