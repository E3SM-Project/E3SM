## This is the automake file to build the PIO C library.
# Ed Hartnett 8/19/17

# The library we are building.
lib_LTLIBRARIES = libpio.la

# These linker flags specify libtool version info.
# See http://www.gnu.org/software/libtool/manual/libtool.html#Libtool-versioning
# for information regarding incrementing `-version-info`.
libpio_la_LDFLAGS = -version-info 1:0:0

# The header file.
include_HEADERS = pio.h

# The library soure files.
libpio_la_SOURCES = bget.c pioc_sc.c pio_darray.c pio_file.c		\
pio_getput_int.c pio_msg.c pio_nc.c pio_rearrange.c pioc.c		\
pioc_support.c pio_darray_int.c pio_get_nc.c pio_lists.c pio_nc4.c	\
pio_put_nc.c pio_spmd.c pio_internal.h bget.h
