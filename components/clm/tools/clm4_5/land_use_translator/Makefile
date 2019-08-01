CC=gcc
CFLAGS=-lnetcdf -lm -I$(INC_NETCDF) -L$(LIB_NETCDF)

land_use_translator: full_updateannuallanduse.c
	$(CC) -o land_use_translator full_updateannuallanduse.c $(CFLAGS)