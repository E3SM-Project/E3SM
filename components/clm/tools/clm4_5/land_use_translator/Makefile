# Set user specified C compiler
ifneq ($(USER_CC),$(null))
    CC := $(USER_CC)
endif

# Check for the netcdf library and include directories
ifeq ($(LIB_NETCDF),$(null))
  LIB_NETCDF := /usr/local/lib
endif

ifeq ($(INC_NETCDF),$(null))
  INC_NETCDF := /usr/local/include
endif

CFLAGS=-lnetcdf -lm -I$(INC_NETCDF) -L$(LIB_NETCDF)

land_use_translator: full_updateannuallanduse.c
	$(CC) -o land_use_translator full_updateannuallanduse.c $(CFLAGS)

clean:
	rm -f land_use_translator