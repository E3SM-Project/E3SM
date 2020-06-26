# Set user specified C compiler
ifneq ($(USER_CC),$(null))
    CC := $(USER_CC)
endif

HN := $(shell hostname)
echo$(HN)

# compiler at NERSC knows where netcdf is because the prg env is set
ifeq ($(HN),cori09)

ifeq ($(LIB_NETCDF),$(null))
  LIB_NETCDF := /usr/local/lib
endif

ifeq ($(INC_NETCDF),$(null))
  INC_NETCDF := /usr/local/include
endif

endif

# compilers at COMPY does not know where netcdf is
# set this for the default gcc/4.8.5 compiler
# if you want to use a different compiler you need to change this
ifeq ($(HN),compy-e.pnl.gov)

LIB_NETCDF := /share/apps/netcdf/4.6.3/gcc/4.8.5/lib

INC_NETCDF := /share/apps/netcdf/4.6.3/gcc/4.8.5/include

endif


CFLAGS=-lnetcdf -lm -I$(INC_NETCDF) -L$(LIB_NETCDF)

land_use_translator: updateannuallanduse_v2.c
	$(CC) -std=c11 -O0 -o land_use_translator updateannuallanduse_v2.c $(CFLAGS)

clean:
	rm -f land_use_translator
