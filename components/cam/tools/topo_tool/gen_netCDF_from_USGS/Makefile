#
# with PGI compiler "-byteswapio" option is needed to read in the binary USGS files correctly
#
EXEDIR = .
EXENAME = create_netCDF
RM = rm

.SUFFIXES:
.SUFFIXES: .F90 .o


FC = lf95

# Check for the NetCDF library and include directories 
 ifeq ($(LIB_NETCDF),$(null))
 LIB_NETCDF := /usr/local/lib
 endif
#
 ifeq ($(INC_NETCDF),$(null))
 INC_NETCDF := /usr/local/include
 endif


#------------------------------------------------------------------------
# AIX
# #------------------------------------------------------------------------
#
 ifeq ($(UNAMES),AIX)
 FC = xlf90
 FFLAGS = -c -I$(INC_NETCDF)
 LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
 .F90.o:
         $(FC) $(FFLAGS) -qsuffix=f=F90 $<
         endif

         .F90.o:
                 $(FC) $(FFLAGS) $<


#
# #------------------------------------------------------------------------
# # LF95
# #------------------------------------------------------------------------
#
 ifeq ($(FC),lf95)
 #
 # Tramhill
 #
   INC_NETCDF :=/usr/local/netcdf-4.1.3-gcc-4.4.4-13-lf9581/include
     LIB_NETCDF :=/usr/local/netcdf-4.1.3-gcc-4.4.4-13-lf9581/lib

       LDFLAGS = -L$(LIB_NETCDF) -lnetcdf -lnetcdff -lcurl -lhdf5 -lhdf5_hl -mcmodel=medium 
         FFLAGS   := -c --trace --trap --wide -CcdRR8 -I$(INC_NETCDF)
           ifeq ($(DEBUG),TRUE)
               #TBH:  this works    FFLAGS += -g --chk --pca
   #TBH:  this FAILS    FFLAGS += -g --chk a,e,s,u,x --pca
    FFLAGS += -g --chk a,e,s,u --pca
  else
    FFLAGS += -O
  endif

endif


.F90.o:
	$(FC) $(FFLAGS) $<

#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------

OBJS := create_netCDF_from_USGS.o shr_kind_mod.o

$(EXEDIR)/$(EXENAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	$(RM) -f $(OBJS)  *.mod $(EXEDIR)/$(EXENAME)

create_netCDF_from_USGS.o: shr_kind_mod.o
