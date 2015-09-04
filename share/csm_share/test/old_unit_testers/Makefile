#-----------------------------------------------------------------------
# This Makefile is for doing csm_share unit testing
#------------------------------------------------------------------------
cpp_dirs := . ../shr ../../utils/mct/mct \
	  ../../utils/mct/mpeu ../../utils/esmf_wrf_timemgr ../../utils/timing \
	  ../../drv/shr
ifneq ($(SPMD),TRUE)
cpp_dirs += ../../utils/mct/mpi-serial
endif
cpp_dirs += ../../utils/pio
cpp_path := $(foreach dir,$(cpp_dirs),-I$(dir)) # format for command line

# Platform specific macros
include make.Macros

space := $(null) $(null)

ifneq ($(ESMF_BLD),$(null))
cpp_dirs += $(ESMF_LIB)
endif

# Expand any tildes in directory names.
VPATH    := $(foreach dir,$(cpp_dirs),$(wildcard $(dir)))
# Change spaces to colons.
VPATH    := $(subst $(space),:,$(VPATH))

.PHONY: debug clean all

all: test_shr_tInterp


ifneq ($(SPMD),TRUE)
   OBJS_NOMPI   := fort.o group.o collective.o comm.o list.o handles.o mpi.o recv.o req.o \
	      		send.o time.o
$(OBJS_NOMPI) shr_mpi_mod.o: mpif.h
mpif.h:
	ln -s ../../utils/mct/mpi-serial/mpif.real4double8.h $@
else
   OBJS_NOMPI   := $(null)
endif
OBJS_MCT := m_Accumulator.o m_AccumulatorComms.o m_AttrVect.o m_AttrVectComms.o \
            m_AttrVectReduce.o m_ConvertMaps.o m_ExchangeMaps.o m_GeneralGrid.o \
            m_GeneralGridComms.o m_GlobalMap.o m_GlobalSegMap.o \
            m_GlobalSegMapComms.o m_GlobalToLocal.o m_MCTWorld.o m_MatAttrVectMul.o \
            m_Merge.o m_Navigator.o m_Rearranger.o m_Router.o m_SparseMatrix.o \
            m_SparseMatrixComms.o m_SparseMatrixDecomp.o m_SparseMatrixPlus.o \
            m_SparseMatrixToMaps.o m_SpatialIntegral.o m_SpatialIntegralV.o \
            m_Transfer.o \
	    m_FcComms.o m_FileResolv.o m_Filename.o m_IndexBin_char.o \
	    m_IndexBin_integer.o m_IndexBin_logical.o m_List.o m_MergeSorts.o \
	    m_Permuter.o m_SortingTools.o m_StrTemplate.o m_String.o m_TraceBack.o \
	    m_chars.o m_die.o m_dropdead.o m_flow.o m_inpak90.o m_ioutil.o m_mall.o \
	    m_mpif.o m_mpif90.o m_mpout.o m_rankMerge.o m_realkinds.o m_stdio.o \
	    m_zeit.o get_zeits.o
OBJS_PIO := alloc_mod.o box_rearrange.o calcdisplace_mod.o iompi_mod.o \
	    ionf_mod.o \
            nf_mod.o pio.o pio_kinds.o pio_mpi_utils.o pio_nf_utils.o \
            pio_msg_callbacks.o pio_msg_getput_callbacks.o pio_msg_mod.o \
	    pio_nf_utils.o pio_quicksort.o pio_spmd_utils.o pio_support.o pio_types.o \
            pio_utils.o piodarray.o piolib_mod.o pionfatt_mod.o \
            pionfget_mod.o pionfput_mod.o pionfread_mod.o pionfwrite_mod.o \
            rearrange.o
OBJS_TIM   := perf_mod.o perf_utils.o GPTLget_memusage.o GPTLprint_memusage.o \
	      GPTLutil.o f_wrappers.o gptl.o gptl_papi.o threadutil.o

OBJS       := test_shr_sys.o shr_sys_mod.o shr_kind_mod.o shr_mpi_mod.o shr_const_mod.o shr_log_mod.o \
	      $(OBJS_NOMPI)
OBJS_FILE  := test_shr_file.o shr_sys_mod.o shr_kind_mod.o shr_file_mod.o shr_mpi_mod.o shr_log_mod.o \
	      $(OBJS_NOMPI)
OBJS_ORB := test_shr_orb.o shr_sys_mod.o shr_kind_mod.o shr_orb_mod.o shr_mpi_mod.o shr_log_mod.o \
	    shr_const_mod.o $(OBJS_NOMPI)
OBJS_STRMS  := test_shr_streams.o shr_kind_mod.o shr_stream_mod.o shr_sys_mod.o \
               shr_file_mod.o shr_string_mod.o shr_timer_mod.o shr_mpi_mod.o    \
               shr_cal_mod.o shr_ncread_mod.o shr_const_mod.o    \
               shr_log_mod.o test_mod.o $(OBJS_NOMPI)
OBJS_SCAM := test_shr_scam.o shr_strdata_mod.o shr_const_mod.o shr_kind_mod.o \
               shr_log_mod.o shr_sys_mod.o shr_file_mod.o shr_stream_mod.o \
               shr_map_mod.o shr_string_mod.o shr_cal_mod.o shr_orb_mod.o \
               shr_tinterp_mod.o shr_dmodel_mod.o shr_mct_mod.o mct_mod.o \
               perf_mod.o pio.o shr_mpi_mod.o seq_flds_mod.o shr_ncread_mod.o \
               shr_scam_mod.o shr_pcdf_mod.o shr_mct_mod.o mct_mod.o shr_timer_mod.o \
	       seq_drydep_mod.o test_mod.o \
               $(OBJS_NOMPI) $(OBJS_MCT) $(OBJS_PIO) $(OBJS_TIM)
OBJS_STIN   := test_shr_tInterp.o shr_kind_mod.o shr_const_mod.o shr_sys_mod.o \
               shr_string_mod.o shr_cal_mod.o shr_log_mod.o shr_orb_mod.o test_mod.o \
               shr_tInterp_mod.o shr_timer_mod.o shr_mpi_mod.o $(OBJS_NOMPI)
OBJS_MPI := test_shr_mpi.o shr_mpi_mod.o shr_kind_mod.o shr_sys_mod.o shr_const_mod.o shr_log_mod.o $(OBJS_NOMPI)

OBJS_LOG  := test_shr_log.o shr_log_mod.o shr_kind_mod.o \
	     test_mod.o shr_sys_mod.o shr_mpi_mod.o $(OBJS_NOMPI)

WRFESMF_OBJS := ESMF_AlarmMod.o ESMF_BaseTimeMod.o ESMF_ClockMod.o ESMF_TimeMod.o \
                Meat.o ESMF_BaseMod.o ESMF_CalendarMod.o ESMF_FractionMod.o   \
                ESMF_TimeIntervalMod.o ESMF_Stubs.o ESMF_Mod.o \
                ESMF_AlarmClockMod.o wrf_error_fatal.o wrf_message.o

ifeq ($(ESMF_BLD),$(null))
   OBJS_STIN   += $(WRFESMF_OBJS)
   OBJS_STRMS  += $(WRFESMF_OBJS)
   OBJS_SCAM   += $(WRFESMF_OBJS)
endif

#
# Executables:
#

debug:
	@echo "VPATH:      " $(VPATH)
	@echo "ESMF_MOD:   " $(ESMF_MOD)
	@echo "ESMF_ARCH:  " $(ESMF_ARCH)
	@echo "FC:         " $(FC)
	@echo "INC_NETCDF: " $(INC_NETCDF)
	@echo "LIB_MPI:    " $(LIB_MPI)
test_shr_sys: $(OBJS)
	$(LD) -o test_shr_sys     $(OBJS) $(LDFLAGS)
test_shr_file: $(OBJS_FILE)
	$(LD) -o test_shr_file    $(OBJS_FILE) $(LDFLAGS)
test_shr_orb: $(OBJS_ORB)
	$(LD) -o test_shr_orb     $(OBJS_ORB) $(LDFLAGS)
test_shr_streams: $(OBJS_STRMS)
	$(LD) -o test_shr_streams $(OBJS_STRMS) $(LDFLAGS)
test_shr_tInterp: $(OBJS_STIN)
	$(LD) -o test_shr_tInterp $(OBJS_STIN) $(LDFLAGS)
test_shr_mpi: $(OBJS_MPI)
	$(LD) -o test_shr_mpi $(OBJS_MPI) $(LDFLAGS)
test_shr_scam: $(OBJS_SCAM)
	$(LD) -o test_shr_scam $(OBJS_SCAM) $(LDFLAGS)
test_shr_log: $(OBJS_LOG)
	$(LD) -o test_shr_log $(OBJS_LOG) $(LDFLAGS)

clean:
	$(RM) -f *.mod *.o *.f *.f90 F mpif.h test_shr_sys test_shr_orb \
		test_shr_file tests_shr_streams tests_shr_tInterp \
                test_shr_mpi libesmf.a test_shr_scam test_shr_log Depends Srcfiles Filepath

#
# Dependencies
#
Depends: Srcfiles Filepath
	./Mkdepends Filepath Srcfiles > $@

paths    := $(subst $(space),"\n",$(cpp_dirs))

Srcfiles: Filepath
	./Mksrcfiles > $@

Filepath:
	@echo -e $(paths) > $@

-include Depends

# ESMF code...
ifeq ($(ESMF_BLD),$(null))

AR  := ar
CPP := cpp

libesmf.a : $(WRFESMF_OBJS)
	$(RM) -f libesmf.a
	$(AR) $(ARFLAGS) libesmf.a $(WRFESMF_OBJS)
	$(RANLIB) libesmf.a

endif
