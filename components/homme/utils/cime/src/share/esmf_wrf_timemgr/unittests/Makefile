
cpp_dirs := . ..
cpp_path := $(foreach dir,$(cpp_dirs),-I$(dir)) # format for command line
# Expand any tildes in directory names. Change spaces to colons.
VPATH    := $(foreach dir,$(cpp_dirs),$(wildcard $(dir)))
VPATH    := $(subst $(space),:,$(VPATH))

#VPATH := .:..


.SUFFIXES: .F90 .o .F .f90

AR := ar
FC := xlf95
FFLAGS := -g -qfullpath -qmaxmem=-1 -O2 -qstrict -qsigtrap=xl__trcedump -Q -qinitauto=7FF7FFFF -qflttrap=ov:zero:inv:en -qspillsize=4000 -qarch=auto -qtune=auto -qsuffix=f=f90:cpp=F90 -I. -I.. -WF,-DHIDE_MPI
LDFLAGS :=

OBJS := ESMF_AlarmMod.o ESMF_BaseTimeMod.o ESMF_ClockMod.o ESMF_TimeMod.o \
        MeatMod.o ESMF_BaseMod.o ESMF_CalendarMod.o ESMF_FractionMod.o   \
        ESMF_TimeIntervalMod.o ESMF_Stubs.o ESMF.o ESMF_ShrTimeMod.o \
	ESMF_AlarmClockMod.o wrf_stuff.o

test: libesmf_time.a test.o
	$(FC) $(LDFLAGS) -o test test.o -L. -lesmf_time

lib: libesmf_time.a

debug: $(OBJS)
      echo "VPATH : $VPATH"
      echo "OBJS  : $OBJ"
      echo "FFLAGS: $FFLAGS"

libesmf_time.a : $(OBJS)
	\rm -f libesmf_time.a
	$(AR) $(ARFLAGS) libesmf_time.a $(OBJS)

.F90.o :
	$(FC) -c $(FFLAGS) $<

clean:
	/bin/rm -f *.o libesmf_time.a *.mod test

# DEPENDENCIES : only dependencies after this line

#$$$ update dependencies!

ESMF_BaseMod.o : ESMF_BaseMod.F90 wrf_stuff.o
ESMF_FractionMod.o: ESMF_FractionMod.F90
MeatMod.o : MeatMod.F90 ESMF_BaseMod.o
ESMF_BaseTimeMod.o : ESMF_BaseTimeMod.F90 ESMF_BaseMod.o
ESMF_CalendarMod.o : ESMF_CalendarMod.F90 ESMF_BaseMod.o ESMF_BaseTimeMod.o
ESMF_Stubs.o : ESMF_Stubs.F90 ESMF_CalendarMod.o ESMF_BaseMod.o
ESMF_ShrTimeMod.o :  ESMF_ShrTimeMod.F90 ESMF_CalendarMod.o ESMF_BaseTimeMod.o ESMF_BaseMod.o
ESMF_TimeIntervalMod.o : ESMF_TimeIntervalMod.F90 ESMF_FractionMod.o
ESMF_TimeMod.o : ESMF_TimeMod.F90 ESMF_ShrTimeMod.o ESMF_Stubs.o ESMF_TimeIntervalMod.o
ESMF_AlarmMod.o : ESMF_AlarmMod.F90 ESMF_BaseTimeMod.o ESMF_TimeMod.o ESMF_TimeIntervalMod.o
ESMF_ClockMod.o : ESMF_ClockMod.F90 ESMF_BaseTimeMod.o ESMF_TimeMod.o ESMF_TimeIntervalMod.o ESMF_AlarmMod.o
ESMF_AlarmClockMod.o : ESMF_AlarmClockMod.F90 ESMF_AlarmMod.o ESMF_ClockMod.o ESMF_TimeIntervalMod.o ESMF_TimeMod.o
ESMF.o : ESMF.F90 ESMF_AlarmMod.o ESMF_BaseTimeMod.o ESMF_ClockMod.o ESMF_TimeMod.o \
        ESMF_BaseMod.o ESMF_CalendarMod.o ESMF_FractionMod.o  ESMF_ShrTimeMod.o  \
        ESMF_TimeIntervalMod.o ESMF_Stubs.o ESMF_AlarmClockMod.o
test.o : test.F90 ESMF.o

