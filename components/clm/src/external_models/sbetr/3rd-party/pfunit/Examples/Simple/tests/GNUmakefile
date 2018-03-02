.PHONY: clean tests
ifneq ($(BASEMK_INCLUDED),YES)
include $(PFUNIT)/include/base.mk
endif

SRCS = $(wildcard *.pf)
OBJS = $(SRCS:.pf=$(OBJ_EXT))

# Support for Windows/CYGWIN
ifeq ($(UNAME),Windows)
	TOP_DIR := ..
	# Since this is a Windows compilation, add the Windows define.
	FFLAGS += -DWindows
	ifeq ($(FC),nagfor)
		FFLAGS += -fpp
	endif
	ifeq ($(FC),ifort)
		FFLAGS += /nologo
	endif
endif

ifndef SKIP_INTENTIONALLY_BROKEN
  FFLAGS += -DINCLUDE_INTENTIONALLY_BROKEN
endif

tests: $(OBJS)

testSuites.inc: $(SRCS)

FFLAGS += -I$(SRC_DIR) -I$(PFUNIT)/mod

ifneq ($(UNAME),Windows)
	ifeq ($(FC),nagfor)
	   FFLAGS += -fpp
	endif
endif

%.F90: %.pf
	$(PFUNIT)/bin/pFUnitParser.py $<  $@

%$(OBJ_EXT): %.F90
	$(FC) -c $(FFLAGS) $(FPPFLAGS) $<

clean: local-E1-clean

local-E1-clean:
	$(RM) *$(OBJ_EXT) *.mod
