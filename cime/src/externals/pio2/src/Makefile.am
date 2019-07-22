# This is part of PIO. It creates the Makefile for the src directory.

# Ed Hartnett

# Does the user want to build fortran?
if BUILD_FORTRAN
FLIB = flib
endif

if USE_GPTL
GPTL = gptl
endif

SUBDIRS = clib ${GPTL} $(FLIB)

EXTRA_DIST = CMakeLists.txt
