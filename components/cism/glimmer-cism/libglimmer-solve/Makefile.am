if EXT_BLAS
 BLAS_DIR =
 BLAS_LIB = 
else
 BLAS_DIR = blas
 BLAS_LIB = blas/libblas.la
endif

SUBDIRS = SLAP $(BLAS_DIR)

lib_LTLIBRARIES = libglimmer-solve.la
CLEANFILES = $(BUILT_SOURCES) f90_dep.mak $(wildcard *.mod)

libglimmer_solve_la_SOURCES = glimmer_sparse_type.F90 glimmer_sparse_slap.F90 \
			      glimmer_sparse.F90 cism_sparse_pcg.F90
libglimmer_solve_la_LIBADD = SLAP/libslap.la $(BLAS_LIB)
libglimmer_solve_la_LDFLAGS=-all-static
nodist_include_HEADERS := $(wildcard *.mod)

AM_FCFLAGS = $(FC_MODINC)../libglimmer

-include f90_dep.mak
include $(top_srcdir)/extra_rules.am
