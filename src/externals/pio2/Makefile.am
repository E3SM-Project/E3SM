# This is part of PIO. It creates the main Makefile.

# Ed Hartnett

# Look in the m4 directory for autotools stuff.
ACLOCAL_AMFLAGS= -I m4

# Does the user want to build fortran?
if BUILD_DOCS
DOC = doc
endif

SUBDIRS = src tests examples ${DOC} scripts

EXTRA_DIST = CMakeLists.txt set_flags.am COPYRIGHT

