
SHELL = /bin/sh

include Makefile.conf

SUBDIRS = $(MPEUPATH) $(MCTPATH)

# TARGETS
subdirs:
	@for dir in $(SUBDIRS); do \
	  cd $$dir;                \
	  $(MAKE);                 \
	  cd $(abs_top_builddir);  \
	done

clean:
	@for dir in $(SUBDIRS); do \
	  cd $$dir;                \
	  $(MAKE) clean;           \
	  cd $(abs_top_builddir);  \
	done

install: subdirs
	@for dir in $(SUBDIRS); do \
	  cd $$dir;                \
	  $(MAKE) install;         \
	  cd $(abs_top_builddir);  \
	done

example: subdirs
	@cd $(EXAMPLEPATH) && $(MAKE)


