
SHELL = /bin/sh

include Makefile.conf

SUBDIRS = $(MPEUPATH) $(MCTPATH)

# TARGETS
subdirs:
	for dir in $(SUBDIRS); do \
	  cd $$dir;               \
	  $(MAKE);                \
	  cd ..;                  \
	done

clean:
	for dir in $(SUBDIRS); do \
	  cd $$dir;               \
	  $(MAKE) clean;          \
	  cd ..;                  \
	done

install: subdirs
	for dir in $(SUBDIRS); do \
	  cd $$dir;               \
	  $(MAKE) install;        \
	  cd ..;                  \
	done

example: subdirs
	cd $(EXAMPLEPATH) && $(MAKE)


