
SHELL = /bin/sh

MPEUPATH = ./mpeu
MCTPATH = ./mct
EXAMPLEPATH = ./ut-mct

SUBDIRS = $(MPEUPATH) $(MCTPATH)

# TARGETS
subdirs: $(SUBDIRS)
	@ argv="$(SUBDIRS)" ; \
	for dir in $$argv; do \
	  cd $$dir;           \
	  $(MAKE);            \
	  cd ..;              \
	done

clean: $(SUBDIRS) 
	@ argv="$(SUBDIRS)" ; \
	for dir in $$argv; do \
	  cd $$dir;           \
	  $(MAKE) clean;      \
	  cd ..;              \
	done

install: subdirs
	@ argv="$(SUBDIRS)" ; \
	for dir in $$argv; do \
	  cd $$dir;           \
	  $(MAKE) install;    \
	  cd ..;              \
	done

ut-mct: subdirs
	cd $(EXAMPLEPATH) && $(MAKE)


