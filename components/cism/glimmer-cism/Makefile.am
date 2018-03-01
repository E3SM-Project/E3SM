SUBDIRS = m4macros scripts libglimmer libglimmer-solve libglimmer-trilinos libglide libglint \
          utils extensions tests example-drivers doc
if BUILD_PYTHON
  SUBDIRS += PyGC
endif

 
EXTRA_DIST = glimmer-cism-config.in
bin_SCRIPTS=glimmer-cism-config

datadir = tests example-drivers

data::
	list='$(datadir)'; for dir in $$list; do \
          (cd $$dir && $(MAKE) $@); \
        done

