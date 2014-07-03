SUBDIRS = EISMINT exact-isothermal higher-order

data::
		list='$(SUBDIRS)'; for dir in $$list; do \
		  (cd $$dir && $(MAKE) $@); \
		done
