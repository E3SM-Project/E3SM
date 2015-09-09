SUBDIRS = eis EISMINT-3 glint simple_glide
tests = eis EISMINT-3

data::
		list='$(tests)'; for dir in $$list; do \
		  (cd $$dir && $(MAKE) $@); \
		done

