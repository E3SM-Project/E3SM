
SHELL = /bin/sh

MPEUPATH = ./mpeu
MCTPATH = ./mct
EXAMPLEPATH = ./ut-mct

SUBDIRS = $(MPEUPATH) $(MCTPATH)

# TARGETS
.PHONY: subdirs $(SUBDIRS) clean

subdirs: $(SUBDIRS)

# RULES
$(SUBDIRS):
	$(MAKE) -C $@

ut-mct: $(SUBDIRS)
	$(MAKE) -C $(EXAMPLEPATH)

clean: 
	@for dir in $(SUBDIRS); do \
	  $(MAKE) -C $$dir clean; \
	done

install: $(SUBDIRS)
	@for dir in $(SUBDIRS); do \
	  $(MAKE) -C $$dir install; \
	done

# DEPENDENCIES
$(MCTPATH): $(MPEUPATH)
