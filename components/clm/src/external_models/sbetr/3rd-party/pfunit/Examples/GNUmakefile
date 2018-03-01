
ALLDIRS = $(shell find . -maxdepth 1 -mindepth 1 -type d -exec basename {} \;)
MPIDIRS = $(wildcard MPI*)
EXDIRS := $(filter-out $(MPIDIRS),$(ALLDIRS))

MPI ?= NO
ifeq ($(MPI),YES)
  EXDIRS := $(ALLDIRS)
endif

all:
	@for dir in $(EXDIRS); do \
	make -C $$dir all; done

tests:
	@for dir in $(EXDIRS); do \
	make -C $$dir tests; done

clean:
	@for dir in $(ALLDIRS); do \
	make -C $$dir clean; done
	rm -f tests.x
