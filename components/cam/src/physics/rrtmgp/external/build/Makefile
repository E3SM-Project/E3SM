#!/bin/make
RTE_DIR = ../rte
RTE_KERNEL_DIR = ../rte/kernels
RRTMGP_DIR = ../rrtmgp
RRTMGP_KERNEL_DIR = ../rrtmgp/kernels
VPATH = $(RTE_DIR):$(RTE_KERNEL_DIR):$(RRTMGP_DIR):$(RRTMGP_KERNEL_DIR)

include Makefile.conf
include Makefile.rules

all: librrtmgp.a librte.a

include $(RTE_DIR)/Make.depends
include $(RRTMGP_DIR)/Make.depends

librte.a: $(RTE_SRC)
	ar -rvs librte.a $(RTE_SRC)

librrtmgp.a: $(RRTMGP_SRC) 
	ar -rvs librrtmgp.a $(RRTMGP_SRC)

clean:
	rm -f *.optrpt *.mod *.o librrtmgp.a librte.a
