#
# Makefile to create HTML documentation of namelists
#
SCRNLDIR := ../../../../../../scripts/doc/modelnl
VPATH    := . $(SCRNLDIR) ../../../../../drv/bld/namelist_files ../../../../../glc/cism/bld/namelist_files \
	    ../../../../../atm/datm/bld/namelist_files ../../../../../../scripts/ccsm_utils/Case.template \
	    ../../../../../../scripts/ccsm_utils/Machines ../../../../../rof/rtm/bld/namelist_files \
	    ../../../bld/namelist_files ../..
SOURCES  := namelist_definition_drv.xml namelist_definition_cism.xml namelist_definition.xml \
	  namelist_definition_rtm.xm namelist_definition_datm.xml config_definition.xml \
	  config_grid.xml config_machines.xml config_compsets.xml ChangeSum
TAGFILE  := clmtag.txt

CWD      := $(shell pwd )
ALLOUT   := $(CWD)/clm_nl_drv.html  $(CWD)/clm_nl_cism.html  $(CWD)/clm_nl_clm.html    $(CWD)/clm_nl_rtm.html \
	    $(CWD)/clm_nl_datm.html $(CWD)/clm_env_case.html $(CWD)/clm_env_build.html $(CWD)/clm_env_pesetup.html \
	    $(CWD)/clm_env_run.html $(CWD)/clm_grid.html     $(CWD)/clm_machines.html  $(CWD)/clm_compsets.html \
	    $(TAGFILE)              $(CWD)/index.html

all: $(ALLOUT)

debug:
	@echo "SOURCES  = $(SOURCES)"
	@echo "VPATH    = $(VPATH)"
	@echo "ALLOUT   = $(ALLOUT)"
	@echo "SCRNLDIR = $(SCRNLDIR)"
	@echo "CWD      = $(CWD)"
	@echo "TAGFILE  = $(TAGFILE)"

.SUFFIXES:
.SUFFIXES: .xml .html .txt

RM := /bin/rm

CTAGNAME = $(shell cat $(TAGFILE) )


$(TAGFILE): ChangeSum
	head -3 $< | tail -1 | awk '{print $$1}' > $@

$(CWD)/index.html: $(TAGFILE) index.cpp
	sed 's/CLMTAGNAME/$(CTAGNAME)/' index.cpp > $@

$(CWD)/clm_nl_drv.html: namelist_definition_drv.xml
	cd $(SCRNLDIR) ; ./nldef2html_drv > $@

$(CWD)/clm_nl_cism.html: namelist_definition_cism.xml
	cd $(SCRNLDIR) ; ./nldef2html_cism > $@

$(CWD)/clm_nl_clm.html: namelist_definition.xml
	cd $(SCRNLDIR) ; ./nldef2html_clm > $@

$(CWD)/clm_nl_rtm.html: namelist_definition_rtm.xml
	cd $(SCRNLDIR) ; ./nldef2html_rtm > $@

$(CWD)/clm_nl_datm.html: namelist_definition_datm.xml
	cd $(SCRNLDIR) ; ./nldef2html_datm > $@

$(CWD)/clm_env_case.html: config_definition.xml
	cd $(SCRNLDIR) ; ./xmldef2html_env_case > $@

$(CWD)/clm_env_build.html: config_definition.xml
	cd $(SCRNLDIR) ; ./xmldef2html_env_build > $@

$(CWD)/clm_env_pesetup.html: config_definition.xml
	cd $(SCRNLDIR) ; ./xmldef2html_env_pesetup > $@

$(CWD)/clm_env_run.html: config_definition.xml
	cd $(SCRNLDIR) ; ./xmldef2html_env_run > $@

$(CWD)/clm_grid.html: config_grid.xml
	cd $(SCRNLDIR) ; ./xmldef2html_grid > $@

$(CWD)/clm_machines.html: config_machines.xml
	cd $(SCRNLDIR) ; ./xmldef2html_machines > $@

$(CWD)/clm_compsets.html: config_compsets.xml
	cd $(SCRNLDIR) ; ./xmldef2html_compsets > $@

clean:
	$(RM) -f $(ALLOUT)

realclean: clean
	$(RM) -f $(TAGFILE)

