#
# Makefile to convert DocBook CLM Users-Guide into html and/or pdf
# (rtf, txt, ps, tex, man, dvi, and texi are also valid docbook formats)
#
VPATH := ../../tools/cprnc . .. ../../bld ../../tools/ncl_scripts \
	../../tools/mksurfdata ../../test/system ../../bld/namelist_files \
	../../bld/config_files ../../tools ../../../../../scripts/ccsm_utils/Tools \
        ../../../../../scripts/ccsm_utils/Tools/lnd/clm/PTCLM/ ../../src/main

PDFUG     := clm_ug.pdf
HTMLUG    := book1.html
DOCBKUG   := clm_ug.xml
CFGLOG    := config_help
CPRLOG    := cprnc_readme
BNMLOG    := buildnml_help
BSTLOG    := build_streams_help
RESLOG    := buildnml_resolutions
USCLOG    := buildnml_usecases
QCKLOG    := quickstart_guide
COPLOG    := filecopies
MKSLOG    := mksurfdata.pl
USRLOG    := quickstart_usrdat
PTCLOG    := ptclm_help
PTCLST    := ptclm_list
TDRLOG    := test_driver.sh
GETREG    := getregional_datasets
CO2DIF    := addco2_datm.buildnml
DATLOG    := build_date
NMLDFTBL  := namelist_definition_table
NMLDLTBL  := namelist_defaults_clm_table
HSFLDTBL  := history_fields_table
CFGDFTBL  := config_definition_table
COMPLIST  := compsets_list_ofIcases.xml
SOURCES   := $(DOCBKUG) $(COMPLIST) $(CFGLOG).xml $(PTCLOG).xml $(BNMLOG).xml \
             $(RESLOG).xml $(USCLOG).xml $(QCKLOG).xml $(COPLOG).xml $(PTCLST).xml \
             $(USRLOG).xml $(GETREG).xml preface.xml custom.xml special_cases.xml \
             tools.xml adding_files.xml single_point.xml addco2_datm.buildnml.xml \
             appendix.xml trouble_shooting.xml ptclm.xml $(BSTLOG).xml \
             $(MKSLOG).xml $(TDRLOG).xml $(DATLOG).xml $(CPRLOG).xml \
	     $(NMLDFTBL).xml $(NMLDLTBL).xml $(CFGDFTBL).xml $(HSFLDTBL).xml

CONVAMP    := sed 's/\&/\&amp;/g'
CONVSIGNS  := sed 's/>/\&gt;/g' | sed 's/</\&lt;/g'

all: $(HTMLUG) $(PDFUG)

.SUFFIXES:
.SUFFIXES: .xml .pdf .html .tlog .diff .xhtml .xsl

RM := /bin/rm

$(CFGLOG).xml: $(CFGLOG).tlog
$(PTCLOG).xml: $(PTCLOG).tlog
$(PTCLST).xml: $(PTCLST).tlog
$(BNMLOG).xml: $(BNMLOG).tlog
$(BSTLOG).xml: $(BSTLOG).tlog
$(RESLOG).xml: $(RESLOG).tlog
$(USCLOG).xml: $(USCLOG).tlog
$(QCKLOG).xml: $(QCKLOG).tlog
$(COPLOG).xml: $(COPLOG).tlog
$(USRLOG).xml: $(USRLOG).tlog
$(GETREG).xml: $(GETREG).tlog
$(MKSLOG).xml: $(MKSLOG).tlog
$(CPRLOG).xml: $(CPRLOG).tlog
$(CO2DIF).xml: $(CO2DIF).diff
$(DATLOG).xml: $(DATLOG).tlog

$(NMLDFTBL).xhtml: namelist_definition.xml namelist_definition.xsl
	xsltproc ../../bld/namelist_files/namelist_definition.xsl $< > tmpFile.txt
	addxhtmlhead.pl tmpFile.txt > $@
	$(RM) tmpFile.txt

$(NMLDLTBL).xhtml: namelist_defaults_clm.xml namelist_defaults.xsl
	xsltproc ../../bld/namelist_files/namelist_defaults.xsl $< > tmpFile.txt
	addxhtmlhead.pl tmpFile.txt > $@
	$(RM) tmpFile.txt

$(CFGDFTBL).xhtml: config_definition.xml config_definition.xsl
	xsltproc ../../bld/config_files/config_definition.xsl $< > tmpFile.txt
	addxhtmlhead.pl tmpFile.txt > $@
	$(RM) tmpFile.txt

history_fields.xml: findHistFields.pl
	cd ../../src/main; ./findHistFields.pl >& tmpFile.txt
	$(RM) tmpFile.txt

$(HSFLDTBL).xhtml: history_fields.xml history_fields.xsl
	xsltproc ../../bld/namelist_files/history_fields.xsl $< > tmpFile.txt
	addxhtmlhead.pl tmpFile.txt > $@
	$(RM) tmpFile.txt

.xhtml.xml:
	xsltproc stylesheethtml2docbook.xsl $< > $@

.tlog.xml:
	$(CONVAMP) $< | $(CONVSIGNS) > tempFile.txt
	limitLineLen.pl tempFile.txt  > $@
	$(RM) tempFile.txt

.diff.xml:
	$(CONVAMP) $< | $(CONVSIGNS) > $@

debug:
	@echo "SOURCES:   $(SOURCES)"
	@echo "CONVAMP:   $(CONVAMP)"
	@echo "CONVSIGNS: $(CONVSIGNS)"

$(COMPLIST): 
	./get_Icaselist.pl > $@

$(HTMLUG): $(SOURCES)
	docbook2html --dsl clm_stylesheet.dsl#html $<

$(PDFUG): $(SOURCES)
	docbook2pdf --dsl clm_stylesheet.dsl#print $<

$(BNMLOG).tlog: build-namelist
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../bld/build-namelist -help >& $@

$(BSTLOG).tlog: build_streams
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../../../../scripts/ccsm_utils/Tools/build_streams --help >& $@

$(DATLOG).tlog:
	@echo "Get current build date"
	date  +%b-%d-%Y >& $@

$(RESLOG).tlog: build-namelist
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../bld/build-namelist -res list >& $@

$(USCLOG).tlog: build-namelist
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../bld/build-namelist -use_case list >& $@

$(CFGLOG).tlog: configure
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../bld/configure -help >& $@

$(PTCLOG).tlog: PTCLM.py
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../../../../scripts/ccsm_utils/Tools/lnd/clm/PTCLM/PTCLM.py --help >& $@

$(PTCLST).tlog: PTCLM.py
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	cd ../../../../../scripts/ccsm_utils/Tools/lnd/clm/PTCLM; \
        PTCLM.py --list >& $(CURDIR)/$@

$(MKSLOG).tlog: mksurfdata.pl
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../tools/mksurfdata/mksurfdata.pl -help >& $@

$(TDRLOG).tlog: test_driver.sh
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../test/system/test_driver.sh -help >& $@

$(QCKLOG).tlog: Quickstart.GUIDE
	cp $<  $@

$(COPLOG).tlog: README.filecopies
	cp $<  $@

$(CPRLOG).tlog: README
	cp $<  $@

$(USRLOG).tlog: Quickstart.userdatasets
	cp $<  $@

$(GETREG).tlog:  getregional_datasets.pl
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $(GETREG) is good and redo your make"
	../../tools/ncl_scripts/getregional_datasets.pl -help >& $@

clean:
	$(RM) -f $(HTMLUG) $(PDFUG) *.tlog $(DATLOG).xml *.xhtml *.tex

realclean: clean
	$(RM) -f f*.html c*.html x*.html a*.html i*.html $(COMPLIST) $(CFGLOG).xml \
        $(BNMLOG).xml $(BSTLOG).xml $(PTCLOG).xml  $(PTCLST).xml \
	$(RESLOG).xml $(USCLOG).xml $(USRLOG).xml $(GETREG).xml $(QCKLOG).xml \
	$(CO2DIF).xml *.tlog $(MKSLOG).xml $(TDRLOG).xml $(DATLOG).xml \
	$(NMLDFTBL).xml $(NMLDLTBL).xml $(CFGDFTBL).xml $(CPRLOG).xml \
        $(COPLOG).xml
