#!/bin/make
#-----------------------------------------------------------------------
# Documentation
all: 
	cd texsrc; make
	make design
	make apis
	
design:
	cd texsrc; make
	make designdvi

apis:
	cd texsrc; make
	make apisdvi

designdvi: mct_design.dvi

apisdvi: mct_APIs.dvi

clean:
	cd texsrc; make clean
	rm -f *.dvi *.log *.bbl *.blg *.aux *.toc

.SUFFIXES: .dvi .tex

.tex.dvi:
	latex $*.tex

#.
