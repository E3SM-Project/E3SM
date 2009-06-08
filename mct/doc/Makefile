#!/bin/make
#-----------------------------------------------------------------------
# Documentation
all: 
	cd texsrc; make
	make apis
	
apis:
	cd texsrc; make
	make apisdvi

apisdvi: mct_APIs.dvi

clean:
	cd texsrc; make clean
	rm -f *.dvi *.log *.bbl *.blg *.aux *.toc

.SUFFIXES: .dvi .tex

.tex.dvi:
	latex $*.tex

#.
