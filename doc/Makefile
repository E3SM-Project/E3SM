#!/bin/make
#-----------------------------------------------------------------------
# Documentation
all: 
	cd texsrc; make
	make users
	make apis
	
users:
	make usersdvi

apis:
	cd texsrc; make
	make apisdvi

usersdvi: MCTUsersGuide.dvi

apisdvi: mct_APIs.dvi

clean:
	cd texsrc; make clean
	rm -f *.dvi *.log *.bbl *.blg *.aux *.toc

.SUFFIXES: .dvi .tex

.tex.dvi:
	latex $*.tex

#.
