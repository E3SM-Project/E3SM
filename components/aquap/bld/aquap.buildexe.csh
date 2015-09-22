#! /usr/bin/env perl
use strict;

chdir "obj";

open(F,">Filepath") or die "Could not open file Filepath to write";
print F "$ENV{CASEROOT}/SourceMods/src.aquap\n";
print F "$ENV{CODEROOT}/atm/cam/src/utils/cam_aqua\n";
print F "$ENV{CODEROOT}/atm/cam/src/utils/cam_aqua/cpl\n";
close(F);

system("gmake complib -j $ENV{GMAKE_J} MODEL=aquap COMPLIB=$ENV{LIBROOT}/libocn.a -f $ENV{CASETOOLS}/Makefile USER_INCLDIR=-I$ENV{EXEROOT}/atm/obj");
