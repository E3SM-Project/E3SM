#!/usr/bin/env python
# script written by Yuri Alexeev to merge .mod/.dat/.run files in one file with XML headers

import sys

fmod = open(sys.argv[1],"r")
fdat = open(sys.argv[2],"r") 
frun = open(sys.argv[3],"r")

print '<document>\n\
<category>minco</category>\n\
<solver>MINLP</solver>\n\
<inputMethod>AMPL</inputMethod>\n\
\n\
<model><![CDATA['

buffer = 1
while buffer:
  buffer =  fmod.read()
  print buffer

print ']]></model>\n\
\n\
<data><![CDATA['

buffer = 1
while buffer:
  buffer =  fdat.read()
  print buffer

print ']]></data>\n\
\n\
<commands><![CDATA['

buffer = 1
while buffer:
  buffer =  frun.read()
  print buffer

print ']]></commands>\n\
\n\
<comments><![CDATA[]]></comments>\n\
\n\
</document>'

fmod.close()
fdat.close()
frun.close()
