#!/bin/csh

if (-e ./ug.pdf) then
  if (-e ./ug.pdf.prev) rm -f ./ug.pdf.prev
  mv ./ug.pdf ./ug.pdf.prev
endif

cp -f index.html index.html.backup
rm *.html
cp -f index.html.backup index.html
docbook2html --dsl stylesheet.dsl ug.xml
docbook2pdf ug.xml
chgrp cgdcsmweb *
chmod -R g+w .

