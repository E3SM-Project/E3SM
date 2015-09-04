#!/bin/csh

if (-e ./ug.pdf) then
  if (-e ./ug.pdf.prev) rm -f ./ug.pdf.prev
  mv ./ug.pdf ./ug.pdf.prev
endif

mv index.html index.html.hold
rm *.html
mv index.html.hold index.html
docbook2html ug.xml
docbook2pdf ug.xml
chgrp cgdcsmweb *
chmod g+w *
