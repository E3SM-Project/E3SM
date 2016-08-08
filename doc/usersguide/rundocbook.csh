#!/bin/csh

if (-e ./ug.pdf) then
  if (-e ./ug.pdf.prev) rm -f ./ug.pdf.prev
  mv ./ug.pdf ./ug.pdf.prev
endif

rm *.html *.htm
docbook2html --dsl stylesheet.dsl ug.xml
docbook2pdf ug.xml

foreach file (*)
    if (-o $file) then
	chgrp cgdcesmweb $file
	chmod g+w $file
    endif
end
