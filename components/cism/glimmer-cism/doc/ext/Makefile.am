SUBDIRS=
if BUILD_EXT_EROSION
 SUBDIRS += erosion
endif
if BUILD_EXT_PHAML
 SUBDIRS += libphaml
endif


BUILT_SOURCES = selected_ext.tex
EXTRA_DIST = ext.tex ext-app.tex selected_ext.tex.in

editext = sed -e 's|@build_extension_erosion[@]|$(build_extension_erosion)|g' \
              -e 's|@build_extension_libphaml[@]|$(build_extension_libphaml)|g'

selected_ext.tex:	selected_ext.tex.in Makefile
		$(editext) '$(srcdir)/$@.in' >$@
