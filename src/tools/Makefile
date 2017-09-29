
all:
	( $(MAKE) build_registry CPPFLAGS="$(CPPFLAGS)" CC="$(CC)" CFLAGS="$(CFLAGS)" )
	( $(MAKE) build_input_gen CPPFLAGS="$(CPPFLAGS)" CC="$(CC)" CFLAGS="$(CFLAGS)" )

build_input_gen:
ifdef TOOL_DIR
	(cp $(TOOL_DIR)/namelist_gen input_gen/)
	(cp $(TOOL_DIR)/streams_gen input_gen/)
else
	(cd input_gen; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CC="$(CC)" CFLAGS="$(CFLAGS)")
endif

build_registry:
ifdef TOOL_DIR
	(cp $(TOOL_DIR)/parse registry/)
else
	(cd registry; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CC="$(CC)" CFLAGS="$(CFLAGS)")
endif

clean:
	(cd input_gen; $(MAKE) clean)
	(cd registry; $(MAKE) clean)
