
all:
	( $(MAKE) build_registry CPPFLAGS="$(CPPFLAGS)" CC="$(CC)" CFLAGS="$(CFLAGS)" )
	( $(MAKE) build_input_gen CPPFLAGS="$(CPPFLAGS)" CC="$(CC)" CFLAGS="$(CFLAGS)" )

build_input_gen:
	(cd input_gen; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CC="$(CC)" CFLAGS="$(CFLAGS)")

build_registry:
	(cd registry; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CC="$(CC)" CFLAGS="$(CFLAGS)")

clean:
	(cd input_gen; $(MAKE) clean)
	(cd registry; $(MAKE) clean)
