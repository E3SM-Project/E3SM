
all: build_input_gen build_registry

build_input_gen:
	(cd input_gen; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CC="$(CC)")

build_registry:
	(cd registry; $(MAKE) CPPFLAGS="$(CPPFLAGS)" CC="$(CC)")

clean:
	(cd input_gen; $(MAKE) clean)
	(cd registry; $(MAKE) clean)
