.SUFFIXES: .c .o

NL_OBJS = namelist_gen.o test_functions.o
ST_OBJS = streams_gen.o test_functions.o
XML_OBJS = ../external/ezxml/ezxml.o

all: ezxml
	($(MAKE) clean)
	($(MAKE) namelist_gen CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)")
	($(MAKE) streams_gen CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)")

ezxml:
	(cd ../external/ezxml; $(MAKE) )

namelist_gen: ezxml $(NL_OBJS) $(XML_OBJS)
	$(CC) $(CPPFLAGS) -I../external/ezxml -o $@ $(NL_OBJS) $(XML_OBJS)

streams_gen: ezxml $(ST_OBJS) $(XML_OBJS)
	$(CC) $(CPPFLAGS) -I../external/ezxml -o $@ $(ST_OBJS) $(XML_OBJS)

clean:
	$(RM) *.o namelist_gen streams_gen
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.c.o:
	$(CC) $(CPPFLAGS) -I../external/ezxml -c $<
