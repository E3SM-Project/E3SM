.SUFFIXES: .c .o

EZXML_PATH= ../../external/ezxml

NL_OBJS = namelist_gen.o test_functions.o
ST_OBJS = streams_gen.o test_functions.o
XML_OBJS = $(EZXML_PATH)/ezxml.o

all: ezxml
	($(MAKE) clean)
	($(MAKE) namelist_gen CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)")
	($(MAKE) streams_gen CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)")

ezxml:
	(cd $(EZXML_PATH); $(MAKE) )

namelist_gen: ezxml $(NL_OBJS) $(XML_OBJS)
	$(CC) $(CPPFLAGS) -I$(EZXML_PATH) -o $@ $(NL_OBJS) $(XML_OBJS)

streams_gen: ezxml $(ST_OBJS) $(XML_OBJS)
	$(CC) $(CPPFLAGS) -I$(EZXML_PATH) -o $@ $(ST_OBJS) $(XML_OBJS)

clean:
	$(RM) *.o namelist_gen streams_gen
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.c.o:
	$(CC) $(CPPFLAGS) -I$(EZXML_PATH) -c $<
