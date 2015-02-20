.SUFFIXES: .c .o

EZXML_PATH= ../../external/ezxml

OBJS = parse.o dictionary.o gen_inc.o fortprintf.o utility.o

all: ezxml parse
	($(MAKE) clean)
	($(MAKE) parse CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)")

ezxml:
	(cd $(EZXML_PATH); $(MAKE) )

parse: ezxml $(OBJS)
	$(CC) $(CPPFLAGS) $(EZXML_PATH)/ezxml.o  -I$(EZXML_PATH) -o $@ $(OBJS)

clean:
	$(RM) *.o parse
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.c.o: ezxml
	$(CC) $(CPPFLAGS) -I$(EZXML_PATH) -c $<
