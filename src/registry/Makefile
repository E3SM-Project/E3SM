.SUFFIXES: .c .o

OBJS = parse.o dictionary.o gen_inc.o fortprintf.o utility.o

all: ezxml parse
	($(MAKE) clean)
	($(MAKE) parse CPPFLAGS="$(CPPFLAGS)" CPPINCLUDES="$(CPPINCLUDES)")

ezxml:
	(cd ../external/ezxml; $(MAKE) )

parse: ezxml $(OBJS)
	$(CC) $(CPPFLAGS) ../external/ezxml/ezxml.o  -I../external/ezxml -o $@ $(OBJS)

clean:
	$(RM) *.o parse
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.c.o: ezxml
	$(CC) $(CPPFLAGS) -I../external/ezxml -c $<
