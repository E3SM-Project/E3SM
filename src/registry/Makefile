.SUFFIXES: .c .o

OBJS = parse.o dictionary.o gen_inc.o fortprintf.o ezxml/ezxml.o

all: parse

ezxml/ezxml.o:
	(cd ezxml; $(CC) -c ezxml.c)

parse: $(OBJS)
	$(CC) -o $@ $(OBJS)

clean:
	$(RM) *.o ezxml/*.o parse
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.c.o:
	$(CC) -c $<
