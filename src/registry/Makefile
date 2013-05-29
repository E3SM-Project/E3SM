.SUFFIXES: .c .o

OBJS = parse.o dictionary.o gen_inc.o fortprintf.o ezxml/ezxml.o

all: parse

ezxml/ezxml.o:
	(cd ezxml; $(CC) -c ezxml.c)

parse: $(OBJS)
	$(CC) -o $@ $(OBJS)

clean:
	$(RM) *.o ezxml/*.o parse

.c.o:
	$(CC) -c $<
