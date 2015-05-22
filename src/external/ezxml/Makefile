.SUFFIXES: .c .o

OBJS = ezxml.o

all: clean
	$(MAKE) library

library: $(OBJS)

clean:
	$(RM) *.o *.i

.c.o: $(CC) $(CFLAGS) $(CPPFLAGS) -c $<
