.SUFFIXES: .c .o

OBJS = ezxml.o

all: $(OBJS)

clean:
	$(RM) *.o *.i

.c.o: $(CC) $(CPPFLAGS) -c $<
