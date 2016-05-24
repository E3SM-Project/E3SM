#INC = -I/BlueLight/ppcfloor/bglsys/include
#INC = -I/bgl/BlueLight/ppcfloor/bglsys/include
INC = 

ifneq ('$(OBJECT_MODE)','32')
  export OBJECT_MODE=64
endif


default : timers.c
	$(CC) $(CFLAGS) -c timers.c -o timers.o
	rm -f libtimers.a
	ar -vq libtimers.a timers.o

blrts : timers.c
	blrts_xlc -c -O -qarch=440 -qmaxmem=-1 -qcpluscmt -qinfo=pro:gen $(INC) timers.c -o timers.o
	rm -f libtimers.a
	ar -vq libtimers.a timers.o 

ibm : timers.c timebase$(OBJECT_MODE).o
	xlc -c -O -qarch=auto -qmaxmem=-1 -qcpluscmt -qinfo=pro:gen $(INC) timers.c -o timers.o
	rm -f libtimers.a
	ar -vq libtimers.a timers.o timebase$(OBJECT_MODE).o

clean :
	rm -f timers.o libtimers.a timebase*.o
