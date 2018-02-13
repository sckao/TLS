CC=g++ -g -fPIC -Wall -O -ansi -D_GNU_SOURCE -O3 -m64 -I/usr/local/include/
ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags`
LDFLAGS+=-lutil -lboost_iostreams -lboost_system -lboost_filesystem -lgsl -lgslcblas -lm
LIB=`${ROOTSYS}/bin/root-config --libs`
GLIB=`${ROOTSYS}/bin/root-config --glibs`

OBJECTS=AnaInput.o MathTools.o Reader.o Study.o Match.o

all:test.exe

AnaInput.o : AnaInput.cc AnaInput.h
	$(CC) -c -o $@ $< $(ROOTFLAG) $(LIB)

MathTools.o : MathTools.cc MathTools.h
	$(CC) -c -o $@ $< $(ROOTFLAG) $(LIB) 

Reader.o : Reader.cc Reader.h
	$(CC) -c -o $@ $< $(ROOTFLAG) $(LIB)

Study.o : Study.cc Study.h
	$(CC) -c -o $@ $< $(ROOTFLAG) $(LIB)

Match.o : Match.cc Match.h
	$(CC) -c -o $@ $< $(ROOTFLAG) $(LIB)

test.exe : main.cc $(OBJECTS)
	$(CC) -o $@ $< $(OBJECTS) $(ROOTFLAG) $(LIB) $(GLIB) $(LDFLAGS)

clean : 
	rm -rf *.o test.exe
