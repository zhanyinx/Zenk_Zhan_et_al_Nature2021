CC= cc
LFLAGS=  -Wall -lgsl -lgslcblas -lm -L/opt/local/lib
CFLAGS=  -Wall -O3 -funroll-all-loops -finline-functions -I/opt/local/include 

.o:     $(CC)  $(CFLAGS) $<


mhistogram: mhistogram.o memory.o misc.o do_mhistogram.o io.o
	$(CC)  mhistogram.o memory.o misc.o do_mhistogram.o io.o -o ../../bin/mhistogram $(LFLAGS)



