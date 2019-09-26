#CC=gcc -fopenmp
CC=cc
CCMP=/usr/lib64/openmpi/bin/mpicc
#LFLAGS -pg
#LFLAGS= -Wall -pg -fopenmp -lm -L/usr/local/lib
LFLAGS= -Wall -pg  -lm -L/usr/local/lib

LGSLFLAGS= -lgsl -lgslcblas 
CFLAGS=  -Wall   -I/opt/local/include -Iinclude

SRCDIR=src
OBJDIR=obj
INCL = -I../include
INCMPIFLAG=-D ACTIVE_MPI
INCSMPFLAG=-D STEMPERING

BINDIR=bin


ifeq ($(version),MPI)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CCMP) $(INCMPIFLAG) $(CFLAGS) -c $< $(INCL) -o $@


montegrappa_MPI: $(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/mc.o $(OBJDIR)/local_move.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/optimizepot.o $(OBJDIR)/montegrappa.o $(OBJDIR)/MPIfunc.o
	$(CCMP) $(INCMPIFLAG) $(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/mc.o $(OBJDIR)/local_move.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/optimizepot.o $(OBJDIR)/montegrappa.o $(OBJDIR)/MPIfunc.o -o $(BINDIR)/montegrappa_mpi $(LFLAGS) 



else ifeq ($(version),STEMPERING) 

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(INCSMPFLAG) $(CFLAGS) -c $< $(INCL) -o $@


montegrappa_stemp:	$(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/mc.o $(OBJDIR)/local_move.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/optimizepot.o $(OBJDIR)/montegrappa.o $(OBJDIR)/stempering.o $(OBJDIR)/memory1.o $(OBJDIR)/adjust_st.o $(OBJDIR)/do_mhistogram.o
	$(CC) $(INCSMPFLAG) $(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/mc.o  $(OBJDIR)/local_move.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/optimizepot.o $(OBJDIR)/montegrappa.o $(OBJDIR)/stempering.o $(OBJDIR)/memory1.o $(OBJDIR)/adjust_st.o $(OBJDIR)/do_mhistogram.o  -o $(BINDIR)/montegrappa $(LFLAGS) $(LGSLFLAGS)

else

$(OBJDIR)/%.o: $(SRCDIR)/%.c 
	$(CC) $(CFLAGS) -c $< $(INCL) -o $@

montegrappa:  $(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/mc.o $(OBJDIR)/local_move.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/optimizepot.o $(OBJDIR)/montegrappa.o 
	$(CC) $(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/mc.o  $(OBJDIR)/local_move.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/optimizepot.o $(OBJDIR)/montegrappa.o -o $(BINDIR)/montegrappa $(LFLAGS)

endif



grappino:    $(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/energy.o $(OBJDIR)/pdb.o $(OBJDIR)/rotamers.o $(OBJDIR)/grappino.o $(OBJDIR)/optimizepot.o	
	$(CC) $(OBJDIR)/geometry.o $(OBJDIR)/io.o $(OBJDIR)/memory.o $(OBJDIR)/misc.o $(OBJDIR)/potential.o $(OBJDIR)/energy.o $(OBJDIR)/pdb.o $(OBJDIR)/rotamers.o $(OBJDIR)/optimizepot.o $(OBJDIR)/grappino.o -o $(BINDIR)/grappino $(LFLAGS)



	
mhistogram:
	cd src/mhistogram; make


clean:
	rm -f $(OBJDIR)/*.o src/mhistogram/*.o  $(BINDIR)/montegrappa*

cleanobj:
	rm -f $(OBJDIR)/*.o 
	rm -f src/mhistogram/*.o



all:
	make cleanobj;
	make version=STEMPERING;
	make cleanobj;
	make grappino;
	make mhistogram;
	make cleanobj;
	make version=MPI;
