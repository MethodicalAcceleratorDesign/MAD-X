#######################################################################

# Makefile for MAD-X development version

#######################################################################

# parameters

FCP=-O3 -fno-second-underscore -funroll-loops
FCM=-O2 -fno-second-underscore -funroll-loops
FCDB=-g -O0 -fno-second-underscore
FC=g77
GCC_FLAGS=-g -Wall -fno-second-underscore
GCCP_FLAGS=-g -O3 -funroll-loops -fno-second-underscore
CC=gcc
LIBX=/usr/X11R6/lib/libX11.a
GLIB=/afs/cern.ch/group/si/slap/lib
GPUB=/afs/cern.ch/group/si/slap/bin


ifeq ($(OSTYPE),darwin)
# allows running of madx under Macinstosh System 10
# -fno-second-underscore  is old, do not use for more recent gnu compilers
# better no optimization (safer for 2.95.2 compiler, should be no problem for 3.1) 
# include headers for gxx11c
  GCCP_FLAGS=-g -Wall -I /usr/X11R6/include/
# use next to get older compiler
  CC=gcc2
# use next to get newer compiler
# CC=gcc3
# better use g77 without optimization
  FCP=-fno-second-underscore -funroll-loops
endif

default: madx

# files

twissp.o: twiss.f
	$(FC) $(FCP) -c -o twissp.o twiss.f

surveyp.o: survey.f
	$(FC) $(FCP) -c -o surveyp.o survey.f

utilp.o: util.f
	$(FC) $(FCP) -c -o utilp.o util.f

dynapp.o: dynap.f
	$(FC) $(FCP) -c -o dynapp.o dynap.f

madxnp.o: madxn.c madxu.c madxe.c madxc.c matchc.c sxf.c madx.h madxl.h madxd.h madxdict.h makethin.c c6t.c c6t.h
	$(CC) $(GCCP_FLAGS) -c -o madxnp.o madxn.c

madx: madxnp.o twissp.o surveyp.o orbfp.o emitp.o utilp.o matchp.o dynapp.o plotp.o ibsdbp.o trrunp.o gxx11.o gxx11c.o
	$(FC) $(FP) -o madx madxm.f madxnp.o twissp.o matchp.o ibsdbp.o plotp.o trrunp.o dynapp.o surveyp.o orbfp.o emitp.o utilp.o gxx11.o gxx11c.o $(LIBX) -lm -lc

ibsdbp.o: ibsdb.f
	$(FC) $(FCP) -c -o ibsdbp.o ibsdb.f

plotp.o: plot.f
	$(FC) $(FCP) -c -o plotp.o plot.f

trrunp.o: trrun.f
	$(FC) $(FCP) -c -o trrunp.o trrun.f

orbfp.o: orbf.fpp
	$(FC) $(FCP) -c -o orbfp.o orbf.fpp

emitp.o: emit.f
	$(FC) $(FCP) -c -o emitp.o emit.f

matchp.o: match.f
	$(FC) $(FCP) -c -o matchp.o match.f

gxx11.o: gxx11.f
	$(FC) $(FCP) -c -o gxx11.o gxx11.f

gxx11c.o: gxx11c.c
	$(CC) $(GCCP_FLAGS) -c -o gxx11c.o gxx11c.c

clean:
	rm -f *.o
	rm -f core
