#######################################################################

# Makefile for MAD-X development version

#######################################################################

# parameters

FCP=-O3 -fno-second-underscore -funroll-loops
FCM=-O2 -fno-second-underscore -funroll-loops
FCDB=-g -O0 -fno-second-underscore
FC=g77
GCC_FLAGS=-g -Wall -fno-second-underscore
GCCP_FLAGS=-O3 -funroll-loops -fno-second-underscore
CC=gcc
LIBX=/usr/X11R6/lib/libX11.a
GLIB=/afs/cern.ch/group/si/slap/lib
GPUB=/afs/cern.ch/group/si/slap/bin

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
