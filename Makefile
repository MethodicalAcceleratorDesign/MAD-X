#######################################################################
#
# Makefile for MAD-X development version
#
#######################################################################

# compilers
CC=gcc
FC=g77
# NAG for testing
#f95=f95
# LF95 for production
f95=lf95

# default fortran compiler options
FCP=-O4 -fno-second-underscore -funroll-loops -I.

# alternative for development and debug
FCM=-O2 -fno-second-underscore -funroll-loops
FCDB=-g -O0 -fno-second-underscore

# default C compiler flag options
GCCP_FLAGS_MPARS=-g -O4 -funroll-loops -fno-second-underscore -D_CATCH_MEM
GCCP_FLAGS=$(GCCP_FLAGS_MPARS) -D_FULL

# alternative for development
GCC_FLAGS=-g -Wall -fno-second-underscore -D_CATCH_MEM -D_FULL

# NAG default f95 compiler options
#f95_FLAGS=-gline -g90 -c -C=all -maxcontin=100 -nan
# NAG alternative
#f95_FLAGS=-c -O4 -maxcontin=100 -w=unused
# LF95 default f95 compiler options
f95_FLAGS= --o1 --tp -c

# NAG f95 compiler options to compile f77 code
#FFLAGS77=-gline -g90 -c -maxcontin=100 -nan
# NAG f95 alternatives for development and debug
#FFLAGS77=-gline -g90 -c -maxcontin=100 -nan -ieee=full
#FFLAGS77=-g90 -c -O4 -maxcontin=100 -w=unused
# LF95 f95 compiler options to compile f77 code
FFLAGS77= --o1 --tp -c

# g77 link options
FP=-static

# NAG f95 link options
#f95_FOPT=
# LF95 f95 link options
f95_FOPT=-static

# libraries
#LIBX="-L/usr/X11R6/lib" -lX11 "-L/usr/lib/" -lgcc
LIBX="-L/usr/X11R6/lib" -lX11 "-L/usr/lib/" -ldl -lpthread

# NAG f95 lib extension
#LIBX_ext= -lgcc
# LF95 f95 lib extension
LIBX_ext=

ifeq ($(OSTYPE),darwin)
# allows running of madx under Macinstosh System 10
# -fno-second-underscore  is old, do not use for more recent gnu compilers
# include headers for gxx11c
  GCCP_FLAGS_MPARS=-g -O4 -funroll-loops -D_CATCH_MEM -I /usr/X11R6/include/
  GCCP_FLAGS=$(GCCP_FLAGS_MPARS) -D_FULL
  FP=
endif

default: madx

# dependencies of madxpf which combines the C-code
madxp.o: madxp.c madxn.c madxu.c madxe.c madxc.c matchc.c matchc2.c sxf.c makethin.c c6t.c madxreg.c madxreg.h madx.h madxl.h madxd.h madxdict.h c6t.h
	$(CC) $(GCCP_FLAGS_MPARS) -c madxp.c

madxpf.o: madxp.c madxn.c madxu.c madxe.c madxc.c matchc.c matchc2.c sxf.c makethin.c c6t.c madxreg.c madxreg.h madx.h madxl.h madxd.h madxdict.h c6t.h
	$(CC) $(GCCP_FLAGS) -c -o madxpf.o madxp.c

# fortran code dependencies on header files fi
twiss_f77.o twiss.o: twiss.F twiss0.fi twissa.fi twissl.fi twissc.fi twissotm.fi track.fi bb.fi name_len.fi twtrr.fi
util_f77.o util.o: util.F twiss0.fi twtrr.fi
dynap_f77.o dynap.o: dynap.F deltra.fi dyntab.fi wmaxmin0.fi tunes.fi
ibsdb_f77.o ibsdb.o: ibsdb.F ibsdb.fi name_len.fi physcons.fi
plot_f77.o plot.o: plot.F plot.fi plot_b.fi plot_c.fi plot_math.fi
sodd_f77.o sodd.o: sodd.F
trrun_f77.o trrun.o: trrun.F twiss0.fi name_len.fi track.fi bb.fi twtrr.fi
emit_f77.o emit.o: emit.F twiss0.fi bb.fi emit.fi twtrr.fi
match_f77.o match.o: match.F name_len.fi match.fi
touschek_f77.o touschek.o: touschek.F touschek.fi name_len.fi physcons.fi
resindex_f77.o resindex.o: resindex.F resindex.fi

# f90 dependencies
a_scratch_size.o: a_scratch_size.f90
b_da_arrays_all.o: a_scratch_size.o b_da_arrays_all.f90
c_dabnew.o: b_da_arrays_all.o c_dabnew.f90
d_lielib.o: c_dabnew.o d_lielib.f90
e_define_newda.o: d_lielib.o e_define_newda.f90
f_newda.o: e_define_newda.o f_newda.f90
g_newLielib.o: f_newda.o g_newLielib.f90
h_definition.o: g_newLielib.o h_definition.f90
i_tpsa.o: h_definition.o i_tpsa.f90
j_tpsalie.o: i_tpsa.o j_tpsalie.f90
k_tpsalie_analysis.o: j_tpsalie.o k_tpsalie_analysis.f90
l_complex_taylor.o: k_tpsalie_analysis.o l_complex_taylor.f90
m_real_polymorph.o: l_complex_taylor.o m_real_polymorph.f90
n_complex_polymorph.o: m_real_polymorph.o n_complex_polymorph.f90
o_tree_element.o: n_complex_polymorph.o o_tree_element.f90
Sa_extend_poly.o: o_tree_element.o Sa_extend_poly.f90
Sb_1_pol_template.o: Sa_extend_poly.o Sb_1_pol_template.f90
Sb_2_pol_template.o: Sb_1_pol_template.o Sb_2_pol_template.f90
Sb_sagan_pol_arbitrary.o: Sb_2_pol_template.o Sb_sagan_pol_arbitrary.f90
Sc_euclidean.o: Sb_sagan_pol_arbitrary.o Sc_euclidean.f90
Sd_frame.o: Sc_euclidean.o Sd_frame.f90
Se_status.o: Sd_frame.o Se_status.f90 a_def_all_kind.inc a_def_sagan.inc \
	a_def_user1.inc a_def_user2.inc a_def_element_fibre_layout.inc
Sf_def_all_kinds.o: Se_status.o Sf_def_all_kinds.f90
Sg_1_fitted.o: Sf_def_all_kinds.o Sg_1_fitted.f90
Sg_1_template_my_kind.o: Sg_1_fitted.o Sg_1_template_my_kind.f90
Sg_2_template_my_kind.o: Sg_1_template_my_kind.o Sg_2_template_my_kind.f90
Sg_sagan_wiggler.o: Sg_2_template_my_kind.o Sg_sagan_wiggler.f90
Sh_def_kind.o: Sg_sagan_wiggler.o Sh_def_kind.f90
Si_def_element.o: Sh_def_kind.o Si_def_element.f90
Sj_elements.o: Si_def_element.o Sj_elements.f90
Sk_link_list.o: Sj_elements.o Sk_link_list.f90
Sl_family.o: Sk_link_list.o Sl_family.f90
Sm_tracking.o: Sl_family.o Sm_tracking.f90
Sn_mad_like.o: Sm_tracking.o Sn_mad_like.f90
So_fitting.o: Sn_mad_like.o So_fitting.f90
Sp_keywords.o: So_fitting.o Sp_keywords.f90
madx_ptc_module.o: Sp_keywords.o madx_ptc_setcavs.o madx_ptc_tablepush.o madx_ptc_module.f90 
madx_ptc_track_run.o: madx_ptc_module.o madx_ptc_track_run.f90
madx_ptc_intstate.o: Sp_keywords.o madx_ptc_intstate.f90
madx_ptc_trackcavs.o: Sp_keywords.o madx_ptc_intstate.o  madx_ptc_setcavs.o madx_ptc_module.o madx_ptc_trackcavs.f90
madx_ptc_setcavs.o  : Sp_keywords.o madx_ptc_intstate.o  madx_ptc_setcavs.f90
madx_ptc_tablepush.o : madx_ptc_tablepush.f90
user2_photon.o: madx_ptc_track_run.o user2_photon.f90 photoni.inc
wrap.o: madx_ptc_module.o wrap.f90
run_madx.o: madx_ptc_module.o run_madx.f90
madx_main.o: run_madx.o madx_main.f90

# implicit rule to compile with C
%.o : %.c
	$(CC) $(GCCP_FLAGS) -c -o $(@) $<

# implicit rule to compile with f77. Append _f77 to distinguish from object code compiled with f95
%_f77.o : %.F
	$(FC) $(FCP) -c -o $(@) $<

# implicit rule to compile f77 code with f95
%.o : %.F
	$(f95) $(FFLAGS77) $<

# implicit rule to compile f90 code with f95
%.o : %.f90
	$(f95) $(f95_FLAGS) $<

#Parser only
mpars: madxm.F madxp.o
	$(FC) $(FP) -o mpars madxm.F madxp.o $(LIBX) -lm -lc

# madx_objectsf77: madxpf.o gxx11c.o  + all *.F except for gxx11ps.F timest.F timex.F (windows special & F90).
# Append f77 to distinguish from objects compiled with f95
madx_objectsf77 = madxpf.o gxx11c.o timel.o $(filter-out gxx11ps_f77.o madxp.o, $(patsubst %.F,%_f77.o,$(wildcard *.F)))

madx: $(madx_objectsf77) ;
	$(FC) $(FP) -o $@ $(madx_objectsf77) $(LIBX) -lgcc -lm -lc

# madx_objectsf95 all *.F without madxm.F, ptc_dummy.F & gxx11ps.F (windows special)
madx_objectsf95 = $(filter-out madxm.o ptc_dummy.o gxx11ps.o madxp.o, $(patsubst %.F,%.o,$(wildcard *.F)))
# madxdev_objects. All *.f90 , some c and F
madxdev_objects = $(patsubst %.f90,%.o,$(wildcard *.f90)) \
	madxpf.o gxx11c.o rplot.o \
	$(madx_objectsf95)
madxdev: $(madxdev_objects)
	$(f95) $(f95_FOPT) -o $@ $(madxdev_objects) $(LIBX) $(LIBX_ext)

clean:
	rm -f *.o
	rm -f *.g90
	rm -f *.mod
	rm -f core
	rm -f *~

info:
	@echo "-------------------------------------"
	@echo  Makefile for madX by Helmut Burkhardt
	@echo "-------------------------------------"
	@echo madx_objectsf77 = $(sort $(madx_objectsf77))
	@echo
	@echo madxdev_objects = $(sort $(madxdev_objects))
	@echo
	@echo default C compiler CC "    " = $(CC)
	@echo GCC_FLAGS "                " = $(GCC_FLAGS)
	@echo GCCP_FLAGS "               " = $(GCCP_FLAGS)
	@echo default Fortran compiler FC  = $(FC)
	@echo FFLAGS77 "                 " = $(FFLAGS77)
	@echo f95 "                      " = $(f95)
	@echo FCP "                      " = $(FCP)
	@echo FCM "                      " = $(FCM)
	@echo FCDB "                     " = $(FCDB)
	@echo LIBX "                     " = $(LIBX)
	@echo GLIB "                     " = $(GLIB)
	@echo GPUB "                     " = $(GPUB)
	@echo the OS is "                " = $(OS)
	@echo the OSTYPE is "            " = $(OSTYPE)
