#######################################################################
#
# Makefile for MAD-X development version
#
#######################################################################

# compilers
CC=gcc
FC=g77
f95=f95

# default fortran compiler options
FCP=-O3 -fno-second-underscore -funroll-loops -I.

# alternative for development and debug
FCM=-O2 -fno-second-underscore -funroll-loops
FCDB=-g -O0 -fno-second-underscore

# default C compiler flag options
GCCP_FLAGS=-g -O3 -funroll-loops -fno-second-underscore -D_CATCH_MEM

# alternative for development
GCC_FLAGS=-g -Wall -fno-second-underscore -D_CATCH_MEM

# default f95 compiler options
f95_FLAGS=-gline -g90 -c -C=all -maxcontin=24 -nan
# alternative
#f95_FLAGS=-c -O4 -maxcontin=24 -w=unused

# f95 compiler options to compile f77 code
#FFLAGS77=-gline -g90 -c -maxcontin=24 -nan -ieee=full
FFLAGS77=-gline -g90 -c -maxcontin=24 -nan
#FFLAGS77=-g90 -c -O4 -maxcontin=24 -w=unused

# libraries
LIBX="-L/usr/X11R6/lib" -lX11 "-L/usr/lib/" -lgcc

ifeq ($(OSTYPE),darwin)
# allows running of madx under Macinstosh System 10
# -fno-second-underscore  is old, do not use for more recent gnu compilers
# include headers for gxx11c
  GCCP_FLAGS=-g -O3 -funroll-loops -D_CATCH_MEM -I /usr/X11R6/include/
endif

default: madx

# dependencies of madxnp which combines the C-code
madxnp.o: madxn.c madxu.c madxe.c madxc.c matchc.c sxf.c madx.h madxl.h madxd.h madxdict.h makethin.c c6t.c c6t.h
	$(CC) $(GCCP_FLAGS) -c -o madxnp.o madxn.c

# fortran code dependencies on header files fi
twiss_f77.o twiss.o: twiss.F twiss0.fi twissa.fi twissl.fi twissc.fi twissotm.fi track.fi bb.fi name_len.fi twtrr.fi
util_f77.o util.o: util.F twiss0.fi twtrr.fi
dynap_f77.o dynap.o: dynap.F deltra.fi dyntab.fi wmaxmin0.fi tunes.fi
ibsdb_f77.o ibsdb.o: ibsdb.F ibsdb.fi name_len.fi physcons.fi
plot_f77.o plot.o: plot.F plot.fi plot_b.fi plot_c.fi plot_math.fi
trrun_f77.o trrun.o: trrun.F twiss0.fi name_len.fi track.fi bb.fi twtrr.fi 
emit_f77.o emit.o: emit.F twiss0.fi bb.fi emit.fi twtrr.fi
match_f77.o match.o: match.F name_len.fi match.fi 
touchek_f77.o touchek.o: touchek.F

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
Sa_extend_poly.o: n_complex_polymorph.o Sa_extend_poly.f90
Sb_1_pol_template.o: Sa_extend_poly.o Sb_1_pol_template.f90
Sb_2_pol_template.o: Sb_1_pol_template.o Sb_2_pol_template.f90
Sc_euclidean.o: Sb_2_pol_template.o Sc_euclidean.f90
Sd_frame.o: Sc_euclidean.o Sd_frame.f90
Se_status.o: Sd_frame.o Se_status.f90
Sf_def_all_kinds.o: Se_status.o Sf_def_all_kinds.f90 \
	a_def_all_kind.inc a_def_user1.inc a_def_user2.inc \
	a_def_element_fibre_layout.inc
Sg_0_fitted.o: Sf_def_all_kinds.o Sg_0_fitted.f90
Sg_1_template_my_kind.o: Sg_0_fitted.o Sg_1_template_my_kind.f90
Sg_2_template_my_kind.o: Sg_1_template_my_kind.o Sg_2_template_my_kind.f90
Sh_def_kind.o: Sg_2_template_my_kind.o Sh_def_kind.f90
Si_def_element.o: Sh_def_kind.o Si_def_element.f90
Sj_elements.o: Si_def_element.o Sj_elements.f90
Sk_link_list.o: Sj_elements.o Sk_link_list.f90
Sl_family.o: Sk_link_list.o Sl_family.f90
Sm_tracking.o: Sl_family.o Sm_tracking.f90
Sn_mad_like.o: Sm_tracking.o Sn_mad_like.f90
So_fitting.o: Sn_mad_like.o So_fitting.f90
Sp_keywords.o: So_fitting.o Sp_keywords.f90
madx_ptc_module.o: Sp_keywords.o madx_ptc_module.f90
wrap.o: madx_ptc_module.o wrap.f90

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

# madx_objects: madxnp.o gxx11c.o  + *.F  Append f77 to distinguish from objects compiled with f95
madx_objectsf77 = madxnp.o gxx11c.o $(patsubst %.F,%_f77.o,$(wildcard *.F))

# madx_objectsf95 all *.F without madxm.F, ptc_dummy.F
madx_objectsf95 = $(filter-out madxm.o ptc_dummy.o, $(patsubst %.F,%.o,$(wildcard *.F)))

madx: $(madx_objectsf77) ; 
	$(FC) $(FP) -o $@ $(madx_objectsf77) $(LIBX) -lm -lc

# madxdev_objects. All *.f90 , some c and F
madxdev_objects = madxm.o $(patsubst %.f90,%.o,$(wildcard *.f90)) \
	madxnp.o gxx11c.o epause.o timel.o usleep.o \
	$(madx_objectsf95)
madxdev: $(madxdev_objects)
	$(f95) $(FOPT) -o $@ $(madxdev_objects) $(LIBX) -lm -lc

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
