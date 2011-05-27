#######################################################################
#
# Makefile for MAD-X development version
#
#######################################################################
#
# Changes on 19.01.2011 by H.Renshall:
# add new flags f95_FLAGSP and Q for separate lf95 compilation of the
# polymorphic code which cannot use the --chk u (undefined) flag.
# Other compilers do not define these flags hence get same as before.
# Link in static libX11 libraries for ARCH=64 builds.
# End of 19.01.2011 changes
#
# C compiler
CC=gcc
# Fortran90 compiler
f95=g95
# architecture 32/64bit
ARCH=32
# Debugging compiler flags
DEBUG=NO
# Online version - mostly SDDS IO
ONLINE=NO
# Memory leak search version
MEMLEAKS=NO
# profiling version
PROFILE=NO
# Piotr's pluggins "root" etc
PLUGIN_SUPPORT=YES
# Alternative DA package in C++
NTPSA=NO

# Mac version
ONMAC=NO

# Temorary Fix for IFORT on Fedora 13/14
IFORTFIX=
# IFORTFIX=-no-ipo

OSTYPE = $(shell uname -s)

ifeq ($(findstring arwin, $(OSTYPE)),arwin)
  ONMAC=YES
#  f95=g95
  f95=gfortran
#  for darwin go now by default to 64 bit executables
  ARCH=64

#!!!!!!!!!!!!!!!!!!!!!!!!!!!
# DON'T TOUCH FOR DARWIN !!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DEBUG=NO
  MEMLEAKS=NO
  PROFILE=NO
endif

#######################################################################
# Compilers for MAD-X
#######################################################################

# CC=gcc

# FORTRAN90 compilers proven to work
#
#             Intel     Lahey   Andy Vaught      GNU       NAG
# Default     Linux                             Darwin
# f95=        ifort     lf95       g95         gfortran    f95

#######################################################################
# Compile flags
#######################################################################

ifeq ($(ARCH),32)
  M32= -m32
else
  ifeq ($(findstring arwin, $(OSTYPE)),arwin)
    M32= -m64
  endif
endif

# Production: C compiler flag options
 GCCP_FLAGS= -g $(M32) -funroll-loops -D_CATCH_MEM -D_WRAP_FORTRAN_CALLS -D_WRAP_C_CALLS -I. -D_FULL
# to turn off fatal error at memory overflow add -D_DONOTCATCHOVERFLOW

# Standard FORTRAN flags
 f95_FLAGS= -c -funroll-loops -I.


#######################################################################
# Link options
#######################################################################

LDOPT=-static $(M32) $(IFORTFIX)

#######################################################################
# Compiler special treatment
#######################################################################


ifeq ($(f95),lf95)
  ifeq ($(ARCH),32)
    f95_FLAGS= --o2 --tp -c -Wa,--32
  else
    f95_FLAGS= --o2 -c
  endif
endif

ifeq ($(f95),f95)
  ifeq ($(ARCH),32)
    f95_FLAGS= -gline -c -Wc,-m32 -Wl,-m32 -abi=32 -maxcontin=100 -ieee=full -D_NAG
    LDOPT= -Bstatic -Wl,-m32 -abi=32
  else
    f95_FLAGS= -gline -c -maxcontin=100 -ieee=full -D_NAG
    LDOPT= -Bstatic
  endif
endif

ifeq ($(f95),ifort)
  f95_FLAGS+= -assume noold_unit_star -D_INTEL_IFORT_SET_RECL
  ifeq ($(ARCH),32)
    f95_FLAGS+= $(M32) -fp-model precise
  endif
endif

ifeq ($(f95),g95)
  ifeq ($(ARCH),32)
    f95_FLAGS+= -Wa,--32 -fno-second-underscore
  else
    f95_FLAGS+= -Wa,--64 -fno-second-underscore
  endif
endif

ifeq ($(f95),gfortran)
  f95_FLAGS+= -fno-range-check
  ifeq ($(ARCH),32)
     f95_FLAGS+=$(M32)
  endif
endif

ifeq ($(DEBUG),YES)
  ifeq ($(f95),lf95)
    # Replace Makefile_develop
    # lf95 compiler options with severe lf95/Fujitsu flags
    # store in separate flags as polymorphic f90 cannot work with --chk u
    ifeq ($(ARCH),32)
      f95_FLAGSP= --info --f95 --lst -V -g --chk aesux  --ap --trace --trap --verbose
      f95_FLAGSQ= --info --f95 --lst -V -g --chk aesux  --ap --trace --trap --verbose
    else
      f95_FLAGSP= --info --f95 --lst -V -g --chk aefosu --ap --trace --trap --verbose
      f95_FLAGSQ= --info --f95 --lst -V -g --chk aefs   --ap --trace --trap --verbose
    endif

    GCCP_FLAGS+= -Wall -pedantic
  endif
  ifeq ($(f95),f95)
    # Replace Makefile_nag
    f95_FLAGS+= -C=all -nan
    GCCP_FLAGS+= -Wall -pedantic
  endif
  ifeq ($(f95),g95)
    # Replace Makefile_nag
    f95_FLAGS+= -ggdb3
    GCCP_FLAGS+= -Wall -pedantic -ggdb3
  endif
  ifeq ($(f95),gfortran)
    f95_FLAGS+= -Wall -pedantic
    GCCP_FLAGS+= -Wall -pedantic
  endif
else
  GCCP_FLAGS+= -O4
  # use level 4 optimization, except for lf95 or on MAC
  ifneq ($(f95),lf95)
    ifeq ($(ONMAC),NO)
      f95_FLAGS+= -O4
    endif
  endif
endif

ifeq ($(ARCH),32)
  LIBX= -L$(PWD)/lib -lX11 -lpthread -lstdc++ -lc -lgcc_eh
else
  ifeq ($(f95),lf95)
    LIBX= -L$(PWD)/lib64 -lX11 -lstdc++ -lgcc_eh
  else
    LIBX= -L$(PWD)/lib64 -lX11 -lpthread -lstdc++ -lc -lgcc_eh
  endif
endif

ifeq ($(findstring arwin, $(OSTYPE)),arwin)
  ONLINE=NO
  LIBX= -L/usr/X11R6/lib -lX11 -L/usr/lib -lpthread -lstdc++
endif

ifeq ($(MEMLEAKS),YES)
  ifeq ($(f95),f95)
    f95_FLAGS+= -C=all -mtrace=size,line
    GCCP_FLAGS+= -Wall -pedantic -D_MEM_LEAKS
    LDOPT+= -mtrace=size,line
  endif
endif

ifeq ($(ONLINE),YES)
  GCCP_FLAGS+= -D_ONLINE
  ifeq ($(f95),ifort)
    LIBX+= -limf
  endif
  ifeq ($(ARCH),32)
    LIBX+= -L$(PWD)/lib -lSDDS1c -lSDDS1 -lrpnlib -lmdbmth -lmdblib -lgsl -lz
  else
    LIBX+= -L$(PWD)/lib64 -lSDDS1c -lSDDS1 -lrpnlib -lmdbmth -lmdblib -lgsl -lz
  endif
endif

ifeq ($(PROFILE),YES)
  f95_FLAGS+= -pg
  GCCP_FLAGS+= -pg
  LDOPT+= -pg
endif

ifeq ($(PLUGIN_SUPPORT),YES)
  GCCP_FLAGS+= -DPLUGIN_SUPPORT
  ifeq ($(f95),lf95)
    LDOPT=-dynamic --export $(M32) -ldl
  endif
  ifeq ($(f95),g95)
    LDOPT=-rdynamic --export $(M32) -ldl
  endif
endif


ifeq ($(findstring arwin, $(OSTYPE)),arwin)
# allows running of madx under Macinstosh System 10
# include headers for gxx11c
  GCCP_FLAGS += -I /usr/X11R6/include/
  ifeq ($(f95),g95)
    f95_FLAGS= -c -funroll-loops -I. -fno-second-underscore $(M32) -O4
  endif
  ifeq ($(f95),gfortran)
    f95_FLAGS += $(M32) -O2
  endif
  LDOPT= $(M32)
endif

default: madx

# dependencies of madxp which combines the C-code
madxp.o: madxp.c madxn.c madxu.c aperture.c madxe.c madxc.c matchc.c matchc2.c sxf.c makethin.c c6t.c madxreg.c madxreg.h madx.h madxl.h madxd.h madxdict.h c6t.h matchptcknobs.h fortran_wrappers.h c_wrappers.h
	$(CC) $(GCCP_FLAGS) -c -o madxp.o madxp.c

# automatically generated code
fortran_wrappers.h fortran_wrappers.c fortran_prototypes.h fortran_wrappers_prototypes.h:
	python wrap_fortran_calls.py # creates fortran_wrappers.c fortran_wrappers.h fortran_prototypes.h fortran_wrappers_prototypes.h

c_wrappers.h c_wrappers.c:
	python wrap_C_calls.py	

fortran_wrappers.o: fortran_wrappers.c
	$(CC) $(GCCP_FLAGS) -c fortran_wrappers.c

c_wrappers.o: c_wrappers.c
	$(CC) $(GCCP_FLAGS) -c c_wrappers.c

matchptcknobs.o: matchptcknobs.h matchptcknobs.c madx.h

# fortran code dependencies on header files fi
util.o: util.f90
dynap.o: util.o dynap.f90
emit.o: util.o emit.f90
ibsdb.o: util.o ibsdb.f90
match.o: util.o match.f90
matchjc.o: util.o matchjc.f90
matchsa.o: util.o matchsa.f90
gxx11.o: util.o gxx11.f90
plot.o: util.o plot.f90
resindex.o: util.o resindex.f90
sodd.o: sodd.f90
survey.o: survey.f90
touschek.o: util.o touschek.f90
trrun.o: util.o trrun.f90
twiss.o: util.o twiss.f90

# f90 dependencies
a_scratch_size.o: a_scratch_size.f90
b_da_arrays_all.o: a_scratch_size.o b_da_arrays_all.f90
ifeq ($(NTPSA),YES)
  c_dabnew_berz.o: b_da_arrays_all.o c_dabnew_berz.f90
  c_tpsa_interface.o: c_dabnew_berz.o c_tpsa_interface.F90
  d_lielib.o: c_tpsa_interface.o d_lielib.f90
  h_definition.o: a_scratch_size.o c_dabnew_berz.o d_lielib.o h_definition.f90 a_def_frame_patch_chart.inc a_def_all_kind.inc a_def_sagan.inc a_def_element_fibre_layout.inc
  tpsa.o: tpsa.cpp tpsa.h
	$(CC) $(GCCP_FLAGS) -c -o tpsa.o tpsa.cpp
  TPSA= tpsa.o
  FILT_TP_OUT= c_dabnew.o
  FILT_TP_OUT_F90=
else
  c_dabnew.o: b_da_arrays_all.o c_dabnew.f90
  d_lielib.o: c_dabnew.o d_lielib.f90
  h_definition.o: a_scratch_size.o c_dabnew.o d_lielib.o h_definition.f90 a_def_frame_patch_chart.inc a_def_all_kind.inc a_def_sagan.inc a_def_element_fibre_layout.inc
  TPSA=
  FILT_TP_OUT= c_dabnew_berz.o
  FILT_TP_OUT_F90= c_tpsa_interface.o
endif
i_tpsa.o: h_definition.o i_tpsa.f90
	$(f95) $(f95_FLAGS)  $(f95_FLAGSQ) i_tpsa.f90
j_tpsalie.o: i_tpsa.o j_tpsalie.f90
	$(f95) $(f95_FLAGS)  $(f95_FLAGSQ) j_tpsalie.f90
k_tpsalie_analysis.o: j_tpsalie.o k_tpsalie_analysis.f90
	$(f95) $(f95_FLAGS)  $(f95_FLAGSQ) k_tpsalie_analysis.f90
l_complex_taylor.o: k_tpsalie_analysis.o l_complex_taylor.f90
	$(f95) $(f95_FLAGS)  $(f95_FLAGSQ) l_complex_taylor.f90
m_real_polymorph.o: l_complex_taylor.o m_real_polymorph.f90
	$(f95) $(f95_FLAGS)  $(f95_FLAGSQ) m_real_polymorph.f90
n_complex_polymorph.o: m_real_polymorph.o n_complex_polymorph.f90
	$(f95) $(f95_FLAGS)  $(f95_FLAGSQ) n_complex_polymorph.f90
o_tree_element.o: n_complex_polymorph.o o_tree_element.f90
	$(f95) $(f95_FLAGS)  $(f95_FLAGSQ) o_tree_element.f90
Sa_extend_poly.o: o_tree_element.o Sa_extend_poly.f90
Sb_sagan_pol_arbitrary.o: Sa_extend_poly.o Sb_sagan_pol_arbitrary.f90
Sc_euclidean.o: Sb_sagan_pol_arbitrary.o Sc_euclidean.f90
Sd_frame.o: Sc_euclidean.o Sd_frame.f90
Se_status.o: Sd_frame.o Se_status.f90 a_def_all_kind.inc a_def_sagan.inc \
	a_def_element_fibre_layout.inc
Sf_def_all_kinds.o: Se_status.o Sf_def_all_kinds.f90 a_def_worm.inc
Sg_sagan_wiggler.o: Sf_def_all_kinds.o Sg_sagan_wiggler.f90
Sh_def_kind.o: Sg_sagan_wiggler.o Sh_def_kind.f90
Si_def_element.o: Sh_def_kind.o Si_def_element.f90
Sk_link_list.o: Si_def_element.o Sk_link_list.f90
Sl_family.o: Sk_link_list.o Sl_family.f90
Sm_tracking.o: Sl_family.o Sm_tracking.f90
Sma0_beam_beam_ptc.o: Sm_tracking.o Sma0_beam_beam_ptc.f90
Sma_multiparticle.o: Sma0_beam_beam_ptc.o Sma_multiparticle.f90
Sn_mad_like.o: Sma_multiparticle.o Sn_mad_like.f90
So_fitting.o: Sn_mad_like.o So_fitting.f90
Sp_keywords.o: So_fitting.o Sp_keywords.f90
Spb_fake_gino_sub.o: Sp_keywords.o Spb_fake_gino_sub.f90
Sq_orbit_ptc.o: Sp_keywords.o Sq_orbit_ptc.f90
Sqb_accel_ptc.o: Sq_orbit_ptc.o Sqb_accel_ptc.f90
Sr_spin.o: Sqb_accel_ptc.o Sr_spin.f90
Sra_fitting.o: Sr_spin.o Sra_fitting.f90
madx_ptc_module.o: util.o Sra_fitting.o madx_ptc_setcavs.o madx_ptc_knobs.o madx_ptc_module.f90
St_pointers.o: Sp_keywords.o madx_ptc_module.o St_pointers.f90
madx_ptc_track_run.o: util.o Sp_keywords.o madx_ptc_module.o madx_ptc_track_run.f90
madx_ptc_intstate.o: Sp_keywords.o madx_ptc_intstate.f90
madx_ptc_trackcavs.o: util.o Sp_keywords.o madx_ptc_intstate.o  madx_ptc_setcavs.o madx_ptc_module.o madx_ptc_trackcavs.f90
madx_ptc_setcavs.o: Sp_keywords.o madx_ptc_intstate.o madx_ptc_setcavs.f90
madx_ptc_script.o: util.o Sp_keywords.o madx_ptc_script.f90
madx_ptc_knobs.o: util.o Sp_keywords.o madx_ptc_intstate.o madx_ptc_knobs.inc madx_ptc_knobs.f90
madx_ptc_eplacement.o: util.o Sp_keywords.o madx_ptc_intstate.o madx_ptc_module.o madx_ptc_eplacement.f90
madx_ptc_normal.o: madx_ptc_module.o madx_ptc_normal.f90
madx_ptc_twiss.o: util.o madx_ptc_module.o madx_ptc_setcavs.o madx_ptc_knobs.o madx_ptc_distrib.o madx_ptc_knobs.inc madx_ptc_distrib.inc madx_ptc_twiss.f90
madx_ptc_distrib.o: util.o madx_ptc_module.o madx_ptc_distrib.inc madx_ptc_distrib.f90

wrap.o: madx_ptc_module.o  madx_ptc_intstate.o \
	madx_ptc_normal.o madx_ptc_twiss.o madx_ptc_distrib.o \
	madx_ptc_setcavs.o madx_ptc_trackcavs.o \
	madx_ptc_knobs.o \
	madx_ptc_script.o St_pointers.o ptc_export_xml.o madx_ptc_eplacement.o \
	wrap.f90
user2_photon.o: madx_ptc_track_run.o user2_photon.f90 photoni.inc
run_madx.o: madx_ptc_module.o run_madx.f90
madx_main.o: run_madx.o madx_main.F90

# implicit rule to compile with C
%.o : %.c
	$(CC) $(GCCP_FLAGS) -c -o $(@) $<

# implicit rule to compile f90 code with f95
%.o : %.f90
	$(f95) $(f95_FLAGS) $(f95_FLAGSP) $<

# implicit rule to compile f90 code with f95
%.o : %.F90
	$(f95) $(f95_FLAGS) $(f95_FLAGSP) $<

# madx_objects  = $(filter-out gxx11psc.o , $(patsubst %.c,%.o,$(wildcard *.c)))
madx_objects = madxp.o gxx11c.o matchptcknobs.o rplot.o fortran_wrappers.o c_wrappers.o $(TPSA)
madx_objects += $(filter-out gxx11ps.o $(FILT_TP_OUT), $(patsubst %.f90,%.o,$(wildcard *.f90)))
madx_objects += $(filter-out $(FILT_TP_OUT_F90),$(patsubst %.F90,%.o,$(wildcard *.F90)))

madx: $(madx_objects)
	$(f95) $(LDOPT) -o madx $(madx_objects) $(LIBX)
#	strip madx

x:
#	$(f95) $(f95_FLAGS) madx_ptc_setcavs.f90
	$(f95) $(f95_FLAGS) madx_ptc_twiss.f90
	$(f95) $(f95_FLAGS) madx_ptc_knobs.f90
	$(f95) $(LDOPT) -o madx $(madx_objects) $(LIBX)

clean:
	rm -f *.o
	rm -f *.g90
	rm -f *.mod
	rm -f *.lst
	rm -f core
	rm -f *~
	python wrap_fortran_calls.py --clean
	python wrap_C_calls.py --clean

info:
	@echo default C compiler CC "    " = $(CC)
	@echo GCC_FLAGS "                " = $(GCC_FLAGS)
	@echo f95 "                      " = $(f95)
	@echo LIBX "                     " = $(LIBX)
	@echo the OS is "                " = $(OS)
	@echo the OSTYPE is "            " = $(OSTYPE)
	@echo madx_objects "             " = $(madx_objects)
