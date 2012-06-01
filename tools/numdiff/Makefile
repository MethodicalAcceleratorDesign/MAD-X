# |
# o---------------------------------------------------------------------o
# |
# | Numdiff makefile
# |
# o---------------------------------------------------------------------o
# |
# | Methodical Accelerator Design
# |
# | Copyright (c) 2011+ CERN, mad@cern.ch
# |
# | For more information, see http://cern.ch/mad
# |
# o---------------------------------------------------------------------o
# |
# | $Id$
# |

# For makefile documentation, please read make/README
# For information and bug report, please contact mad@cern.ch
#

###################
# Project settings

PROJECT := numdiff

#################
# Build settings
#

# architecture bit: detect/32/64 (default is detect)
ARCH    := detect

# debugging mode: yes/no (default is no)
DEBUG   := no

# profiling mode: yes/no (default is no)
PROFILE := no

#############################
# Compilers/Linkers settings
# see make/compiler.* for supported compilers
# GNU=yes   sets CC=gcc,     CXX=g++,     FC=gfortran (default)
# Intel=yes sets CC=icc/icl, CXX=icc/icl, FC=ifort    (use icl on Windows)

# C compiler (default is gcc)
CC  := gcc

# C++ compiler (default is g++)
CXX := g++

# Fortran compiler (default is gfortran)
FC  := gfortran

# Linker (default is fortran compiler, deferred)
LD     = $(CC)
LDNAME = $(CCNAME)

####################
# Includes settings

# makefiles to include (in this order, all optional)
FILE_PRE  := Makefile_pre          # user's preprocessing extra stuff
FILE_VER  := VERSION               # setup project VERSION and VERSION_DATE
FILE_CPP  := Makefile_cpp          # setup project specific CPPFLAGS defines
FILE_C    := Makefile_c            # setup project CC_SRC and CC_HDR
FILE_CXX  := Makefile_cxx          # setup project CXX_SRC and CXX_HDR
FILE_F90  := Makefile_f90          # setup project FC_SRC
FILE_LIB  := Makefile_lib          # setup project LIBS paths and libs
FILE_SYS  := Makefile_sys          # setup project System specific stuffs
FILE_TEST := Makefile_test         # setup project tests
FILE_POST := Makefile_post         # user's postprocessing extra stuff

####################
# Makefile includes

makedir := ../../make
include $(makedir)/make.inc

# end of makefile
