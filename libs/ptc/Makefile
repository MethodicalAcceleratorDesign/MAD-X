# |
# o---------------------------------------------------------------------o
# |
# | PTC makefile
# |
# o---------------------------------------------------------------------o
# |
# | Polymorphic Tracking Code from Etienne Forest
# |
# | For more information, see http://cern.ch/mad
# |
# o---------------------------------------------------------------------o
# |
# | $Id$
# |

# For makefile documentation, please read make/README
# For information and bug report, please contact mad@cern.ch

###################
# Project settings

PROJECT := libptc

#################
# Build settings
#

# architecture bit: detect/32/64 (default is detect)
ARCH    := detect

# debugging mode: yes/no (default is no)
DEBUG   := no

# profiling mode: yes/no (default is no)
PROFILE := no

# link with static libs: yes/no (default is no)
STATIC  := no

# make shared lib: yes/no (default is yes)
SHARED  := yes

#############################
# Compilers/Linkers settings
# see make/compiler.* for supported compilers
# GNU=yes   sets CC=gcc,     CXX=g++,     FC=gfortran (default)
# Intel=yes sets CC=icc/icl, CXX=icc/icl, FC=ifort    (use icl on Windows)

# Fortran compiler (default is gfortran)
FC  := gfortran

# Linker (default is Fortran compiler, deferred)
LD   = $(FC)

# Archiver (default is ar)
AR  := ar

####################
# Includes settings

FILE_F90 := Makefile_f90

####################
# Makefile includes

makedir := ../../make
include $(makedir)/make.inc

# end of makefile
