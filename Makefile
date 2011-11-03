# |
# o---------------------------------------------------------------------o
# |
# | MAD makefile
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

###################
# Project Settings

PROJECT := madx
VERSION := 5.00.08
VERDATE := 2011-11-01

#
# Build settings
#

# architecture 32/64bit (default is 64)
ARCH      := 64

# debugging mode (default is no)
DEBUG     := no

# profiling mode (default is no)
PROFILE   := no

# plugging support (default is no)
PLUGIN    := no

# online version - includes SDDS I/O (default is no)
ONLINE    := no

# alternative DA package in C++ (default is yes)
NTPSA     := yes

# memory leaks debug (default is no)
MEMLEAKS  := no

#
# Compilers settings
#

# C compiler (default is gcc)
CC  := gcc

# C++ compiler (default is g++)
CXX := g++

# Fortran compiler (default is ifort)
FC  := ifort

#
# Files settings
#

FILES_C   := Makefiles_c
FILES_CXX := Makefiles_cpp
FILES_F90 := Makefiles_f90

#
# Platform settings
#

# auto detect
OSTYPE := $(shell uname -s)
OSVERS := $(shell uname -r)
OSARCH := $(shell uname -m)

# build directory
OBJDIR := $(OSTYPE)$(ARCH)

# non-build targets
NOBUILD := clean% info%

###########
# Includes

# make utilities
include make/make.lib

# files selection
-include $(FILES_C)            # setup CC_SRC,   CC_HDR
-include $(FILES_CXX)          # setup CXX_SRC,  CXX_HDR 
-include $(FILES_F90)          # setup FC_SRC

# compilers and linker specific settings
include make/compiler.$(CC)
include make/compiler.$(CXX)
include make/compiler.$(FC)
-include make/linker.$(LD)      # optional

# system specific settings
-include make/system.$(OSTYPE)  # optional

# compilers, linker and depend rules
include make/compiler.rules
include make/linker.rules
include make/depend.rules

# cleaning & debugging
include make/clean.rules
include make/infos.rules

# end of makefile
