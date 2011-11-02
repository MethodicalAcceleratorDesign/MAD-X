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

# must be first
include make/make.lib

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
# Platform settings
#

# auto detect
OSTYPE := $(shell uname -s)
OSVERS := $(shell uname -r)
OSARCH := $(shell uname -m)

# build directory
OBJDIR := $(CURDIR)/$(OSTYPE)$(ARCH)

###########
# Includes

# system specific settings (optional)
-include make/system.$(OSTYPE)

# compilers specific settings
include make/compiler.$(CC)
include make/compiler.$(FC)
include make/compiler.rules

# files dependencies
include make/files.c
include make/files.cpp
include make/files.f90

# cleaning & debugging
include make/clean.rules
include make/infos.rules

######################
# End of MAD makefile

