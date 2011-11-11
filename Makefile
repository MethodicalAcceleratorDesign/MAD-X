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

# For this makefile documentation, please read make/README
# For information and bug report, please contact mad@cern.ch

###################
# Project settings

PROJECT := madx

# online version - includes SDDS I/O (default is no)
ONLINE    := no

# alternative DA package in C++ (default is yes)
NTPSA     := yes

# memory leaks debug (default is no)
MEMLEAKS  := no

#################
# Build settings
#

# architecture 32/64bit (default is 64)
ARCH    := 64

# debugging mode (default is no)
DEBUG   := no

# profiling mode (default is no)
PROFILE := no

# link with static libs (default is no)
STATIC  := no

# plugin dynamic loading (default is no)
PLUGIN  := no

#############################
# Compilers/Linkers settings
#

# C compiler (default is gcc)
CC  := gcc

# C++ compiler (default is g++)
CXX := g++

# Fortran compiler (default is ifort)
FC  := ifort

# Linker (default is fortran compiler)
LD  := $(FC)

####################
# Includes settings

# makefiles to include (all optional)
FILE_VER := VERSION               # setup project VERSION and VERSION_DATE
FILE_CPP := Makefile_cpp          # setup project specific CPPFLAGS defines
FILE_C   := Makefile_c            # setup project CC_SRC and CC_HDR
FILE_CXX := Makefile_cxx          # setup project CXX_SRC and CXX_HDR
FILE_F90 := Makefile_f90          # setup project FC_SRC
FILE_LIB := Makefile_lib          # setup project LIBS paths and libs
FILE_SYS := Makefile_sys          # setup project System specific stuffs
FILE_USR := Makefile_usr          # user's extra stuffs

####################
# Makefile includes

makedir := make
include $(makedir)/make.inc

# end of makefile
