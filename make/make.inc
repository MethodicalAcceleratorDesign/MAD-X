# |
# o---------------------------------------------------------------------o
# |
# | MAD makefile - includes
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

# make utilities
include $(makedir)/make.lib

# make config
include $(makedir)/make.cfg

# user's preprocessing extra stuff
-include $(FILE_PRE)

# make actions
include $(makedir)/make.act

# project files
-include $(FILE_VER)
-include $(FILE_CPP)

# c files
ifneq ($(CCNAME),)
$(call exists, $(makedir)/compiler.$(CCNAME), \
unable to include settings for compiler CC=$(CCNAME) (make [options] infoinc); \
check CC or CCNAME for correct compiler configuration file)
include $(makedir)/compiler.$(CCNAME)
-include $(FILE_C)
endif

# c++ files
ifneq ($(CXXNAME),)
$(call exists, $(makedir)/compiler.$(CXXNAME), \
unable to include settings for compiler CXX=$(CXXNAME) (make [options] infoinc); \
check CXX or CXXNAME for correct compiler configuration file)
include $(makedir)/compiler.$(CXXNAME)
-include $(FILE_CXX)
endif

# fortran files
ifneq ($(FCNAME),)
$(call exists, $(makedir)/compiler.$(FCNAME), \
unable to include settings for compiler FC=$(FCNAME) (make [options] infoinc); \
check FC or FCNAME for correct compiler configuration file)
include $(makedir)/compiler.$(FCNAME)
-include $(FILE_F90)
endif

# lib, sys files
-include $(FILE_LIB)
-include $(FILE_SYS)

# linker files
ifneq ($(LDNAME),)
$(call exists, $(makedir)/linker.$(LDNAME), \
unable to include settings for linker LD=$(LDNAME) (make [options] infoinc); \
check LD or LDNAME for correct linker configuration file)
include $(makedir)/linker.$(LDNAME)
endif

# system files
-include $(makedir)/system.$(OSTYPE)

# tests files
ifneq ($(NDNAME),)
$(call exists, $(makedir)/tester.$(NDNAME), \
unable to include settings for tester ND=$(NDNAME) (make [options] infoinc); \
check ND or NDNAME for correct tester configuration file)
include $(makedir)/tester.$(NDNAME)
-include $(FILE_TEST)
endif

# user's postprocessing extra stuff
-include $(FILE_POST)

# compilers, linker, depend and test rules
include $(makedir)/compiler.rules
include $(makedir)/linker.rules
include $(makedir)/tester.rules
include $(makedir)/depend.rules

# cleaning & debugging
include $(makedir)/clean.rules
include $(makedir)/info.rules

# end of makefile
