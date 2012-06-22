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

# make utilities (no rules!)
include $(makedir)/make.lib

# make config (no rules!)
include $(makedir)/make.cfg

# user's preprocessing extra stuff
-include $(FILE_PRE)

# actions
include $(makedir)/make.act

# system
$(call exists, $(makedir)/system.$(OSTYPE), \
unable to include settings for system OSTYPE=$(OSTYPE) (make [options] infoinc); \
operating system not supported)
include $(makedir)/system.$(OSTYPE)

# C compiler
ifneq ($(CCNAME),)
$(call exists, $(makedir)/compiler.$(CCNAME), \
unable to include settings for compiler CC=$(CCNAME) (make [options] infoinc); \
check CC or CCNAME for correct compiler configuration file)
include $(makedir)/compiler.$(CCNAME)
endif

# C++ compiler
ifneq ($(CXXNAME),)
$(call exists, $(makedir)/compiler.$(CXXNAME), \
unable to include settings for compiler CXX=$(CXXNAME) (make [options] infoinc); \
check CXX or CXXNAME for correct compiler configuration file)
include $(makedir)/compiler.$(CXXNAME)
endif

# Fortran compiler
ifneq ($(FCNAME),)
$(call exists, $(makedir)/compiler.$(FCNAME), \
unable to include settings for compiler FC=$(FCNAME) (make [options] infoinc); \
check FC or FCNAME for correct compiler configuration file)
include $(makedir)/compiler.$(FCNAME)
endif

# linker
ifneq ($(LDNAME),)
$(call exists, $(makedir)/linker.$(LDNAME), \
unable to include settings for linker LD=$(LDNAME) (make [options] infoinc); \
check LD or LDNAME for correct linker configuration file)
include $(makedir)/linker.$(LDNAME)
endif

# archiver
ifneq ($(ARNAME),)
$(call exists, $(makedir)/linker.$(ARNAME), \
unable to include settings for archiver AR=$(ARNAME) (make [options] infoinc); \
check AR or ARNAME for correct linker configuration file)
include $(makedir)/linker.$(ARNAME)
endif

# tester
ifneq ($(NDNAME),)
$(call exists, $(makedir)/tester.$(NDNAME), \
unable to include settings for tester ND=$(NDNAME) (make [options] infoinc); \
check ND or NDNAME for correct tester configuration file)
include $(makedir)/tester.$(NDNAME)
endif

# project files
-include $(FILE_VER)
-include $(FILE_CPP)
-include $(FILE_C)
-include $(FILE_CXX)
-include $(FILE_F90)
-include $(FILE_LIB)
-include $(FILE_SYS)
-include $(FILE_TEST)

# compilers, linker, depend and test rules
include $(makedir)/compiler.rules
include $(makedir)/linker.rules
include $(makedir)/tester.rules
include $(makedir)/depend.rules

# cleaning & debugging
include $(makedir)/clean.rules
include $(makedir)/info.rules

# user's postprocessing extra stuff
-include $(FILE_POST)

# end of makefile
