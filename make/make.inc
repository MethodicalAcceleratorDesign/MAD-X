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
include make/make.lib

# make config
include make/make.cfg

# compilers and linker specific settings
-include make/compiler.$(CC)
-include make/compiler.$(CXX)
-include make/compiler.$(FC)
-include make/linker.$(LD)

# system specific settings
-include make/system.$(OSTYPE)

# project files
-include $(FILE_VER)
-include $(FILE_CPP)
-include $(FILE_C)
-include $(FILE_CXX)
-include $(FILE_F90)
-include $(FILE_LIB)
-include $(FILE_SYS)

# compilers, linker and depend rules
include make/compiler.rules
include make/linker.rules
include make/depend.rules

# cleaning & debugging
include make/clean.rules
include make/infos.rules

# end of makefile
