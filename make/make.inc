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

# compilers and linker specific settings
-include $(makedir)/compiler.$(CC)
-include $(makedir)/compiler.$(CXX)
-include $(makedir)/compiler.$(FC)
-include $(makedir)/linker.$(LD)

# system specific settings
-include $(makedir)/system.$(OSTYPE)

# project files
-include $(FILE_VER)
-include $(FILE_CPP)
-include $(FILE_C)
-include $(FILE_CXX)
-include $(FILE_F90)
-include $(FILE_LIB)
-include $(FILE_SYS)

# compilers, linker and depend rules
include $(makedir)/compiler.rules
include $(makedir)/linker.rules
include $(makedir)/depend.rules

# cleaning & debugging
include $(makedir)/clean.rules
include $(makedir)/infos.rules

# end of makefile
