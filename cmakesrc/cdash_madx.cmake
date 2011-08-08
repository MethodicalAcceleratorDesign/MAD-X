# This is a script for testing the source and submitting
# your results to a common server (cdash). This currently
# only works inside CERN. You can see the server at:
# http://137.138.26.237/cdash/index.php?project=Mad-X
# 
# How to:
#  - check out the source from svn in a temporary directory, from the url 
#     http://svnweb.cern.ch/guest/madx/trunk/madX/
#  - set the CTEST_SOURCE_DIRECTORY accordingly below
#  - set a useful test site name (e.g. "your name"."machine type")
#  - set a useful build name (e.g. architecture, os and compiler used)
#  - Run this script with the command:
#          ctest -S cdash_madx.cmake

# Necessary edits:
SET(CTEST_SITE "myname.maymachine")
set(CTEST_BUILD_NAME "SLC5-64bit-gfortran")
# Your source should be checked out from svn into this directory:
set(CTEST_SOURCE_DIRECTORY "/path/to/source/madX/")
# and compilation will be done in this directory (no need to edit, must be unique for every configuration):
set(CTEST_BINARY_DIRECTORY "${CTEST_SOURCE_DIRECTORY}/build")

# Optional edits:
ctest_start(Experimental)
set(cfg_options
 -DCMAKE_BUILD_TYPE=Release
 )

# Do not edit (unless you know cmake/ctest):

set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_UPDATE_COMMAND "svn")

ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
 
ctest_update()
 
ctest_configure(OPTIONS "${cfg_options}")
ctest_build() 
ctest_test()
# coverage test doesn't work at the moment..
#ctest_coverage()
ctest_submit()
