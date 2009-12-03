#!/usr/bin/python
# THIS CODE NOT USED YET. INSTEAD MADBUILDPY.PL IS USED

import time
import optparse

# MISSING TIME MEASUREMENT



web = "/afs/cern.ch/user/n/nougaret/www/mad/"

class Repository:
    repoDir = ":gserver:isscvs.cern.ch:/local/reps/madx"
    def __init__(self):
        pass
    def checkout(self):
        command = "cvs -d " + repoDir + " checkout -r" + releaseTag + "madX" 
        os.system(command)

# first set environment variables for lf95 and NAG compilers
# this is necessary for the acron job which does not perform a login
# that would set the variables transparently.

def main():
    usage = "%prog [options]"
    parser = optparse.OptionParser(usage)
    parser.add_option('--release','-r',help="release tag",action="store_true")
    parser.add_option('--dev','-d',help="construct executable with Makefile_develop",action="store_true")
    parser.add_option('--nag','-n',help="construct executable with Makefile_nag",action="store_true")
    (options,args) = parser.parse_args()

    if not options.release:
        raise("except a release tag to be specified")
    
    os.environ['PATH'] = os.environ['PATH'] +\
                         ":/afs/cern.ch/sw/fortran/nag/f95.5.361/bin:/afs/cern.ch/sw/fortran/lahey/lf9562/bin"
    os.environ['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH'] + ":/afs/cern.ch/sw/fortran/lahey/lf9562/lib"

    makefiles = ['Makefile']

    if options.nag:
        os.environ['NAG95_ROOT'] = "/afs/cern.ch/sw/fortran/nag/f95.5.361" # flexlm license manager for NAG compiler
        os.environ['LM_LICENSE_FILE'] = "/afs/cern.ch/sw/fortran/nag/f95.5.361/license.dat"
        makefiles.append('Makefile_nag')

    if options.dev:
        makefiles.append('Makefile_develop')

    repo = Repository()
    repo.checkout()
    
    for m in makefiles:
        if m == Makefile:
            invocation = "Makefile f95=/opt/intel/Compiler/11.0/081/bin/ia32/ifort DEBUG=NO madx"
        elif m == Makefile_develop:
            invocation = "Makefile f95=lf95 DEBUG=YES madx" # the Lahey Fortran compiler
        elif m == Makefile_nag:
            invocation = "Makefile f95=f95 DEBUG = YES madx" # the NAG compiler

        output = "./makeResult"
        os.system(invocation+">& "+output)

if __name__ == "__main__":
    main()

