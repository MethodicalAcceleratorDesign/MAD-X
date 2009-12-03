#!/usr/bin/python
# THIS CODE NOT USED YET. INSTEAD MADBUILDPY.PL IS USED

import time
import optparse
import os
import re
import shutil

# MISSING TIME MEASUREMENT

web = "/afs/cern.ch/user/n/nougaret/www/mad/"
currentDir = os.getcwd()
extractDir = currentDir+'/MadCvsExtract'

class Repository:
    repoDir = ":gserver:isscvs.cern.ch:/local/reps/madx"
    def __init__(self):
        pass
    def checkout(self,releaseTag):
        command = "cvs -d " + Repository.repoDir + " checkout -r" + releaseTag + " madX"
        os.system(command)

# first set environment variables for lf95 and NAG compilers
# this is necessary for the acron job which does not perform a login
# that would set the variables transparently.

def main():
    usage = "%prog [options]"
    parser = optparse.OptionParser(usage)
    parser.add_option('--release','-r',help="release tag",dest="releaseTag")
    parser.add_option('--dev','-d',help="construct executable with Makefile_develop",action="store_true")
    parser.add_option('--nag','-n',help="construct executable with Makefile_nag",action="store_true")
    (options,args) = parser.parse_args()

    if not options.releaseTag:
        raise("except a release tag to be specified")
    else: # is the release tag well formed?
        releasePattern = re.compile(r'^madX\-\d+_\d{2}_\d{2}$')
        if not releasePattern.match(options.releaseTag):
            raise("release tag is ill-formed: it should be like 'madX-4_01_01' instead of '"+options.releaseTag+"'")

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

    # if the directory in which to extract the CVS does not yet exist, then create it
    if not os.path.exists(extractDir):
        os.mkdir(extractDir)
    os.chdir(extractDir)
    repo = Repository()
    repo.checkout(options.releaseTag)
    os.chdir('./madX')
 
    for m in makefiles:
        # do a make clean and remove any preexisting target
        if os.path.exists('./madx'):
            os.remove('./madx')
        os.system('make clean')
        
        print("now to compile "+m+" in " + os.getcwd())
        if m == 'Makefile':
            invocation = "make f95=/opt/intel/Compiler/11.0/081/bin/ia32/ifort DEBUG=NO madx"
        elif m == 'Makefile_develop':
            invocation = "make f95=lf95 DEBUG=YES madx" # the Lahey Fortran compiler
        elif m == 'Makefile_nag':
            invocation = "make f95=f95 DEBUG = YES madx" # the NAG compiler

        output = "./makeResult"
        #os.system(invocation)
        os.system(invocation+">& "+output)
        # if the target madx is absent, compilation failed
        if os.path.exists('./madx'):
            # success or warning
            shutil.copyfile('./madx','./madx_'+m)
            print("succesfully compiled for "+m)
        else: # failure
            print("failed to compile for "+m)

if __name__ == "__main__":
    main()

