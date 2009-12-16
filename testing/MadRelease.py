#!/usr/bin/python

# TODO handle dev versions as well

import os
import shutil
import socket
import optparse
import re

cvsDir = ':gserver:isscvs.cern.ch:/local/reps/madx'
    
def main():
    # expect this script to run on pcslux99
    host = socket.gethostname()
    if not host == 'pcslux99':
        raise("expect this script to run on pcslux99")
    
    # get the release _compulsory_ option
    usage = "%prog options"
    parser = optparse.OptionParser(usage)
    parser.add_option("-r","--release",help="the release tag in the form madX_xx_yy where xx and yy are two-digit numbers",dest="release")
    (options,args) = parser.parse_args()
    if not options.release:
        raise("expect a release to be passed as compulsory option")
    # check well formedness of the release tag
    pattern = re.compile(r'^madX_\d{2}\-\d{2}$')
    m = pattern.match(options.release)
    if not m:
        raise("expect release tag to be of the form madX_xx-yy where xx and yy are two-digit numbers")
    

    

    
    initialDir = os.getcwd()
    checkout = 'cvs -d '+cvsDir+' checkout madX'
    try:
        shutil.rmtree('./madX') # clean-up the directory
    except:
        pass
    os.system(checkout)
    os.chdir('./madX')
    make = 'make'
    os.system(make)
    
    # now copy to the AFS web-folder

    # back to the inital dir
    os.chdir(initialDir)


if __name__ == "__main__":
    main()
