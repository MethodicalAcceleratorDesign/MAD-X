#!/usr/bin/python

import re
import os
import sys
import traceback
import shutil

import socket

# public:

def releaseForWindows(tag): # to be called by MadTrigRelease.pl
    try:
        client = windowsReleaseClient()
        client.release(tag)
        return 'success'
    except:
        return 'failure'


        
# private:

class windowsReleaseClient:

    def __init__(self):
        pass
    
    def extractCVS(self,tag): # populate /user/nougaret/MAD-X-WINDOWS

        currentDir = os.getcwd()

        os.chdir('/user/nougaret/MAD-X-WINDOWS')

        shutil.rmtree('madX',ignore_errors=True) # clean-up
        cmd = 'cvs -d :gserver:isscvs.cern.ch:/local/reps/madx '+\
                  'checkout -r '+tag+' '+\
                  'madX'
        print cmd
        os.system(cmd)
        os.chdir(currentDir) # back to initial location

    def remoteCompile(self):
        thisLinuxHost = socket.gethostname()
        windowsHost = 'abpc10788'             # The remote host
        toWindowsHostPort = 7070              # The same port as used by the server
        #fromWindowsHostPort = 7071
        s = None
        for res in socket.getaddrinfo(windowsHost, toWindowsHostPort,\
                                      socket.AF_INET, socket.SOCK_STREAM):
            af, socktype, proto, canonname, sa = res
            try:
                s = socket.socket(af, socktype, proto)
            except socket.error, msg:
                s = None
                continue
            try:
                s.connect(sa)
            except socket.error, msg:
                s.close()
                s = None
                continue
            break
        if s is None:
            print 'could not open socket'
            sys.exit(1)
        s.send(thisLinuxHost+' asks: Compile MAD for Windows!')
        data = s.recv(1024) # blocking on reception (what about the time out???)
        s.close()
        print 'Received', repr(data)
    

    def checkCompileOutcome(self):
        pass

    def generateHtmlOutput(self):
        pass

    def release(self,tag):
        # first extract the CVS on NFS
        self.extractCVS(tag)
        # then remote invoke compilation on the Windows machine
        self.remoteCompile()
        # control the outcome of the compilation
        self.checkCompileOutcome()
        # generate HTML information and link to the executable
        self.generateHtmlOutput()

        
if __name__ == '__main__':

    if len(sys.argv) != 2:
        raise('expect exactly one argument: release tag') 
    else:
        tag = sys.argv[1]
        print "/"+tag+"/"
        pattern = re.compile(r'^madX\-(\d+)_(\d+)_(\d+).+$') # pro, dev ok
        m = pattern.match(tag)
        if not m:
            raise('argument release-tag is ill-formed')
    
        outcome = releaseForWindows(tag) # the main function
        print outcome
