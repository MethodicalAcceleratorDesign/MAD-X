#!/usr/bin/python

import os
import re

import Notify

global representative
representative = "madxd.h"

def byDecreasingReleaseNumber(a,b): # passed to the Python sort list function
    aNumbers = a.split('_')
    bNumbers = b.split('_')
    counter = 0
    result = 0
    for aNumber in aNumbers:
        if int(aNumber) < int(bNumbers[counter]):
            result = 1
            return result
        if int(aNumber) > int(bNumbers[counter]):
            result = -1
            return result
        counter += 1
    
    return result

class MadTestTrigger:
    
    def __init__(self):
        self.whatToDo = "don't know"
        self.release = "don't know"
        self.debugMode = False

    def setDebugMode(self,mode):
        self.debugMode=mode

    def getWhatToDo(self):
        return self.whatToDo

    def getRelease(self):
        return self.release

    def run(self,rootDir):
        extractDir = self.extractCVS(rootDir)
        [ releases, tests ]= self.collectReleases(extractDir)
        # sort releases and tests by increasing numbers
        
        releases.sort(byDecreasingReleaseNumber) # in place sorting, with supplied comparison function
        tests.sort(byDecreasingReleaseNumber)

        # debug:
        #for release in releases:
        #    print "release " + release

        lastRelease = releases[0]
        beforeLastRelease = releases[1]

        lastTest = tests[0]

#        print "lastTest="+lastTest
        # completion
        self.release = "madX-"+lastRelease
        if lastTest == lastRelease:
            self.whatToDo = "do-nothing"
        else:
            # SHOULD TAG THE CVS HERE (if mode != debug)
            # should also record work
            global representative
            newTestTag = "test-"+lastRelease
            command = "cvs tag " + newTestTag + " " + representative
            if not self.debugMode:
                os.system(command)
            else:
                print "skip command '"+command+"' in debug-mode"
            self.whatToDo = "run-test"
        
    def extractCVS(self,rootDir):
        # do we need to check-out or could we rely on FesaBuild.pl instead?
        cvsRootDir = ":gserver:isscvs.cern.ch:/local/reps/madx"
        command = "cvs -d "+cvsRootDir+" checkout madX"
        os.chdir(rootDir)
        extractDir = "MadCvsExtract_trig_test"
        try:
            os.mkdir(extractDir)
        except OSError:
#            print "OS error (should be instance errno 17)"
            os.system('rm -rf ./'+extractDir) # dangerous
            os.mkdir(extractDir)
        os.chdir(extractDir)
        os.system(command)
        os.chdir('./madX')
        theDir = os.getcwd()
        return theDir

    def collectReleases(self,extractDir):

        releases = []
        tests = []

        releaseRevision = {}
        testRevision = {}
        
        os.chdir(extractDir)
        global representative
        log = os.system('cvs log '+representative+' > tempfile')
        tempfile = open('tempfile','r')
        lines = tempfile.readlines()

        for line in lines:

            pattern = re.compile(r'^[\s\t]*madX\-(\d+)_(\d+)_(\d+)[\s\t]*:[\s\t]*([\d\.]+)[\s\t]*$')
            m = pattern.match(line)
            if m:
                release = m.group(1)+"_"+m.group(2)+"_"+m.group(3)
                releaseRevision[release] = m.group(4)
                releases.append(release)
                
            patternTest = re.compile(r'^[\s\t]*test\-(\d+)_(\d+)_(\d+)[\s\t]*:[\s\t]*([\d\.]+)[\s\t]*$')
            m = patternTest.match(line)
            if m:
                test = m.group(1)+"_"+m.group(2)+"_"+m.group(3)
                testRevision[test] = m.group(4)
                tests.append(test)

        return [ releases, tests ]
        


if __name__ == "__main__":
    print("test program")
    testTrigger = MadTestTrigger()
    testTrigger.setDebugMode(True)
    testTrigger.run("./") # local directory
    what = testTrigger.getWhatToDo()
    release = testTrigger.getRelease()
    print("for release "+release+", "+what)
