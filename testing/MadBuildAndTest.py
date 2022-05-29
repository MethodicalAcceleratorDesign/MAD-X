#!/usr/bin/python
from Notify import notify
from MadTestExceptions import MadException
import MadTrigTest

import threading
import time, datetime
import os
import sys
import re
import traceback
import optparse

global debugMode

threadList = []

rootDir = "/afs/cern.ch/user/n/nougaret/scratch1/mad-automation"

class MadBuildAndTestException(Exception):
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class updateTokensThread(threading.Thread):
    def __init__(self):
        threading.Thread.__init__(self)
        self.killed = False
        threadList.append(self)
        self.fileTokens = open('MadBuildAndTest.trace','w')
    def setDebugMode(self,mode):
        self.debugMode = mode
    def run(self):
        global debugMode
        
        origin = datetime.datetime.now() # return string
        duration = 0 # init
        keytabDir = '/extra/home/nougaret/MadTestScriptAuth'
        keytabFile = keytabDir + '/' + 'mykeytab'
        outfile = rootDir+'/kinitReply.txt'

        initialCmd = '/usr/sue/bin/kinit -k -t ' + keytabFile + ' nougaret' \
                     + ' > '+outfile

        while True:
            now = datetime.datetime.now()            
            try:
            
                os.system('rm '+outfile)
                os.system(initialCmd) # forget kinit -R: only worked 5+1 days
                # by default unless kinit -r 10days could be used initially
                # which is not the case with an acron => instead rely on
                # keytable
                f = open(outfile,'r')
                msg = ''
                lines = f.readlines()
                for line in lines:
                    msg = msg + line + '\n'
                f.close()
                os.system('rm '+outfile)

                # must provide path to cope with reduced acron environment
                os.system('/usr/kerberos/bin/klist > ' + outfile)             
                f = open(outfile,'r')             
                lines = f.readlines()               
                msg2 = ''
                for line in lines:
                    msg2 = msg2 + line + '\n'
                f.close()              
                os.system('rm '+outfile)

                # refresh AFS ticket
                # must provide path to cope with reduced acron environment
                os.system('/usr/bin/aklog')

                if debugMode:
                    notify('jean-luc','Kerberos / AFS ticket','ran for ' \
                           +str((now-origin).days)+' days\n\n'\
                           + 'kinit reply:\n\n'+ msg+'\n\nklist output:\n\n'+\
                           msg2)
                
            except:
                traceback.print_tb(self.fileTokens)
                notify('jean-luc','Automatic message (PROBLEM)',\
                       'failed to refresh tokens!')
            for sec in range(1,21600+1): # sleep 6 hours
            #for sec in range(1,7200+1): # sleep 2 hours               
                time.sleep(1) # 1 sec
                if self.killed == True:
                    return
            
    def kill(self):
        self.killed = True # work-around: not yet able to kill Python threads
        close(self.fileTokens)

os.chdir(sys.path[0]) # change to Python script's directory
# ???
currentDir= os.getcwd()

# temporary: for the Intel compiler, lines split after the 80th character
# unless we set up a specific environment variable to set the line length
# note that this variable will be transmitted to the child processes
# launched through os.sytem() calls.
# all 'madx < input > output' calls where madx is compiled with the ifort
# Intel compiler will take FORT_FMT_RECL into account.
os.putenv('FORT_FMT_RECL','256') # lines split after 256 chars instead of
# the default 80

reportFile = open(rootDir+"/MadBuildAndTest_Report.txt","w")
now = time.ctime()
reportFile.write("MadBuildAndTest.py report from "+now+"\n")

# option parsing
try:
#if True:
    parser = optparse.OptionParser()
    parser.add_option("-a","--automatic", action="store_true", \
                      dest = "automatic", default=False, \
                      help="wake-up once a day and check for trigger tag in the CVS")
    parser.add_option("-m","--manual", \
                      action="store_true", dest="manual", \
                      help="start straight away, looking at the specified tag")
    parser.add_option("-t","--tag", \
                      action="store", type="string", dest="specified_ver", \
                      help="specify tag madX-XX_YY_ZZ")
    parser.add_option("-g","--debug", action="store_true", dest="debug",\
                      default=False,\
                      help="no commit and file written. no message sent.")
    parser.add_option("-s","--silent", action="store_true", dest="silent",\
                      default=False,\
                      help="no e-mail sent to module keepers")

    # 4 december 2009
    parser.add_option("-d","--dev",action="store_true",help="build and test also with the 'Makefile_develop' Lahey Fortran compiler")
    parser.add_option("-n","--nag",action="store_true",help="build and test also with the 'Makefile_nag' NAG f95 compiler")
    
    (options,args) = parser.parse_args()

    additionalCompilers = ''
    if options.nag:
        additionalCompilers += ' --nag '
    if options.dev:
        additionalCompilers += ' --dev '

    # check if all the command line arguments could be processed as options

    if len(args)>0:
        raise MadBuildAndTestException("arguments not understood as valid options")

    # check it all selected options are compatible
    if options.automatic and options.manual:
        raise MadBuildAndTestException("options -a and -m are mutually exclusive")

    global debugMode
    debugMode = options.debug

    if options.manual:
        # is it a well formed release tag?
        runTest = 'run-test'
        pattern = r'^madX-\d_\d{2}_\d{2}$'
        m = re.match(pattern,options.specified_ver)
        if m:
            # extract the version from the CVS at the specified tag in a local mirror under ./MadCvsExtract
            try:
                # delete any pre-existing directory
                #os.chdir(rootDir)
                # 23 june 2009: all the following can be removed as it's MadBuildPy.pl which creates ./MadCvsExtract
                #os.system('rm -rf ./temporary_MadCvsExtract')
                #os.mkdir('./temporary_MadCvsExtract')
                #os.chdir('./temporary_MadCvsExtract')
                #command = 'cvs -d :gserver:isscvs.cern.ch/local/reps/madx checkout -r ' + options.specified_ver + ' madX'
                #notify('jean-luc','CVS extract command','under "'+os.getcwd()+'", issue: '+command)
                #os.system(command)
                
                release = options.specified_ver

                # did the checkout succeed?
                # 23 june 2009: actually this CVS extract is useless, as MadBuilPy.pl carries-out its on extract
                # => such check-out could be removed eventually. but for test purpose, it still should work !!!

                # must go back to rootDir
                #os.chdir(rootDir)
                #os.system('rm -rf ./MadCvsExtract') # clean-up. directory will be recreated by MadBuildPy.pl
                
            except:
                print("failed to extract the specified release from the CVS ("+options.specified_ver+")")
                notify('jean-luc','failed to extract '+options.specified_ver,'failed to extract specified release from CVS')
            pass
        else:
            raise MadBuildAndTestException("incorrect version specified")
    else:
        # retreive the release number from the CVS
        # (idem as MadTrigTest.pl)
        trigTest = MadTrigTest.MadTestTrigger()
        trigTest.setDebugMode(options.debug) # if false, won't tag the CVS with 'test-XX_YY_ZZ'
        trigTest.run(rootDir)
        runTest = trigTest.getWhatToDo()
        release = trigTest.getRelease() # the tag with which to extract the CVS
        os.chdir(rootDir) # changed in between
     
#else:    
except:
    print("failed to parse the command line. exit.")
    (type,value,tb) = sys.exc_info()
    traceback.print_tb(tb)
    sys.exit()

th = updateTokensThread()
if options.debug:
    th.setDebugMode(True)
else:
    th.setDebugMode(False)
th.start()

try:
    if runTest == 'do-nothing':
        reportFile.write("No new release detected => no need to run "+\
                         "the test-suite\n")
        global debugMode
        if debugMode:
            notify('jean-luc','no new release detected','last release = '+release)
        
        #reportFile.close()
    elif not runTest == 'run-test':
        reportFile.write("don't know what to do!\n")
        #reportFile.close()
    else:
        reportFile.write("now about to release "+release+"\n")
        if options.manual:
            notify('jean-luc','manual release','new release = '+release)
        else:
            notify('admin','new release detected','new release = '+release)       

        # leave the try statement (go to finally if any)
        #th.kill()
        #sys.exit()

        # no work report for the time being
        #os.chdir(rootDir)
        #os.system("./MadWorkReport.pl")

        reportFile.write("entering MadBuild.pl, with release-tag "+\
                     release+"\n")
        os.chdir(rootDir)
        os.system("rm -rf ./MadSvnExtract")
        try:
            os.system("rm -rf ./tempfile")
        except:
            pass
        # 4 decembre 2009: invoke MadBuild.py instead of MadBuildPy.pl
        # os.system("./MadBuildPy.pl "+ release + " > ./tempfile") # MadBuildPy.pl instead of MadBuild.pl
        try:
            os.system("./MadBuild.py "+additionalCompilers + " --release " + release + "> ./tempfile")
            notify('jean-luc','here','managed to launch '+"./MadBuild.py "+additionalCompilers+  " --release " + release + "> ./tempfile")           
        except:
            notify('jean-luc','here','failed to launch '+"./MadBuild.py "+additionalCompilers+ release + "> ./tempfile")
            raise("forward exception")
        tempfile = open('./tempfile','r')
        templines = tempfile.readlines()
        if templines[0]=='False': # the return status indicating a compilation failure
            if options.manual:
                notify('jean-luc','stop test','because compilation did not succeed.')
            else:
                notify('admin','stop test','because compilation did not succeed.')
            # shall we kill the thread as well?
            th.kill()
            sys.exit()
        
        reportFile.write("MadBuildPy.pl completed\n")

        # no MadTest for the time-being
        reportFile.write("entering MadTest.pl\n")
        os.chdir(rootDir)

        # for the time being, invoke debug mode for target cororbit only
        #os.system("./MadTest.pl ./MadCvsExtract/madX")
        # WARNING: MadTestPy.pl should replace MadTest.pl to avoid sending
        # e-mail - instead, notification should take place within Python.
        if options.debug: # or global debugMode
            notify('jean-luc','Start test','invoke ./MadTest.py '+additionalCompilers)
        # os.system("./MadTestPy.pl ./MadCvsExtract/madX debug=c6t")
        if options.silent:
            # December 2nd, 2009: invoke the latest MadTest.py
            #os.system("./MadTestPy.pl ./MadCvsExtract/madX silent") # no e-mail
            os.system("./MadTest.py "+additionalCompilers+" --silent") # assumes the path to the MadCvsExtract/madX is known ...
            # currently the '--silent' option is recognized but not taken into account by MadTest.py
        else:
            # December 2nd, 2009: invoke the latest MadTest.py
            #os.system("./MadTestPy.pl ./MadCvsExtract/madX")
            os.system("./MadTest.py "+additionalCompilers) # assumes the path  to the MadCvsExtract/madX is known by MadTest.py
        reportFile.write("MadTestPy.py completed\n")

        if options.manual:
            notify('jean-luc','Completed test','MadTest.py completed')
        else:
            notify('admin','Completed test','MadTest.py completed')
            
        # final debug test to check we still have access to AFS
        # make sure still the script still has access to AFS
        # by listing the examples directory
        directory =  '/afs/cern.ch/user/n/nougaret/scratch1/'+\
                    'mad-automation/TESTING/'
        command = 'ls '+directory+ ' > '+directory+'theFile'
        # apparently os.system can't handle relative path names!!!!!!!????
        #print command
        os.system(command)
        time.sleep(5)
        afsDebugFile = open(directory+'theFile','r')
        lines = afsDebugFile.readlines()
        msg =''
        for line in lines:
            msg = msg + line + '\n'
        if debugMode:
            notify('jean-luc','Very Last Message',\
                   'Contents of '+directory+'\n'+\
               msg)              
        os.system('rm '+directory+'theFile')
        # end of final test


        #   reportFile.close()

except:
    reportFile.write("exception caught when trying Python scripts in:\n"+\
                     currentDir+"\n")
    traceback.print_exc(file=reportFile)

    notify('jean-luc','Failure',\
           "exception caught when trying the successive Python scripts in:\n"+\
           currentDir)

th.kill()
reportFile.close()
