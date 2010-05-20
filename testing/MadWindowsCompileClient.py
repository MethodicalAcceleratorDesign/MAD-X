#!/usr/bin/python

import re
import os
import sys
import traceback
import shutil
import time

import socket

# public:

def releaseForWindows(tag): # to be called by MadTrigRelease.pl
    try:
        client = windowsReleaseClient()
        client.release(tag)
        return 'success'
    except:
        exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        traceback.print_tb(exceptionTraceback)
        return 'failure'

# private:

def getMonthAsString(monthInt):
    monthStrs = ['Jan','Fev','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    monthStr = monthStrs[monthInt-1]     
    return monthStr

class windowsReleaseClient:

    def __init__(self):
        self.megabytes = 0 # the size of the binary
        self.version = 'undefined' # defined in release
        self.startTimeStr = '' # defined in remoteCompile()
        self.endTimeStr = ''# defined in remoteCompile()
        pass
    
    def extractSVN(self,tag): # populate /user/nougaret/MAD-X-WINDOWS

        currentDir = os.getcwd()

        os.chdir('/user/nougaret/MAD-X-WINDOWS')

        shutil.rmtree('madX',ignore_errors=True) # clean-up
        #cmd = 'cvs -d :gserver:isscvs.cern.ch:/local/reps/madx '+\
        #          'checkout -r '+tag+' '+\
        #          'madX'
        cmd = 'svn co svn+ssh://svn.cern.ch/reps/madx/tags/' +\
              tag +'/madX'+ ' madX'        
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
        print('send compilation request to the Windows host...')
        startTime = time.localtime()
        self.startTimeStr = str(startTime.tm_mday)+' '+getMonthAsString(startTime.tm_mon)+\
                      ' '+str(startTime.tm_year)+ ' '+str(startTime.tm_hour)+':'+str(startTime.tm_min)
        data = s.recv(1024) # blocking on reception (what about the time out???)
        s.close()
        endTime = time.localtime()
        self.endTimeStr = str(endTime.tm_mday)+' '+getMonthAsString(endTime.tm_mon)+\
                      ' '+str(endTime.tm_year)+ ' '+str(endTime.tm_hour)+':'+str(endTime.tm_min)        
        print 'Received', repr(data)

    def checkCompileOutcome(self):
        # look into MAD-X-WINDOWS and check whether a fresh madx.exe is present
        try:
            os.system('rm ./outfile') # if any
        except:
            pass
        os.system('ls -l /user/nougaret/MAD-X-WINDOWS/madX/madx.exe > ./outfile')
        f = open('./outfile','r')
        lines = f.readlines()
        f.close()
        os.system('rm ./outfile')
        if len(lines) != 1:
            raise('error trying to access ./outfile for checking madx.exe')
        else:
            singleLine = lines[0]
            pattern = re.compile(r'(\d+)[\s\t]+([\w\W]{3})[\s\t]+(\d{1,2})[\s\t]+(\d{2}):(\d{2})')
            m = pattern.search(singleLine)
            if m:
                size = float(m.group(1))
                month = m.group(2)
                day = m.group(3)
                hour = m.group(4) # can be year if from last year, in which case 200X instead of xx:yy
                min = m.group(5)
                self.megabytes = size/1000000.0 # for subsequent use when generating HTML page
                print("size="+str(self.megabytes)+", month="+month+", day="+day+",hour="+hour+", min="+min)
                # compare with local time - should not be older than 5 minutes
                now = time.localtime()
                #monthStrs = ['Jan','Fev','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
                #monthStr = monthStrs[now.tm_mon-1]
                monthStr = getMonthAsString(now.tm_mon)
                print("now.tm_mon="+monthStr+", now.tm_day="+str(now.tm_mday)+\
                      ",now.tm_hour="+str(now.tm_hour)+", min="+str(now.tm_min))  
                print('monthStr='+monthStr+', month='+month+', now.tm_mday='+str(now.tm_mday)+', day='+day)
                if (monthStr==month) and (now.tm_mday==int(day)):
                    # check if generated less than 5 minutes ago
                    nowTotalMinutes = 60*(now.tm_hour) + now.tm_min # minutes elapsed since now
                    totalMinutes = 60*int(hour) + int(min) # minutes elapsed since madx.exe created
                    print("nowTotalMinutes="+str(nowTotalMinutes)+", totalMinutes="+str(totalMinutes))
                    print("differences in minutes="+str(nowTotalMinutes-totalMinutes))
                    if (nowTotalMinutes - totalMinutes)>5:
                        raise('madx.exe is older than 5 minutes old. abort.')
                    else:
                        # done
                        pass
                else:
                    raise('madx.exe is too old. abort.')
                print("now year="+str(now.tm_year)+", month="+str(now.tm_mon)+\
                      ", day="+str(now.tm_mday)+", hour="+str(now.tm_hour)+", min="+str(now.tm_min))
            else:
                raise('fail to match the date')
        pass
        print("completed checkCompileOutcome")
       # could also gather the executable size here
    

    def generateHtmlOutput(self):
        # move MAD-X fresh executable to the AFS web folder
        executablesAfsWebFolder = '/afs/cern.ch/user/n/nougaret/www/mad/windows-binaries'        
        # (it has already been put into MAD-X-WINDOWS/madX by the Windows host)
        os.system('cp /user/nougaret/MAD-X-WINDOWS/madX/madx.exe '+\
                 executablesAfsWebFolder + '/madx.exe' ) # hard-coded!
        # now produce the HTML page (may be static actually)
        htmlFile = executablesAfsWebFolder + '/executables.htm'
        contents =''
        contents += "<p>Version "+ self.version+" compiled with Intel ifort and Microsoft Visual C++:</p>\n";
        contents += "<table width=\"75%\" border=\"0\">\n";
        oddOrEven = 'even' # to colorize successive lines differently

	contents += "<tr class=\"odd\"><td>Download</td><td><a href=\"./madx.exe\">madx.exe</a></td><td>("+\
                    str(self.megabytes)+" Megabytes)</td><td>for the latest version.</td></tr>\n";
    
        contents += "</table>\n"
        contents += "<p>Version 3.04.53 compiled with Lahey Fortran and accepting sequences with BV flag, as until March 2009:</p>\n"
        contents += "<table width=\"75%\" border=\"0\">\n"
        contents += "<tr class=\"even\"><td>Download</td><td><a href=\"./madx-old.exe\">madx-old.exe</a></td><td>(2.6132 Megabytes)</td><td>for the archived version, without PTC.</td></tr>\n"
        contents += "<tr class=\"odd\"><td>Download</td><td><a href=\"./madxp-old.exe\">madxp-old.exe</a></td><td>(6.7554 Megabytes)</td><td>for the archived version, including PTC.</td></tr>\n"
        contents += "</table>\n"

        # create web page in the correct AFS web folder location
        html = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">'
        html += "<html>\n"
        html += "<head>\n"
        html += "<title>MAD-X downloadable executables</title>\n"
        html += "<link rel=stylesheet href=\"../MadTestWebStyle.css\" type=\"text/css\">" # CSS one level up
        html += "</head>\n"
        html += "<!-- generated by Windows compilation script -->\n"
        html += "<body>\n"
        html += "<p>Windows compilation started "+self.startTimeStr+", ended "+self.endTimeStr+"</p>\n"
        html += contents
        html += "</body>\n"
        html += "</html>\n"

        htmlFile = executablesAfsWebFolder + '/executables.htm'
        f = open(htmlFile,'w')
        f.write(html)
        f.close()
    

    def release(self,tag):
        pattern = re.compile(r'^madX\-(\d+)_(\d+)_(\d+).*$') # pro, dev ok
        m = pattern.match(tag)
        if m:
            major = m.group(1)
            medium = m.group(2)
            minor = m.group(3)
            self.version = major + '.' + medium + '.' + minor
        else:
            raise('ill-formed release tag')
        # first extract the SVN on NFS
        self.extractSVN(tag)
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
