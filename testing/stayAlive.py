#!/usr/bin/python

import threading
from Notify import notify
import time
import os
import sys
import re

class updateTokensThread(threading.Thread):
    def __init__(self):
        threading.Thread.__init__(self)
    def run(self):
        startTime = time.ctime()
        while True:
            try:
                # provide path to cope with reduced acron environment
                os.system('/usr/sue/bin/kinit -k -t /extra/home/nougaret/MadTestScriptAuth'+
                          '/mykeytab nougaret')
                os.system('/usr/bin/aklog')
                # must provide path to cope with reduced acron environment
                os.system('/usr/sue/bin/klist > ' + outfile)             
                f = open(outfile,'r')             
                lines = f.readlines()               
                msg2 = ''
                for line in lines:
                    msg2 = msg2 + line + '\n'
                f.close()              
                os.system('rm '+outfile)               
            except:
                notify('jean-luc','failure','failed to refresh tokens!')
                for sec in range(1,30): # sleep 30 seconds
                    time.sleep(1) # 1 sec
                    if self.kill == True:
                        return
            now = time.ctime()
            notify('jean-luc','alive','still alive after '+ +'hours')
    def kill(self):
        self.kill = True # work-around: not yet able to kill Python threads

print 'now about to launch thread'
th = updateTokensThread()
th.start()
print 'threaded launched'

while True:
    pass

th.kill()
