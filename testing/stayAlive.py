#!/usr/bin/python
import threading
import os
import time

class updateTokensThread(threading.Thread):
    def __init__(self):
        threading.Thread.__init__(self)
    def run(self):
        while True:
            try:
                os.system('/usr/sue/bin/kinit -k -t /extra/home/nougaret/MadTestScriptAuth'+
                          '/mykeytab nougaret')
                os.system('/usr/bin/aklog')           
                outfile = "stayAlive.out"
                os.system('/usr/kerberos/bin/klist > ' + outfile)            
                f = open(outfile,'r')             
                lines = f.readlines()               
                msg2 = ''
                for line in lines:
                    msg2 = msg2 + line + '\n'
                f.close()              
                os.system('rm '+outfile)
                print('klist contents:' + msg2)
            except:
                print('failed to refresh tokens!')
            for sec in range(1,30): # sleep 30 seconds
                time.sleep(1) # 1 sec
                if self.kill == True:
                    return
    def kill(self):
        self.kill = True # work-around: not yet able to kill Python threads

th = updateTokensThread()
th.start()

while True:
    pass

th.kill()
