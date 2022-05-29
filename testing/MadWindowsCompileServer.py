#!/usr/bin/python

import socket

import shutil
import re
import os
import sys

HOST = ''                 # Symbolic name meaning all available interfaces
PORT = 7070              # Arbitrary non-privileged port
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.bind((HOST, PORT))
s.listen(1)
while True:
    os.chdir('C:\\') # return to C: root directory
    print('Now accepting connections ...')
    conn, addr = s.accept()
    print 'Connected by', addr
    if True: # do the following only once ... => could remove the test altogether 
        data = conn.recv(1024)
        if not data: break
        
        pattern = re.compile(r'([\w\d\-\_]+) asks: Compile MAD for Windows!')
        m = pattern.match(data)
        if m:
            # copy contents of /user/nougaret/MAD-X-WINDOWS/madX into C:/madXCompilationSandbox, after cleaning the later
            currentDir = os.getcwd()
            print("now in "+currentDir)
            try:
                shutil.rmtree('C:\\madXCompilationSandbox')
            except:
                print("no compilation sandbox had to be deleted")
            shutil.copytree('Y:\\MAD-X-WINDOWS\\madX','C:\\madXCompilationSandbox')
        
            # invoke the compilation
            print("remote invocation of MAD-X compilation on Windows")
            os.chdir('C:\madXCompilationSandbox')
            currentDir = os.getcwd()
            print("now in "+currentDir)
            os.system('MakefileIntel.bat')
            # copy executable to NFS for subsequent transfer to AFS web folder
            shutil.copy('C:\\madXCompilationSandbox\madx.exe',\
                        'Y:\\MAD-X-WINDOWS\\madX\madx.exe')
            clientHost = m.group(1)
            # send acknowldegement to the client
            # jluc rajoute a l-instant the following
            conn.send('Compilation completed\n')
            conn.close()
            
#            clientPort = 7071 # agreed-upon with the server
#            sAcknowledge = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
#            sAcknowledge.connect((clientHost, clientPort))
#            sAcknowledge.send('Compilation completed\n')
            # data = sAcknowledge.recv(1024)
#            sAcknowledge.close()
#            print 'Received', repr(data)
            # now should leave this inner loop and wait for a new connection ...
        else:
            print('unexpected data packet received.')
            
    # conn.send(data)
    
        conn.close()

