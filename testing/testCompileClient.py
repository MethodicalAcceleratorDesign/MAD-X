#!/usr/bin/python
import os
import shutil

tag = "madX-4_01_23"

currentDir = os.getcwd();
os.chdir('/user/nougaret/MAD-X-WINDOWS')
shutil.rmtree('madX',ignore_errors=True);
cmd = 'svn co svn+ssh://svn.cern.ch/reps/madx/tags/' + tag +'/madX'+ ' madX'
os.system(cmd)
os.chdir(currentDir) # back to initial location
