#!/usr/bin/python
import os

tag = "madX_4-01-23"

currentdir = os.getcwd();
os.chdir('/user/nougaret/MAD-X-WINDOWS')
shutil.rmtree('madX',ignore_errors=True);
cmd = 'svn co svn+ssh://svn.cern.ch/reps/madx/tags/' + tag
os.system(cmd)
os.chdir(currentDir) # back to initial location
