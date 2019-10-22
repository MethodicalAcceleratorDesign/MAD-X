import numpy as np
import os,shutil
import matplotlib.pyplot as plt

import sys
import subprocess


sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")
try:
     from metaclass import *
except:
     from metaclass25 import *

ang = np.linspace(0,2*np.pi, 20)
r = 0.01
ra = 0.01
x = r*np.cos(ang)
y = r*np.sin(ang)
x[8]=0
length = 20
f = open("myge.dat", "w")
for i in range(0,len(x)):
	f.write(str(x[i]) + " "+ str(y[i]) +"\n")
f.close()
f = open("internal.dat", "w")
for i in range(0,len(x)):
	f.write(str(x[i]) +" ,")

f.write("\n")
for i in range(0,len(x)):
	f.write(str(y[i]) +" ,")
f.close()



xo=[]
yo=[]
xi=[]
yi=[]



xp = np.linspace(-ra,ra, length)
yp = np.linspace(-ra,ra, length)
#x1 = [0.011, 0.0075, 0.0025,-0.0025]
#y1 = [0.00, 0.0075, 0.0025,0.0025]
with open('mask.madx', 'r') as file :
  filedatain = file.read()

for i in range(0,len(xp)):
	for j in range(0,len(yp)):
		filedata = filedatain.replace('xini', str(xp[i]))
		filedata = filedata.replace('yini', str(yp[j]))
		print("qqqqqqqqqq", xp[i], yp[j])
		with open('torun.madx', 'w') as file:
			file.write(filedata)
		#os.system('/home/tobias/codes/MAD-X/madx64 torun.madx')
		subprocess.check_call('/home/tobias/codes/MAD-X/madx64 torun.madx', shell=True)
		track = twiss('track.outone')
		#p = subprocess.Popen('/home/tobias/codes/MAD-X/madx64 ','torun.madx')
  		#p.wait()
		print("aaaaaaaa",len(track.X) )
		if(len(track.X) == 1):
			xo.append(xp[i])
			yo.append(yp[j])
		else:
			xi.append(xp[i])
			yi.append(yp[j])



print(xo)





plt.plot(x,y)
plt.plot(xo,yo, 'r*')
plt.plot(xi,yi, 'go')
plt.show()


#u = np.loadtxt("myge.dat")
#plt.plot(u[0:2:-1], u[1:2:-1])
#plt.show()
print(u)
