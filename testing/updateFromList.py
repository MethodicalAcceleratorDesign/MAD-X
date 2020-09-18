from os import listdir
import os
from os.path import isfile, join
from shutil import copyfile
'''
 + test-ptc-twiss-1                                   (0.00 s) -  2/ 4 : FAIL
 + test-ptc-twiss-2                                   (0.21 s) -  1/ 3 : FAIL
 + test-ptc-twiss-3                                   (0.00 s) -  2/ 4 : FAIL
 + test-ptc-twiss-4                                   (0.00 s) -  0/ 1 : FAIL
 + test-ptc-twiss-old1                                (0.00 s) -  3/ 6 : FAIL
 + test-ptc-twiss-old2                                (0.00 s) -  3/ 6 : FAIL
 + test-ptc-twiss-old3                                (0.00 s) -  3/ 7 : FAIL
 + test-ptc-twiss-old4                                (0.00 s) -  3/ 6 : FAIL
 + test-ptc-twiss-old5                                (0.03 s) -  3/ 5 : FAIL
 + test-ptc-twiss-old6                                (0.11 s) -  2/ 5 : FAIL
 + test-ptc-twiss-old7                                (0.01 s) -  3/ 8 : FAIL
 + test-ptc-twiss-5Dt                                 (0.01 s) -  2/ 5 : FAIL
 + test-ptc-twiss-56Dt                                (0.00 s) -  0/ 3 : FAIL
 + test-ptc-twiss-56Dtl                               (0.00 s) -  0/ 3 : FAIL
 + test-ptc-twiss-6D                                  (0.00 s) -  0/ 3 : FAIL
 + test-ptc-twiss-6D-ALS                              (0.01 s) -  0/ 3 : FAIL
 + test-ptc-twiss-56Dt-ini_map_man                    (0.02 s) -  0/ 5 : FAIL
 + test-ptc-twiss-56Dt-ini_mtx_man                    (0.00 s) -  0/ 3 : FAIL
 + test-ptc-twiss-56Dt-ini_mtx_tbl                    (0.01 s) -  1/ 5 : FAIL
 + test-ptc-normal                                    (0.15 s) -  4/ 5 : FAIL
 + test-ptc-twiss-normal-5D                           (0.01 s) -  4/ 7 : FAIL
 + test-ptc-twiss-normal-6D                           (0.05 s) -  4/ 7 : FAIL
 + test-ptc-track-5                                   (0.03 s) -  4/ 6 : FAIL
 + test-ptc-track-6Dtl-acd                            (0.04 s) -  3/ 5 : FAIL
'''

#copyfile(src, dst)
mypath_orig = '../tests/'
file1 = open('../myfail.txt', 'r') 
Lines = file1.readlines() 
loc = os.getcwd()
count = 0
# Strips the newline character 
for line in Lines: 
	os.chdir(loc)
	mypath = mypath_orig + line.strip() + "/"
	print(mypath)
	#mypath = mypath +"test-emit/"
	print(mypath)
	os.chdir(mypath)
	onlyfiles = [f for f in listdir(".") if isfile(join(".", f))]
	#for i in range(0, len(onlyfiles)):
	#	onlyfiles[i] = [mypath + onlyfiles[i]]

	for f in onlyfiles:
		if(f.endswith('.ref')):
			oname = f[:-4]
			print(oname, f)
			if(isfile(oname)):
				copyfile(oname, f)
		if(f.endswith(".out")):
			cpname = f[:-4] + ".ref"
			copyfile(f, cpname)

	#print(onlyfiles)
