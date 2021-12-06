from os import listdir
import os
from os import listdir
from os.path import isfile, join

from os.path import isfile, join
from shutil import copyfile, move

def my_listdirs(folder):
    return [d for d in os.listdir(folder) if os.path.isdir(os.path.join(folder, d))]

def change_out(mypath, f_name):
	old_f = mypath + "/" + f_name + '.out'
	new_f = mypath + "/" + f_name + '.ref'
	if(os.path.exists(old_f) and os.path.exists(new_f)):
		copyfile(old_f, new_f)
def update_cfg_twiss(out_file):
	cfg_file = out_file + ".cfg"
	cfg_file_tmp = out_file + ".cfg.tmp"
	f_cfg = open(cfg_file, "r")
	f_cfg_tmp = open(cfg_file_tmp, "w")
	lines = f_cfg.readlines()
	for line in lines:
		isskip = "skip" in line
		isnumb = "43" in line or "44" in line or "47" in line or "48" in line
		if(not isskip and not isnumb):
			f_cfg_tmp.write(line)
	f_cfg_tmp.write("48-50     *     skip    # twissHead")
	f_cfg_tmp.close()
	move(cfg_file_tmp, cfg_file)
def write_new_columns(ref_file):
	#print(ref_file)
	out_file = ref_file[0:-4]
	tmp_file = ref_file[0:-4] + ".tmp"
	if(os.path.exists(out_file)):
		f_out = open(out_file, "r")
	else:
		f_out = open(out_file + ".out", "r")
	f_ref = open(ref_file, "r")
	f_tmp = open(tmp_file, "w")

	out_lines = f_out.readlines()
	lines = f_ref.readlines()
	count = 0
	isUpdated = "SYNCH_6" in lines[42]
	print (lines[42], isUpdated) 	
	for line in lines:
		count = count + 1
		if(count == 43 and not isUpdated):
			f_tmp.write(out_lines[42])
			f_tmp.write(out_lines[43])
			f_tmp.write(out_lines[44])
			f_tmp.write(out_lines[45])
		f_tmp.write(line)
	f_out.close()
	f_ref.close()
	f_tmp.close()
	move(tmp_file, ref_file)
	update_cfg_twiss(out_file)
def find_twiss_files(mypath):
	onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
	for mf in onlyfiles:
		mf = mypath + "/"  + mf 
		if(".ref" in mf):
			f = open(mf, "r")
			lines = f.readlines()
			count = 0
			for line in lines:
				count = count + 1			
				if ("SYNCH_1" in line and not "SYNCH_2" in line):
					#addcolumns(
					print("TWISS file is:", mf)
					write_new_columns(mf)
					break

	

#copyfile(src, dst)
mypath = '../tests/'	
allfolders = my_listdirs(mypath)
print(allfolders)
for folder in allfolders:
#for folder in ['test-makethin-3']:
	if("share" not in folder and "lhc" not in folder):
		mypath = '../tests/'		
		mypath = mypath + folder
		print("myyyy", mypath)	
		change_out(mypath, folder)
		find_twiss_files(mypath)
	
		

