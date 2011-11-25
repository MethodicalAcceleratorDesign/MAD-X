import os,sys
import optparse
import numpy
import platform

class InputError(Exception):
    pass

def tfs(inputfile):
    '''
    Give file path as input, returns tfs table as output
    '''
    summary={}
    PYVER=int(platform.python_version().split('.')[0])
    if not os.path.isfile(inputfile):
        if os.path.isfile(inputfile+'.tfs'):
            inputfile+='.tfs'
        elif os.path.isfile(inputfile+'.TFS'):
            inputfile+='.TFS'
        else:
            raise ValueError("ERROR: "+inputfile+" is not a valid file path")
    if PYVER==3:
        f=open(inputfile,'r')
    else:
        f=file(inputfile,'r')
    l=f.readline()
    while(l):
        if l.strip()[0]=='@':
            _addParameter(summary,l)
        if l.strip()[0]=='*': # beginning of vector list...
            names=l.split()[1:]
            table=_read_table(f,names)
        l=f.readline()
    return table, summary
 
# very simple requirement at the moment..
def compare(arr1,arr2):
    return numpy.average(arr1-arr2)< 1.0/len(arr1)
    
    
# Give two tables as input, returns 0 if they are
# similar enough...
def check_tables(t1,t2, add_keys):
    needed_keys = ['s']
    needed_keys.extend(add_keys)
    for k in needed_keys:
        if k not in t1 or k not in t2:
            print("Could not find key %s" %k)
            return 1
    if (t1['s'] != t2['s']).any():
        raise ValueError("ERROR: position array (s-coordinates) changed")
    for k in add_keys:
        if not compare(t1[k],t2[k]):
            print("Too large difference for key "+k)
            return 5

def _addParameter(summary,line):
    '''
     helper function for tfs()
    '''
    lname=line.split()[1].lower()
    if line.split()[2]=='%le':
        summary[lname]=float(line.split()[3])
    if line.split()[2][-1]=='s':
        summary[lname]=line.split('"')[1]
    if line.split()[2]=='%d':
        summary[lname]=int(line.split()[3])


def _read_table(fstream,names):
    '''
     helper function for tfs()
    '''
    l=fstream.readline()
    types=[]
    table={}
    for n in names:
        table[n.lower()]=[]
    while(l):
        if l.strip()[0]=='$':
            types=l.split()[1:]
        else:
            for n,el in zip(names,l.split()):
                table[n.lower()].append(el)
        l=fstream.readline()
    for n,typ in zip(names,types):
        if typ=='%le':
            table[n.lower()]=numpy.array(table[n.lower()],dtype=float)
        elif typ=='%d':
            table[n.lower()]=numpy.array(table[n.lower()],dtype=int)
        elif typ=='%s':
            for k in range(len(table[n.lower()])):
                table[n.lower()][k]=table[n.lower()][k].split('"')[1]
    return table

# Checks that two tfs tables are sufficiently similar..
if __name__=="__main__":
    usage = '%prog olt.tfs new.tfs'
    parser = optparse.OptionParser(usage)
    (options,args) = parser.parse_args()
    if len(args)!=2:
        print(len(args))
        parser.print_help()
        print("\n      -- Wrong input -- \n")
        sys.exit(1)
    for f in args:
        if not os.path.isfile(f):
            parser.print_help()
            print("%s is not a valid file path" % f)
            sys.exit(1)
    print("Checking %s against %s" % (args[0],args[1]))
    ret=check_tables(tfs(args[0])[0],tfs(args[1])[0],['betx','bety'])
    if not ret:
        print("\n\t--- Result produced in this test seems to be of good quality ---\n")
    sys.exit(ret)
    
            
    
