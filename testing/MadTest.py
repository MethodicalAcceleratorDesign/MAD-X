#!/usr/bin/python

# currently works with following invocation:
# MadTest.py -k -m Makefile -t aperture
# MadTest.py -k -m Makefile -t ptc_twiss -c example1.madx
# MadTest.py -k -m Makefile -t ptc_twiss
import optparse
import os
import shutil
import re
import Notify

class State: # the state of a test
    pass

init = State()
incomplete = State() # missing resources
ready = State() # all resources present
running = State()
aborted = State()
completed = State()

rootDir = "/afs/cern.ch/user/n/nougaret/scratch0/mad-automation/madX-examples/REF"

topDir = "/afs/cern.ch/user/n/nougaret/scratch0/mad-automation"
testingDir = topDir + "/TESTING"
repositoryDir = topDir + "/madX-examples"
madxDir = topDir + "/MadCvsExtract/madX"

htmlRootDir = "/afs/cern.ch/user/n/nougaret/www/mad"
mainHtmlPage = htmlRootDir + "/test_new.htm"

makefiles = ['Makefile','Makefile_develop','Makefile_nag']
    
# all files with absolute path

class Repository: # wrapper on CVS or SVN repository for examples
    def __init__(self):
        pass
    def checkout(self):
        # clean-up the directory in which the examples's repository is to be extracted
        try:
            shutil.rmtree(repositoryDir)
        except:
            pass # directory absent
        currentDir = os.getcwd()
        os.mkdir(repositoryDir)
        os.chdir(topDir)
        checkoutCommand = "cvs -d :gserver:isscvs.cern.ch:/local/reps/madx-examples checkout madX-examples";
        os.system(checkoutCommand)
        os.chdir(currentDir)

class Resource:
    def __init__(self,name,source,destinations): # the name is not unique, as this is within scope of a testcase
        # for a given resource, there are 3 destinations, i.e. one for each Makefile
        self.name = name
        self.source = source # the expanded filename of the file from CVS
        self.destinations = destinations # the expanded filename of the file in local work directory
        if options.verbose:
            print("create resource '"+name+"' retreived from '"+source+"'")

    def lookForDependantResources(test,filename,lvl=0): # returns a list of filenames with absolute path

        resources = [] # list to be returned
        
        if lvl==0:
            source = filename
            name = source[source.rfind('/')+1:]
            destinations = []
            for m in makefiles:
                destination = source.replace(rootDir,testingDir+'/'+m) # replace prefix between source and destination
                destinations.append(destination)            
            resource = Resource(name,source,destinations)  # primary resource is the MAD-X input file itself
            resources.append(resource)
            
        try:
            f = open(filename,'r')
            # regex to input a file (not anchored as an instruction may precede)
            commonPatterns = (re.compile(r'readmytable,?[\s\t]*file[\s\t]*=[\s\t]*[\"\']?([\w\.\_\d/]+)[\"\']?[\s\t]*;'),\
                              re.compile(r'call,?[\s\t]*file[\s\t]*=[\s\t]*[\"\']?([\w\.\_\d/]+)[\"\']?[\s\t]*;'),\
                              re.compile(r'readtable,?[\s\t]*file[\s\t]*=[\s\t]*[\"\']?([\w\.\_\d/]+)[\"\']?[\s\t]*;'))
            for line in f.readlines():
                line = line.rstrip('\n')
                if len(line.lstrip())>=1 and line.lstrip()[0]=='!': # if line is a comment discard it
                    continue
                k = line.find('=')
                if not k == -1:
                    # put to lower case everything in front of the equal sign to simply patterns
                    line = line[0:k].lower() +"="+ line[k+1:]
                for p in commonPatterns:
                    m = p.search(line)
                    if m:
                        source = m.group(1)
                        if not source[0] == '/': # relative path to be converted to absolute
                            dirPrefix = filename[0:filename.rfind('/')]                      
                            source = dirPrefix+'/'+source
                            source = os.path.normpath(source) # normalize the pathname
                        else:
                            pass

                        # the short name of the resource is obtained by removing the complete path prefix
                        name = source[source.rfind('/')+1:]
                        # for each source, there are three destinations as there are 3 different makefiles
                        destinations = []
                        for m in makefiles:
                            destination = source.replace(rootDir,testingDir+'/'+m) # replace prefix between source and destination
                            destinations.append(destination)
                          
                        resource = Resource(name,source,destinations)
                        resources.append(resource)
                        
                        # now dig further
                        additionalResources = Resource.lookForDependantResources(test,source,lvl+1)                     
                        for a in additionalResources:
                            resources.append(a)

            f.close()
        except: # failed to open the file: this may be a temporary file to be created at run time, any way one should
            # mark the resource as missing
            print("caught exception")
            test.missing_resources.append(filename)
            if test.state == init:
                test.state = incomplete;

        return resources
    lookForDependantResources = staticmethod(lookForDependantResources)

class Target: # a named "module" holding on or several tests
    targets = []
    targetsDict = {} # for fast access (but random order)
    def __init__(self,name):
        self.name = name
        self.tests = []
        Target.targets.append(self)
        Target.targetsDict[name] = self
        if options.verbose:
            print("create target "+name)
    def registerTest(targetname,test):
        registered = False # default
        for t in Target.targets:
            if t.name == targetname:
               registered = True
        if not registered:
                target = Target(targetname)
                target.tests.append(test)
        else:
            for t in Target.targets:
                if t.name == targetname:
                    t.tests.append(test)
    registerTest = staticmethod(registerTest)

class Test:
    tests = []

    def __init__(self,name,program,input,output,subdirectory):

        self.name = name # the target
        self.program = program
        self.input = input
        self.output = output
        self.resources = []
        self.missing_resources = []
        self.state = init # initial state
        
        Target.registerTest(name, self)

        if subdirectory == 0:
            self.testcaseDir = "test_"+str(len(Target.targetsDict[name].tests))
        else:
            self.testcaseDir = "test_"+str(len(Target.targetsDict[name].tests))+"_"+subdirectory
        
        if not subdirectory == 0:
            self.subdirectory = subdirectory
        else:
            self.subdirectory = ''
            
        Test.tests.append(self) # irrespective of whether the tests has access to all resources it requires

        if self.collectResources():
            if options.verbose:
                print("create test '"+name+"' with program '"+program+"', input '"+input+"', and output '"+output+"'")

        else:
            if options.verbose:
                print("fail to create test '"+name+"' with program = '"+program+"', input = '"+\
                      input+"', and output = '"+output+", subdir='",subdirectory)
            self.state = incomplete
            # issue a warning stating that the file has incomplete resources


    def collectResources(self):
        for entry in os.walk(rootDir,topdown=True): # walk file and directory structure under root

            foundDirectory = False

            if not self.subdirectory == '':
                if rootDir+'/'+self.name+'/'+self.subdirectory == entry[0]:
                    foundDirectory = True
            else:
                if rootDir+'/'+self.name == entry[0]:
                    foundDirectory = True
                                
            if foundDirectory:
                if self.input in entry[2]: # got the primary input madx file
                    # now check that all included files are also available (resources)
                    self.resources = Resource.lookForDependantResources(self,entry[0]+'/'+self.input) # provide expanded filename
                    return True # found the input madx file
        return False

    def copyResourcesToTestingDir(self):
        for r in self.resources:
            if options.verbose:
                print("now to copy '"+r.source+"'")
            for d in r.destinations: # a single resource has several destinations, i.e. one per Makefile
                destinationDir = d[:d.rfind('/')]
                if not os.path.exists(destinationDir): # output directory does not exist => create it before copying file
                    os.makedirs(destinationDir) # also create all intermediate-level directories to contain the leaf directory
                        
                shutil.copyfile(r.source,d)

    def run(self):
        currentDir = os.getcwd()
        for i,m in enumerate(makefiles):
            # retreive the name of the directory in which the madx input is located (first resource)            
            script = self.resources[0].destinations[i]
            scriptDir = script[:script.rfind('/')]
            os.chdir(scriptDir)
            command = madxDir+'/madx_'+m+' <'+self.input +'>'+self.output
            print("now to execute "+command+" under "+scriptDir)
            os.system(command)

            # now create two subdirectories, input and output to store the outcome
            os.mkdir('./input')
            os.mkdir('./output')
            class IOType:
                pass
            input = IOType()
            output = IOType()

            for f in os.listdir(scriptDir):
                if os.path.isdir(f):
                    continue
                type = output # default                
                for r in self.resources:
                    print("compare "+scriptDir+"/"+f+" with "+r.destinations[i])
                    if scriptDir+"/"+f == r.destinations[i]: # this is an (input) resource file
                        type = input
                if type == input:
                    shutil.move(f,'./input')
                else:
                    shutil.move(f,'./output')

            
            os.chdir(currentDir) # back to the initial work directory
            
class Tester:



    def __init__(self):
        pass

    def run(self):

        # first extract examples from the repository
        if not options.keep_data:
            rep = Repository()
            rep.checkout()

        testinfoPattern = re.compile(r'^[\s\t]*(.+)[\s\t]*<[\s\t]*(.+)[\s\t]*>[\s\t]*([\w\.\d]+)[\s\t]*,?[\s\t]*'+\
                                     '(subdirectory=)?[\s\t]*(.*)[\s\t]*$')
        # e.g. ./madx < touschek.lhcinjection.madx >  touschek.lhcinjection.out, subdirectory=LHC_injection
        # or ./madx < lep.madx  >  lep.out

        # process TestScenario.xml
        os.system('xsltproc --stringparam what list_targets ProcessScenario.xsl'+\
                  ' TestScenario.xml > ./outfile')
        f = open('./outfile','r')
        targets = f.readlines()
        f.close()
        os.remove('./outfile')
        for target in targets:
            target = target.rstrip('\n')
            
            if options.singleTarget and not options.singleTarget == target:
                continue # skip this test
            
            # extract detailed information about each test
            os.system('xsltproc --stringparam what list_tests'+\
                      ' --stringparam target '+target+\
                      ' ProcessScenario.xsl TestScenario.xml > ./outfile')
            f = open('./outfile','r')
            testinfos = f.readlines()
            f.close()
            os.remove('./outfile')
            for testinfo in testinfos:
                testinfo = testinfo.rstrip('\n')
                
                m = testinfoPattern.match(testinfo)
                if m:
                    program = m.group(1).rstrip()
                    input = m.group(2).rstrip()
                    output = m.group(3).rstrip()
                    if m.lastindex == 5:
                        subdirectory = m.group(5)
                    else:
                        subdirectory = 0

                    if options.singleCase and not options.singleCase  == input:
                        continue
                
                    test = Test(target,program,input,output,subdirectory)

                else:
                    if options.verbose:
                        print("failed to parse line "+testinfo)

            # cleanup the TESTING directory structure
            try:
                shutil.rmtree(testingDir)
            except:
                pass # directory absent

            # now populate testingDir with all sources and associated resources
            for t in Test.tests:
                t.copyResourcesToTestingDir()

            # now run the tests
            for t in Test.tests:
                t.run()

            # notify module keepers if required
            if options.notify:
                Notify.notify("jean-luc","test completion","test completed.") # for the time-being

        page = WebPage(mainHtmlPage)
        page.output()
        if not options.quiet:
            page.display()

class WebPage:
    
    def __init__(self,name):
        self.name = name
        self.contents = ''

    def header(self):
        self.contents += '<DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">\n'
        self.contents += '<html>\n'
        self.contents += '<head>\n'
        self.contents += '<title>MAD testing main page</title>'
        self.contents += '<link rel=stylesheet href="./MadTestWebStyle.css" type="text/css">'
        self.contents += '</head>'
        
    def body(self):
        self.contents += '<body>\n'
        self.contents += '<table>\n'
        for target in Target.targets:
            self.contents += '<tr class="test_target"><td colspan="2"><div align="center"><strong>'+target.name+'</strong></div></td></tr>'  
            for i,t in enumerate(target.tests):
                self.contents += '<tr class="test_case"><td width=\"80%\">'+t.testcaseDir+\
                                 ': '+t.program +'&lt;'+t.input+'&gt;'+t.output+\
                                 '</td><td width=\"20%\"><table width=\"100%\" style=\"text-align: center\"><tr>'+\
                                 'DEVRES'+'NAGRES'+\
                                 '</tr></table></td></tr>\n'; 
                for r in t.resources:
                    self.contents += '<tr><td>'+r.name+'</td></tr>'
        self.contents += '<table>\n'
        self.contents += '<body>\n'
    def footer(self):
        self.contents += '</html>\n'

    def output(self):
        self.header()
        self.body()
        self.footer()
        f = open(self.name,'w')
        for c in self.contents:
            f.write(c)
        f.close()
        
    def display(self):
        os.system('firefox ' + self.name+ '&')

if __name__ == "__main__":
    usage = "%prog [options]"
    parser = optparse.OptionParser(usage)
    parser.add_option("--verbose","-v",help="display messages useful for debug", action="store_true")
    parser.add_option("--target","-t",help="select a specific target to run the test",dest="singleTarget")
    parser.add_option("--case","-c",help="select a specific test-case to run the test, in combination with --target",\
                      dest="singleCase")
    parser.add_option("--keep_data","-k",help="keep old data without extracting repository",action="store_true")
    parser.add_option("--quiet","-q",help="does not produce a web page",action="store_true")
    parser.add_option("--makefile","-m",help="select Makefile (Makefile, Makefile_develop or nag Makefiles), none if unspeficied"\
                      ,dest="makefile")
    parser.add_option("--notify","-n",help="notify module keepers",action="store_true")
    
    (options, args) = parser.parse_args()

    if options.singleCase and not options.singleTarget:
        print("option --case assumes option --target is selected as well")

    if options.makefile:
        found = False # default
        for m in makefiles:
            if m == options.makefile:
                makefiles = [m]
                found = True
        if not found:
            raise('option --makefile or -m expect a valid makefile')
                
    tester = Tester()
    tester.run()
    if options.verbose:
        print("program completed.")

