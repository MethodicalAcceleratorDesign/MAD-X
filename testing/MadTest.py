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
import time
import datetime

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
        # check the resource indeed exists. If not this a temp file meant to be created at runtime
        if not os.path.exists(source):
            self.type = 'runtime'
            if options.verbose:
                print("WARNING resource file "+source+" is probably meant to be created at runtime")
        else:
            self.type = 'static'
        self.destinations = destinations # the expanded filename of the file in local work directory
        if options.verbose:
            print("create resource '"+name+"' retreived from '"+source+"'")

    def lookForDependantResources(test,filename,lvl=0,leafDir=0):
        # returns a list of filenames with absolute path

        resources = [] # list to be returned
        
        if lvl==0:
            source = filename
            name = source[source.rfind('/')+1:]
            parts = source.split('/')
            leafDir =  parts[-2]
            print("leafDir="+leafDir)
            # the leaf directory that contain the example and that should be renamed by test.testcaseDir
            destinations = []
            for m in makefiles:
                destination = source.replace(rootDir,testingDir+'/'+m) # replace prefix between source and destination
                match = re.match(r'^test_\d+(_.+)?$',test.testcaseDir)
                if not match:
                    raise("pattern should always match")
                if match.lastindex == 1: # testcaseDir like test_1_LHC
                    destination = destination.replace('/'+leafDir+'/','/'+test.testcaseDir+'/') # testcaseDir like test_1_LHC                    
                else:
                    destination = destination.replace('/'+leafDir+'/','/'+leafDir+'/'+test.testcaseDir+'/') # testcaseDir like test_1
                    
                print("destination="+destination+" for source="+source+", leafDir="+leafDir)
                destinations.append(destination)            
            resource = Resource(name,source,destinations)  # primary resource is the MAD-X input file itself
            resources.append(resource)
            
        try:
            f = open(filename,'r')
            # regex to input a file (not anchored as an instruction may precede)
            commonPatterns = (re.compile(r'readmytable,?[\s\t]*file[\s\t]*=[\s\t]*[\"\']?([\w\.\_\-\d/]+)[\"\']?[\s\t]*;'),\
                              re.compile(r'call,?[\s\t]*file[\s\t]*=[\s\t]*[\"\']?([\w\.\_\-\d/]+)[\"\']?[\s\t]*;'),\
                              re.compile(r'readtable,?[\s\t]*file[\s\t]*=[\s\t]*[\"\']?([\w\.\_\-\d/]+)[\"\']?[\s\t]*;'))
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
                            match = re.match(r'^test_\d+(_.+)?$',test.testcaseDir)
                            if not match:
                                raise("pattern should always match")
                            if match.lastindex == 1: # testcaseDir like test_1_LHC
                                destination = destination.replace('/'+leafDir+'/','/'+test.testcaseDir+'/') # testcaseDir like test_1_LHC
                            else:
                                destination = destination.replace('/'+leafDir+'/','/'+leafDir+'/'+test.testcaseDir+'/') # testcaseDir like test_1
                            destinations.append(destination)
                          
                        resource = Resource(name,source,destinations)
                        resources.append(resource)
                        
                        # now dig further
                        additionalResources = Resource.lookForDependantResources(test,source,lvl+1,\
                                                                                 leafDir)                     
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

class Output:
    def __init__(self,name,outcome,pageLink):
        self.name = name
        self.outcome = outcome
        self.weblink = pageLink

class Test: # a test case
    tests = []

    def __init__(self,name,program,input,output,subdirectory):

        self.name = name # the target
        self.program = program
        self.input = input # the main input file
        self.output = output # the main output file
        self.resources = [] # secondary inputs
        self.outputs = [] # secondary outputs
        self.missing_resources = []
        self.state = init # initial state
        
        Target.registerTest(name, self)

        if subdirectory == 0 or subdirectory == '':
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
            if r.type == 'runtime':
                # this resource is meant to be created at run time => do not attempt to copy as it does not exists yet
                continue
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
            print("looking for /test_1/: script is '"+script+"'")            
            scriptDir = script[:script.rfind('/')]
            self.topDir = scriptDir # for future reuse when comparing the output with the reference
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
                    # print("compare "+scriptDir+"/"+f+" with "+r.destinations[i])

                    if scriptDir+"/"+f == r.destinations[i]: # this is an (input) resource file
                        type = input
                if type == input:
                    shutil.move(f,'./input')
                else:
                    shutil.move(f,'./output')
            
            os.chdir(currentDir) # back to the initial work directory

    def compareOutputWithReference(self):
        # only for the main Makefile
        # must pick-up all files under the /output directory and compare them with the reference
        outputDir = self.topDir+'/output'
        print "outputDir="+outputDir
        files = os.listdir(outputDir)
        
        for fname in files: # the short name of the file without its path prefix
            try:
                os.remove('./tempfile')
            except:
                pass

            createdFile = outputDir + '/' + fname  # fname with adequate prefix
            referenceFile = self.resources[0].source[:self.resources[0].source.rfind('/')] +'/'+ fname # fname with adequate prefix -> take the madx input host directory in the CVS and add the filename

            # make sure the reference file exists, otherwise we shall omit to look for differences
            if os.path.exists(referenceFile):

                #print("name of createdFile="+createdFile)
                #print("name of referenceFile="+referenceFile)

                # specific case when the HTML file name is of the form XX.map or XX.map.htm
                # webserver will fail to display the HTML although one can open it from the webfolder...
                # to overcome this limitation, we need to juggle with the HTML file name
                fname = fname.replace('map','maAap')
            
                htmlFile = "./temp.html" # output HTML file, to be delivered to the web site...
                weblink = "./DiffResult_" + self.name + "_" + self.testcaseDir + "_" + fname + ".htm" # again test.name stands for the target
                htmlFile = htmlRootDir+"/details/"+weblink

                os.system('./MadDiff.pl '+createdFile+' '+referenceFile+' ' +htmlFile+' > ./tempfile')

                tf = open("./tempfile","r")
                outcome = tf.readlines()[0] # single line

                if outcome == 'failure':
                    pass
                    #print("this is a failure")
                elif outcome == 'warning':
                    pass
                    #print("this is a warning")
                elif outcome == 'quasi-success':
                    pass
                    #print("this is a quasi-sucess")
                elif outcome == 'success':
                    pass
                    #print("this is a success")
                else:
                    raise("should never reach this point, outcome = "+outcome)

                tf.close()

                os.remove('./tempfile')

            else:
                outcome = 'omit'

            # store the short name of the output file, together with the outcome of the comparison

            # skip .ps and .eps files
            if fname[-3:]=='.ps' or fname[-4:]=='.eps':
                pass # skip
            else:
                out = Output(fname,outcome,weblink)
                self.outputs.append(out)

                # store output file information in the test object for subsequent reuse
                # by the web page output
            

            
class Tester:



    def __init__(self):
        pass

    def run(self):

        page = WebPage(mainHtmlPage)
        
        # first extract examples from the repository
        if not options.keep_data:
            rep = Repository()
            rep.checkout()

        testinfoPattern = re.compile(r'^[\s\t]*(.+)[\s\t]*<[\s\t]*(.+)[\s\t]*>[\s\t]*([\w\.\d\-\_]+)[\s\t]*,?[\s\t]*'+\
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

            if options.skippedTarget:
                if target == options.skippedTarget:
                    print("skip "+options.skippedTarget)
                    continue
                
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

        try:
            shutil.rmtree(testingDir)
        except:
            pass # directory absent

        page.output() # refresh the web page

        # now populate testingDir with all sources and associated resources
        for t in Test.tests:
            t.copyResourcesToTestingDir()


        # now run the tests
        for t in Test.tests:

            t.state = running
            page.output() # to mark the current test as running
            t.run()
            t.compareOutputWithReference()
            t.state = completed # or aborted?
            page.output() # refresh the web page

        # notify module keepers if required
        if options.notify:
            Notify.notify("jean-luc","test completion","test completed.") # for the time-being

            page.output() # refresh the web page for the last time


class WebPage:
    
    def __init__(self,name):
        self.name = name

    def header(self):
        self.contents += '<DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">\n'
        self.contents += '<html>\n'
        self.contents += '<head>\n'
        # remove the auto refresh and rely on the user to refresh on demand
        #self.contents += '<meta http-equiv="refresh" content="5" >' # page reloads itself every 5 seconds
        # this tag should only be present while the test is running, and be absent when completed.
        self.contents += '<title>MAD testing main page</title>'
        self.contents += '<link rel=stylesheet href="./MadTestWebStyle.css" type="text/css">'
        self.contents += '</head>'
        
    def body(self):
        self.contents += '<body>\n'
        self.contents += (datetime.datetime.now()).ctime()
        self.contents += '<table>\n'
        for target in Target.targets:
            self.contents += '<tr class="test_target"><td colspan="2"><div align="center"><strong>'+target.name+'</strong></div></td></tr>\n'  
            for i,t in enumerate(target.tests):
                self.contents += '<tr class="test_case"><td width=\"80%\">'+t.testcaseDir+\
                                 ': '+t.program +'&lt;'+t.input+'&gt;'+t.output+\
                                 '</td><td width=\"20%\"><table width=\"100%\" style=\"text-align: center\"><tr>'+\
                                 'DEVRES'+'NAGRES'+\
                                 '</tr></table></td></tr>\n';
                if t.state == init or t.state == incomplete:                
                    for r in t.resources:
                        self.contents += '<tr><td>'+r.name+'</td></tr>\n'
                elif t.state == running:
                    self.contents += '<tr><td width="80%"><center><b>The script is now busy running this test.</center></td><td width="20%"><center><img src="underConstruction.gif" width="134" height="139"></b></center></td></tr>\n' # image showing running job
                elif t.state == completed or t.state == aborted:
                    for o in t.outputs:
                        self.contents += '<tr class="'+o.outcome+'"><td width=\"80%\">'+o.name+\
                                         '<td width="30%"><a href="./details/'+o.weblink+'">'+o.outcome+'</a></td></tr>\n'

                else:
                    raise("should never reach this point")
                    
                    
        self.contents += '<table>\n'
        self.contents += '<body>\n'
    def footer(self):
        self.contents += '</html>\n'

    def output(self):
        self.contents = ""
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
    parser.add_option("--skip","-s",help="skip a particular target causing trouble, to proceed with debugging",dest='skippedTarget')
    
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

