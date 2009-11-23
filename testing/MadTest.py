#!/usr/bin/python
import optparse
import os
import re

class State: # the state of a test
    pass

init = State()
incomplete = State() # missing resources
ready = State() # all resources present
running = State()
aborted = State()
completed = State()

rootDir = "/afs/cern.ch/user/n/nougaret/scratch0/mad-automation/madX-examples/REF"    

# all files with absolute path

def lookForDependantResources(test,filename,lvl=0): # returns a list of filenames with absolute path
    resources = [] # list to be returned
    try:
        f = open(filename,'r')
        # regex to input a file (not anchored as an instruction may precede)
        commonPatterns = (re.compile(r'readmytable,?[\s\t]*file[\s\t]*=[\s\t]*[\"\']?([\w\.\_\d/]+)[\"\']?[\s\t]*;'),\
                          re.compile(r'call,?[\s\t]*file[\s\t]*=[\s\t]*[\"\']?([\w\.\_\d/]+)[\"\']?[\s\t]*;'),\
                          re.compile(r'readtable,?[\s\t]*file[\s\t]*=[\s\t]*[\"\']?([\w\.\_\d/]+)[\"\']?[\s\t]*;'))
        for line in f.readlines():
            line = line.rstrip('\n')
            if len(line.lstrip())>1 and line.lstrip()[0]=='!': # if line is a comment discard it
                continue
            k = line.find('=')
            if not k == -1:
                # put to lower case everything in front of the equal sign to simply patterns
                #print "before:"+line
                line = line[0:k].lower() +"="+ line[k+1:]
                #print "after:"+line
            for p in commonPatterns:
                m = p.search(line)
                if m:
                    resource = m.group(1)
                    if not resource[0] == '/': # relative path to be converted to absolute
                        dirPrefix = filename[0:filename.rfind('/')]                      
                        resource = dirPrefix+'/'+resource
                        resource = os.path.normpath(resource) # normalize the pathname
                    else:
                        pass
                        #print "'"+resource[0:2]+"'"
                    if options.verbose:
                        print("at "+str(lvl)+", resource '"+m.group(1)+"' has absolute path: '"+resource+"'")
                    resources.append(resource)
                    # now dig further
                    additionalResources = lookForDependantResources(test,resource,lvl+1)
                    for a in additionalResources:
                        resources.append(a)
        f.close()
    except: # failed to open the file: this may be a temporary file to be created at run time, any way one should
        # mark the resource as missing
        test.missing_resources.append(filename)
        if test.state == init:
            test.state = incomplete;
    return resources

class Test:
    tests = []

    def __init__(self,name,program,input,output,subdirectory=0):
        self.name = name
        self.program = program
        self.input = input
        self.output = output
        self.resources = []
        self.missing_resources = []
        self.state = init # initial state
        
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
                    self.resources = lookForDependantResources(self,entry[0]+'/'+self.input) # provide expanded filename
                    return True # found the input madx file
        return False
        
class Tester:

    makefiles = ['mk1','mk2','mk3']

    def __init__(self):
        pass

    def run(self):

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

        page = WebPage('testpage.htm')
        page.output()
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
        self.contents += '</head>'
        
    def body(self):
        self.contents += '<body>\n'
        self.contents += '<table>\n'
        for t in Test.tests:
            self.contents += '<tr><td>'+t.name+'</td><td>'+\
                             t.program+'&lt;'+t.input+'&gt;'+t.output+'</td></tr>\n'
            for r in t.resources:
                self.contents += '<tr><td>'+r+'</td></tr>'
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
        os.system('firefox ' + self.name)

if __name__ == "__main__":
    usage = "%prog [options]"
    parser = optparse.OptionParser(usage)
    parser.add_option("--verbose","-v",help="display messages useful for debug", action="store_true")
    parser.add_option("--target","-t",help="select a specific target to run the test",dest="singleTarget")
    parser.add_option("--case","-c",help="select a specific test-case to run the test, in combination with --target",\
                      dest="singleCase")
    (options, args) = parser.parse_args()
    if options.singleCase and not options.singleTarget:
        print("option --case assumes option --target is selected as well")
    tester = Tester()
    tester.run()
    if options.verbose:
        print("program completed.")

