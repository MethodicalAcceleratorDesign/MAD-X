#!/usr/bin/python


import datetime
import optparse
import os
import stat
import re
import shutil


from Notify import notify


class Logger: # in more recent versions, Python features a logger
    def __init__(self,name):
        self.name = name
        if os.path.exists(self.name):
            os.remove(self.name)

    def log(self,msg):
        f = open(self.name,"a") # append
        f.write(msg+'\n')
        f.close() # close after each log message

website =   "/afs/cern.ch/user/n/nougaret/www/mad"
reportDir = "/afs/cern.ch/user/n/nougaret/scratch1/mad-automation"
# 12 april 2009 fix the directory in which to extract MAD-X
#currentDir = os.getcwd()
currentDir = reportDir
extractDir = currentDir+'/MadSvnExtract'

notify('jean-luc','MadSvnExtract','MadSvnExtract will be created as ' + extractDir)

class Repository:
    def __init__(self):
        pass
    def checkout(self,releaseTag):
#        command = "cvs -d " + Repository.repoDir + " checkout -r" + releaseTag + " madX"
        command = "svn co svn+ssh://svn.cern.ch/reps/madx/tags/" + releaseTag + " " + extractDir
        os.system(command)

class Build:
    builds = []
    def __init__(self,name,command,outcome):
        self.name = name # 'Makefile', 'Makefile_dev' or 'Makefile_nag'
        self.command = command # the build command
        self.outcome = outcome
        Build.builds.append(self)

# first set environment variables for lf95 and NAG compilers
# this is necessary for the acron job which does not perform a login
# that would set the variables transparently.

# 27 september: this script is no longer called from an acron job, therefore
# no longer need to reset the environment variables within the scripts.
# instead we inherit everything from ~/.cshrc on pcslux99

def main():

# temporary: for the Intel compiler, lines split after the 80th character
# unless we set up a specific environment variable to set the line length
# note that this variable will be transmitted to the child processes
# launched through os.sytem() calls.
# all 'madx < input > output' calls where madx is compiled with the ifort
# Intel compiler will take FORT_FMT_RECL into account.
#    os.environ['FORT_FMT_RECL']='256' # lines split after 256 chars instead of
# THIS DOES NOT SEEM TO WORK 13 april 2010 - force the variable when invoking the makefile instead
# the default 80    
         
    LOG_FILENAME = reportDir + "/MadBuild_Report.txt"
    logger = Logger(LOG_FILENAME)
    logger.log('first logging message')    
    logger.log('start logging at '+(datetime.datetime.now()).ctime())
                 
    compilationOK = True # default (assumes all the compilations proceed till completion)

    globalStartTime =  (datetime.datetime.now()).ctime()
    usage = "%prog [options]"
    parser = optparse.OptionParser(usage)
    parser.add_option('--release','-r',help="release tag",dest="releaseTag")
    parser.add_option('--dev','-d',help="construct executable with Makefile_develop",action="store_true")
    parser.add_option('--nag','-n',help="construct executable with Makefile_nag",action="store_true")
    (options,args) = parser.parse_args()

    if options.dev and not options.nag:
        logger.log("option --dev selected")
    elif options.nag and not options.dev:
        logger.log("option --nag selected")
    elif options.dev and options.nag:
        logger.log("options --dev and --nag selected")
    else:
        logger.log("only compile madx_Makefile")

    if not options.releaseTag:
        raise("except a release tag to be specified")
    else: # is the release tag well formed?
        releasePattern = re.compile(r'^madX\-\d+_\d{2}_\d{2}(\_dev)?$')
        if not releasePattern.match(options.releaseTag):
            raise("release tag is ill-formed: it should be like 'madX-4_01_01' instead of '"+options.releaseTag+"'")

    # 23 juin 2010
    #os.environ['PATH'] = os.environ['PATH'] +\
    #                     ":/afs/cern.ch/sw/fortran/nag/f95.5.361/bin:/afs/cern.ch/sw/fortran/lahey/lf9562/bin"
    #
    #if os.getenv('LD_LIBRARY_PATH')==None: # would cause a key error at runtime
    #    os.environ['LD_LIBRARY_PATH'] = ":/afs/cern.ch/sw/fortran/lahey/lf9562/lib"
    #else:
    #    os.environ['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH'] + ":/afs/cern.ch/sw/fortran/lahey/lf9562/lib"

    makefiles = ['Makefile']

    if options.nag:
        # 23 juin 2010
        #        os.environ['NAG95_ROOT'] = "/afs/cern.ch/sw/fortran/nag/f95.5.361" # flexlm license manager for NAG compiler
        #os.environ['LM_LICENSE_FILE'] = "/afs/cern.ch/sw/fortran/nag/f95.5.361/license.dat"
        # 27 september 2010: rely on the path and env set in ~/.cshrc instead of setting environment variables
        #os.environ['NAG95_ROOT'] = "/usr/local/lib/NAG_Fortran"
        #os.environ['NAG_KUSARI_FILE'] = os.environ['NAG95_ROOT']+"/license.dat"
        makefiles.append('Makefile_nag')

        #os.environ['PATH'] = ".:/usr/local/lib/NAG_FORTRAN:"+os.environ['PATH']

    if options.dev:
        makefiles.append('Makefile_develop')

    # if the directory in which to extract the CVS already exists, then delete it
    if os.path.exists(extractDir):
        shutil.rmtree(extractDir)
    # now create the directory in which to extract the CVS    
    #os.mkdir(extractDir)
    #os.chdir(extractDir)
    repo = Repository()
    repo.checkout(options.releaseTag) # creates or overwrites extractDir
    os.chdir(extractDir)
    os.chdir('./madX')

    for m in makefiles:

        # do a make clean and remove any preexisting target
        if os.path.exists('./madx'):
            os.remove('./madx')
            logger.log('remove preexisting ./madx executable')
        else:
            logger.log('no preexisting ./madx executable')
        os.system('make clean')
        logger.log("invoked 'make clean' at time "+(datetime.datetime.now()).ctime())

        #print("now to compile "+m+" in " + os.getcwd())
        if m == 'Makefile':
            #invocation = "make f95=/opt/intel/Compiler/11.0/081/bin/ia32/ifort DEBUG=NO FORT_FMT_RECL=256 madx"
            invocation = "make f95=ifort DEBUG=NO FORT_FMT_RECL=256 madx"
        elif m == 'Makefile_develop':
            invocation = "make f95=lf95 DEBUG=YES madx" # the Lahey Fortran compiler
        elif m == 'Makefile_nag':
            invocation = "make f95=f95 DEBUG=YES madx" # the NAG compiler
        else:
            raise("should never reach this point")
        print "make invocation=",invocation
        #notify('jean-luc','Point 3.2', 'makefile='+m) 
        output = "./makeResult"
        #os.system(invocation)
        startTime = (datetime.datetime.now()).ctime()
        #notify('jean-luc','now in directory', 'os.getcwd()='+os.getcwd()+" at time"+(datetime.datetime.now()).ctime())              
        #notify('jean-luc','Point 3.3', 'now to issue command='+invocation+">& "+output+" at time"+(datetime.datetime.now()).ctime())
        logger.log('now located in '+os.getcwd())
        logger.log('now to invoke command='+invocation+">& "+output+" at time "+(datetime.datetime.now()).ctime())

        os.system(invocation+">& "+output)

        endTime = (datetime.datetime.now()).ctime()

        logger.log("compilation command completed at time : " + endTime)

        f = open(output,'r')
        lines = f.readlines()
        f.close()
        #notify('jean-luc','Point 4', '') 
        # if the target madx is absent, compilation failed
        if os.path.exists('./madx'):
            # success or warning
            outcome = 'success' # for the time-being
            outcome = lookForWarnings(lines)
            executable = './madx_'+m # new name of the executable
            shutil.copyfile('./madx',executable)
            if os.path.exists(executable):
                logger.log('.madx has been copied into '+executable+" under "+os.getcwd())
            else:
                logger.log("can't find executable "+executable+" which must have just been created")
            os.chmod(executable,stat.S_IRWXU) # read, write, and execute permission by owner
            logger.log("a ./madx executable was succesfully created at "+(datetime.datetime.now()).ctime())
            # very important: grant execution rights to the script
            #print("succesfully compiled for "+m)
        else: # failure
            outcome = 'failure'
            compilationOK = False # overwrite
            logger.log("missing ./madx executable")
            #print("failed to compile for "+m)

        b = Build(m,invocation,outcome)
        #notify('jean-luc','Point 5', 'makefile='+m) 
        page = BuildWebPage(website+'/build_'+m+'_madx.htm',m,lines,outcome,startTime,endTime)
        page.output()

    #notify('jean-luc','Point 6', '') 
    globalEndTime = (datetime.datetime.now()).ctime()       
    # create the main HTML page from which the three above HTML pages can be reached    
    page = MainWebPage(website+'/build.htm',globalStartTime,globalEndTime,Build.builds)
    page.output()

    # Finally output the return status on stdout to be processed by the calling MadBuildAndTest.py
    if not compilationOK:
        notify('jean-luc','compilation failed','failed to compile release = '+options.releaseTag)
    else:
        notify('jean-luc','compilation succeeded','managed to compile last release = '+options.releaseTag)
    print compilationOK

    notify('jean-luc','reached the end of MadBuild.py','compilationOK=' + str(compilationOK))



def lookForWarnings(lines): 
    warningPattern = re.compile(r'[\s\t]*[\w\d_\-:]+[\s\t]+[Ww]arning(s?)')
    for line in lines:
        m = warningPattern.search(line)
        if m:
            if m.lastindex == 1:
                return('warning')
    return('success')
        
    

class MainWebPage: # the main web page that lists all the outcomes for the successive builds
    
    def __init__(self,name,startTime,endTime,builds):
        self.name = name
        self.startTime = startTime
        self.endTime = endTime
        self.builds = builds # a list of Build objects

    def header(self):
        self.contents += '<DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">\n'
        self.contents += '<html>\n'
        self.contents += '<head>\n'
        self.contents += '<title>MAD build result page</title>\n'
        self.contents += '<link rel=stylesheet href="./MadTestWebStyle.css" type="text/css">\n'
        self.contents += '</head>\n'
        
    def body(self):
        self.contents += '<body>\n'
        self.contents += 'Build started '+ self.startTime + ', ended  ' + self.endTime
        self.contents += '<table>\n'      
        for b in self.builds:
            weblink = "build_"+b.name+"_madx.htm"
            self.contents += '<tr class="'+b.outcome+'"><td>'+b.name+'</td><td>'+b.command+'</td><td><a href="'+weblink+'">'+b.outcome+'</a></td></tr>\n'
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

class BuildWebPage(MainWebPage):
    def __init__(self,name,makefile,lines,returnStatus,startTime,endTime):
        self.name = name
        self.returnStatus = returnStatus
        self.buildResult = lines
        self.startTime = startTime
        self.endTime = endTime
        self.makefile = makefile

    def header(self):
        self.contents += '<DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">\n'
        self.contents += '<html>\n'
        self.contents += '<head>\n'
        self.contents += '<title>MAD build main page</title>\n'
        self.contents += '<link rel=stylesheet href="./MadTestWebStyle.css" type="text/css">\n'
        self.contents += '</head>\n'

    def body(self):
        self.contents += '<body>\n'
        self.contents += '<p>Build started '+self.startTime+', ended '+self.endTime+'</p>\n'
        self.contents += '<table width="75%" border="0">\n'
        # top coloured banner
        self.contents += '<tr class='+self.returnStatus+'><td width="80%">Outcome of build process for \''+self.makefile+'\'</td><td width="20%">'+\
                         self.returnStatus+"</td></tr>\n"
        for l in self.buildResult:
            self.contents += '<tr><td colspan="2">'+l.rstrip("\n")+'</td></tr>\n'
        self.contents += "</table>\n"
        self.contents += '</body>\n'        
        
if __name__ == "__main__":
    main()

