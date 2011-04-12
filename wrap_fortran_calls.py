#!/usr/bin/python
# this script must work on Linux and Windows

import os
import re
import optparse
import sys

class subroutine: # Fortran subroutine related information
    subroutines = [] #  ordered
    subroutinesCalledFromC = []
    def __init__(self,file,name,argString,options):           
        self.name = name
        self.arguments = [] # ordered as in signature which may differ from order of parameter types in Fortran code
        self.argType = {} # dictionary
        self.fortranFile = file
        self.isCalledFromC = False # default
        argString = argString.replace(' ','') # remove blanks
        if argString == '':
            self.nArguments = 0
        else:
            parts = argString.split(',')
            self.nArguments = len(parts)
            for argName in parts:
                self.arguments.append(argName)
                self.argType[argName] = 'undefined_type' # default
        self.argString = argString

        redefined = False # by defaulr
        for i,sub in enumerate(subroutine.subroutines):
            if sub.name == name: # prevent double definitions
                if options.verbose:
                    print("warning: overwrite wrapper for subroutine '"+name+"' in file '"+file+\
                          "', previous definition from '" +\
                          sub.fortranFile +\
                          "' will be lost.")
                redefined = True
        if redefined:
            subroutine.subroutines[i] = self
        else:
            subroutine.subroutines.append(self)
        
    def calledFromC(self):
        self.isCalledFromC = True
        alreadyStored = False # default
        for i,sub in enumerate(subroutine.subroutinesCalledFromC):
            if sub.name == self.name:
                alreadyStored = True
        if alreadyStored:
            pass
        else:
            subroutine.subroutinesCalledFromC.append(self)
    def parsed(self): # returns true when all arguments have been recognized in terms of type and name
        for a in self.arguments:
            if self.argType[a] == 'undefined_type': # the default
                return False
        return True
    def registerArgumentType(self,cType,cName):
        if cName in self.argType:
            self.argType[cName] = cType
        else:
            pass # ignore this variable which is actually not a parameter of the function
    def generateCcode(self,file):
        file.write("/* Wrap '"+self.name+"' defined in '"+self.fortranFile+"' */\n")
        file.write('void '+self.name+'_wrapper(')
        if self.nArguments>0:
            for a in self.arguments:
                if self.name == "res_index" and a == "indexa":
                    file.write("int indexa[4][1000],")
                else:
                    if not a == self.arguments[-1]: # not yet the last one
                        file.write(self.argType[a]+' '+a+',')
                    else:
                        file.write(self.argType[a]+' '+a+'){\n')
        else:
            file.write('){\n')

        file.write('\tfflush(stdout);\n')
        file.write('\t'+self.name+'_(')

        if self.nArguments>0:
            for a in self.arguments:
                if not a == self.arguments[-1]: # not yet the last one
                    file.write(a+',')
                else:
                    file.write(a+');\n')
        else:
            file.write(');\n')

        file.write('\tcall_fortran_flush_("",0);\n')
        file.write("}\n")

    def generateWrapperHeaderCode(self,file):
        file.write("#define " + self.name+'_ ' + self.name + '_wrapper\n')

    def generatePrototypeCode(self,file):
        file.write("/* Wrap '"+self.name+"' defined in '"+self.fortranFile+"'*/\n")
        file.write("void "+self.name+"_(")
        if self.nArguments > 0:
            for a in self.arguments:
                if self.name == "res_index" and a == "indexa":
                    file.write("int indexa[4][1000],")
                else:
                    if not a == self.arguments[-1]: # not yet the last one
                        file.write(self.argType[a]+' '+a+',')
                    else:
                        file.write(self.argType[a]+' '+a+');\n')
        else:
            file.write(');\n')
            
    def generateWrapperPrototypeCode(self,file):
        file.write("/* Wrap '"+self.name+"' defined in '"+self.fortranFile+"'*/\n")
        file.write("void "+self.name+"_wrapper(")
        if self.nArguments > 0:
            for a in self.arguments:
                if self.name == "res_index" and a == "indexa":
                    file.write("int indexa[4][1000],")
                else:
                    if not a == self.arguments[-1]: # not yet the last one
                        file.write(self.argType[a]+' '+a+',')
                    else:
                        file.write(self.argType[a]+' '+a+');\n')
        else:
            file.write(');\n')        


def main():
    # process the command line arguments and options
    usage = '%prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-c','--clean',help='deletes all files created by this script',action='store_true')
    parser.add_option('-v','--verbose',help='display messages that may be useful for debugging',action='store_true')
    parser.add_option('-o','--outdir',
            type=str,
            default='',
            help='Output directory for the generated code')
    (options,args) = parser.parse_args()
    
    prefix = 'fortran'
    # YIL Added optional output directory
    if options.outdir:
        prefix = options.outdir+'/'+prefix
    
    wrapper_filename = prefix + '_wrappers.c'
    wrapper_header_filename = prefix + '_wrappers.h'
    prototypes_filename = prefix + '_prototypes.h'
    wrappers_prototypes_filename = prefix + '_wrappers_prototypes.h'

    targetFiles = [wrapper_filename, wrapper_header_filename, prototypes_filename,wrappers_prototypes_filename]

    if options.clean:
        for file in targetFiles:
            try:
                os.remove(file)
                print("delete " + file)            
            except:
                pass
        sys.exit()
        

    Cfiles = []
    Hfiles = []
    F90files = [] # Fortran files
    FIfiles = [] # Fortran included files

    allFiles = os.listdir('.')
    if options.outdir:
      allFiles.extend(os.listdir(options.outdir))

    for file in allFiles:
        if file[-2:] == '.c':
            Cfiles.append(file)
        elif file[-2:] == '.h':
            Hfiles.append(file)
        elif file[-4:] == '.f90':
            F90files.append(file)
        elif file[-3:] == '.fi':
            FIfiles.append(file)

    Csourcefiles = Cfiles + Hfiles
    FortranSourceFiles = F90files + FIfiles

    maybeFortranCalls = {}

    for file in Csourcefiles:
        f = open(file,'r')

        # preprocessing: get rid of C-style /*...*/ comments
        bigString = ''
        class State:
            pass

        # first step expand on one single line
        for line in f.readlines():
            bigString += line

        commentPattern = re.compile('/\*.+\*/')
        # eat all patterns where they appear
        bigString = re.sub(commentPattern,'',bigString) # remove all comments

        lines = bigString.split('\n')
            
        # special case where the function is called through a C macro of the form: #define f f_
        patternDefine = re.compile(r'^#define[\s\t]*([\w\d\_]+)[\s\t]+([\w\d\_]+).*$')
        patternCall = re.compile(r'[\s\t]+([\w\d_]+)_\((.+)') # handle case where second parenthesis runs over several lines
        for line in lines:
            m = patternDefine.match(line)
            if m:
                str1 = m.group(1)
                str2 = m.group(2)
                if (str1+'_')==str2:
                    # register the routine as being a Fortran call from C
                    maybeFortranCalls[str1] = str1
            else:
                m = patternCall.search(line)
                if m:
                    str1 = m.group(1)
                    maybeFortranCalls[str1]=str1
        f.close()

    for file in FortranSourceFiles:
        f = open(file,'r')
        #if options.verbose:
        #    print("currently scanning file "+file);
        subroutineResIndex =re.compile(r'^[\s\t]*subroutine[\s\t]+res_index\(skew,mynorder,myn1,myn2,indexa,mynres\)')
        subroutinePattern = re.compile(r'^[\s\t]*[sS][uU][bB][rR][oO][uU][tT][iI][nN][eE][\s\t]*(.+)\((.*)\)$')
        # WARNING: sometime we can have: subroutine SUBNAME only! (i.e. no parenthesis afterwards)
        subroutineShortPattern = re.compile(r'^[\s\t]*[sS][uU][bB][rR][oO][uU][tT][iI][nN][eE][\s\t]*(.+)$') # PROBABLY USELESS
        class State:
            pass

        # preprocessing: concatenate instructions that run across successive lines
        concatenatedLines = []
        firstLine = State()
        nextLine = State()
        state = firstLine
        for line in f.readlines():
            line = line.rstrip('\n\t ')
            line = line.lstrip('\t ')
            if len(line)<1:
                continue
            if state == firstLine:
                if line[-1] == '&':
                    concatenatedLine = line.rstrip(' &')
                    nextState = nextLine
                else:
                    concatenatedLine = line # untouched
                    concatenatedLines.append(concatenatedLine)                
                    nextState = firstLine
            elif state == nextLine:
                if line[0] == '&' or True: # 12 March 2009 line continuation change: the following line does not need to start with a '&'
                    if line[-1] == '&': # continues further
                        concatenatedLine += line.lstrip('&').rstrip(' &')
                        nextState = nextLine
                    else:
                        concatenatedLine += line.lstrip('&')
                        concatenatedLines.append(concatenatedLine)                    
                        nextState = firstLine
                else:
                    concatenatedLine = line # untouched
                    concatenatedLines.append(concatenatedLine)                    
                    nextState = firstLine
            state = nextState
            
        lookForSubroutine = State()
        lookForArguments = State()
        state = lookForSubroutine
        
        for line in concatenatedLines:
            # remove Fortran style comments
            if line.lstrip(' ')[0]=='!':
                continue # skip this line
                
            if line == concatenatedLines[0] and not state == lookForSubroutine:
                raise("failed to parse a subroutine's argument in the previous file") 
            
            if state == lookForSubroutine:

                # very special case: res_index in resindex.f90 relies on an include file resindex.fi
                # for its parameters. Let's handle it manually:
                m = subroutineResIndex.match(line)
                if m:
                    # print("found resindex")
                    res_index = subroutine(file,"res_index","skew,mynorder,myn1,myn2,indexa,mynres",options)
                    res_index.calledFromC()
                    # "int* skew,int* mynorder,int* myn1,int* myn2,int indexa[4][1000],int* mynres"                
                    # now define all parameters
                    res_index.registerArgumentType('int*','skew')
                    res_index.registerArgumentType('int*','mynorder')
                    res_index.registerArgumentType('int*','myn1')
                    res_index.registerArgumentType('int*','myn2')
                    res_index.registerArgumentType('int','indexa') # THIS ONE SHOULD BE "int indexa[4][1000]" => handled in codegen routines
                    res_index.registerArgumentType('int*','mynres')
                    continue # next in the for-loop

                # Some very special cases of subroutines that cause problem with argument parsing and that we want to discard alltogether

    #            skippedSubroutine = re.compile('subroutine trrun')
    #            m = skippedSubroutine.match(line);
    #            if m:
    #                continue # next in the for loop
                trrunSubroutine = re.compile('subroutine trrun') # let's try not to skip it actually
                m = trrunSubroutine.match(line);
                if m:
                    trrun = subroutine(file,"trrun","flag,turns,orbit0,oneturnmat,ibuf1,ibuf2, buf1, buf2,buf_dxt,buf_dyt,buf3,buf4,buf5,e_flag,ibuf3,buf6",options)
                    trrun.calledFromC()
                    trrun.registerArgumentType('int*','flag')
                    trrun.registerArgumentType('int*','turns')
                    trrun.registerArgumentType('double*','orbit0')
                    trrun.registerArgumentType('double*','oneturnmat')
                    trrun.registerArgumentType('int*','ibuf1')
                    trrun.registerArgumentType('int*','ibuf2')
                    trrun.registerArgumentType('double*','buf1')
                    trrun.registerArgumentType('double*','buf2')
                    trrun.registerArgumentType('double*','buf_dxt')
                    trrun.registerArgumentType('double*','buf_dyt')
                    trrun.registerArgumentType('double*','buf3')
                    trrun.registerArgumentType('double*','buf4')
                    trrun.registerArgumentType('double*','buf5')
                    trrun.registerArgumentType('int*','e_flag')
                    trrun.registerArgumentType('int*','ibuf3')
                    trrun.registerArgumentType('double*','buf6')
                    continue # next in the for loop




                
                        
                m = subroutinePattern.match(line) # no longer line.lower()
                if m:             
                    name = m.group(1)
                        
                    arguments = m.group(2)
                    newSubroutine = subroutine(file,name,arguments,options) # BUG: arguments names are lowered
                    
                    for k,v in maybeFortranCalls.iteritems(): # check if this Fortran subroutine is called from C   
                        if k == name:
                            newSubroutine.calledFromC()
                            if options.verbose:
                                print("found subroutine '"+newSubroutine.name+"', which is called from C")                        
                            if not arguments.lstrip(' \t').rstrip(' \t') == "":
                                state = lookForArguments # now ready to retreive the argument's types and names one-by-one (if any!)
                            else:
                                pass # this routine has no argument

                else: # try with the other pattern where we have no parenthesized arguments
                    m = subroutineShortPattern.match(line)
                    if m:
                        name = m.group(1)
                        arguments = ''
                        newSubroutine = subroutine(file,name,arguments,options)

                        for k,v in maybeFortranCalls.iteritems(): # check if this Fortran subroutine is called from C   
                            if k == name:
                                newSubroutine.calledFromC()
                                if options.verbose:
                                    print("found subroutine '"+newSubroutine.name+"', which is called from C")
                                if not arguments.lstrip(' \t').rstrip(' \t') == "":
                                    state = lookForArguments # ready to retreive argument's types and names one-by-one (if any!)
                                else:
                                    pass # this routine has no argument                    

            elif state == lookForArguments:


                # now extract the subroutine's signature
                
                typePattern = re.compile(r'^[\s\t]*([rR][eE][aA][lL]\(?[\w\d\(\)]*\)?|[dD][oO][uU][bB][lL][eE] [pP][rR][eE][cC][iI][sS][iI][oO][nN]|'+\
                                        '[iI][nN][tT][eE][gG][eE][rR]|[cC][hH][aA][rR][aA][cC][tT][eE][rR][*\(]?\d*\)?|'+\
                                        '[cC][oO][mM][pP][lL][eE][xX]|[lL][oO][gG][iI][cC][aA][lL]\(?[\w\d]*\)?)[\s\t]+(.+)$')
                m = typePattern.match(line) # assume types are always in lower cases, while parameters can be mixed cases
                if m:

                    argType = m.group(1)
                    argString = m.group(2)
                    if argType == "double precision":
                        cType = "double*"
                    elif argType == "real":
                        cType = "float*"
                    elif argType == "real(kind(1d0))":
                        cType = "double*"
                    elif argType == "integer":
                        cType = "int*"
                    elif argType[0:len("logical")] == "logical": # accounts for logical as well as logical(lp) etc...
                        cType = "int*"
                    elif argType[0:len("character")] == "character": # accounts for character as well as character*16 etc...
                        cType = "char*"
                    else:
                        raise("unknown type encountered:" + argType)

                    argString = argString.replace(' ','')
                    argString = argString.replace('::','') # sometimes, we get type :: variable

                    # in Fortran, arrays appear as a(b) or a (b,c) where the type of the array is retreived from a, hence let's discard parenthesis
                    consumeParenthesis = State()
                    normal = State()
                    mode = normal
                    cleanString = ''
                    for c in argString:
                        if c=='(':
                            mode = consumeParenthesis
                        elif c==')':
                            mode = normal
                            continue
                        if mode == normal:
                            cleanString += c
                            
                    argString = cleanString


                    # remove fortran style comment at the end of the line
                    if not argString.find('!')==-1:
                        argString = argString[0:argString.find('!')]

                    args = argString.split(',') # split series of parameter definitions sharing the same type
                    for arg in args:
                        cName = arg
                        newSubroutine.registerArgumentType(cType,cName) # the name should be already known by subroutine
                        # otherwise, this is not a function's parameter and it must be ignored

                else:
                    m = subroutinePattern.match(line)
                    if m:
                        print("WARNING: encountered new subroutine definition of '"+\
                            m.group(1)+\
                            "' in file '"+\
                            file+\
                            "' before being able to parse parameters of "+newSubroutine.name)
                    
                if newSubroutine.parsed():
                    state = lookForSubroutine # now ready to look for the next subroutine in the current file
                    
                    
            else:
                raise("FSM state must be either 'lookForSubroutine' or 'lookForArguments'")
                

        f.close()
    
    fWrapperCode = open(wrapper_filename,"w")
    fWrapperCode.write("/* set of "+str(len(subroutine.subroutinesCalledFromC))+" wrappers to synchronize FORTRAN and C stdout buffers */\n")
    fWrapperCode.write("/* when crossing the border upon calling FORTRAN from C. */\n\n")
    fWrapperCode.write("#include <stdio.h>\n")
    fWrapperCode.write('#include "fortran_prototypes.h"\n\n')
    fWrapperCode.write("extern void call_fortran_flush_(char *, int);\n\n")

    fWrapperHeader = open(wrapper_header_filename,"w")
    fWrapperHeader.write(\
        '#ifndef _FORTRAN_WRAPPERS_H\n'+\
        '#define _FORTRAN_WRAPPERS_H\n'+\
        '/* redirect FORTRAN calls to wrappers that synchronize FORTRAN and C stdout buffering */\n'+\
        '/* when crossing the border upon calling FORTRAN from C. */\n'+\
        '#include "fortran_wrappers_prototypes.h"\n');


    fPrototypes = open(prototypes_filename,"w")
    fPrototypes.write(\
        '/* to avoid warnings of implicit declarations from fortran_wrappers.c */\n'+\
        '#ifndef _FORTRAN_PROTOTYPES_H\n'+\
        '#define _FORTRAN_PROTOTYPES_H\n')

    fWrappersPrototypes = open(wrappers_prototypes_filename,"w")
    fWrappersPrototypes.write(\
        '/* to avoid warnings of implicit declaration from code calling the functions below */\n'+\
        '#ifndef _FORTRAN_WRAPPERS_PROTOTYPES_H\n'+\
        '#define _FORTRAN_WRAPPERS_PROTOTYPES_H\n')

    if options.verbose:
        print("number of subroutines: "+str(len(subroutine.subroutines)))
        print("number of subroutines called from C: "+str(len(subroutine.subroutinesCalledFromC)))

    alreadySeen = []      
    for  s in subroutine.subroutines:
        for a in alreadySeen:
            if s == a:
                print "warning saw '"+s.name+"' twice"
        alreadySeen.append(s)

    number = 0
    for s in subroutine.subroutines:
        if s.parsed() and s.isCalledFromC:
            number += 1
            s.generateCcode(fWrapperCode)
            s.generateWrapperHeaderCode(fWrapperHeader)
            s.generatePrototypeCode(fPrototypes)
            s.generateWrapperPrototypeCode(fWrappersPrototypes)
        elif (not s.parsed()) and s.isCalledFromC:
            print("skipped "+s.name+" defined in "+ s.fortranFile + ", due to undefined arguments:")
            for arg in s.arguments:
                print("\targ '"+arg+"' is of type '"+s.argType[arg]+"'")

    if options.verbose:
        print("number of wrapped subroutines: "+str(number))

    fWrapperCode.close()

    fWrapperHeader.write('#endif\n')
    fWrapperHeader.close()

    fPrototypes.write('#endif\n')
    fPrototypes.close()

    fWrappersPrototypes.write('#endif\n')
    fWrappersPrototypes.close()


    if options.verbose:
        print('completed.')
if __name__ == '__main__':
    main()
