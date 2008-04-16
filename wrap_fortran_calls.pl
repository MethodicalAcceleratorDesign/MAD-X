#!/usr/bin/perl
my $os;
if ($#ARGV>-1){
	if (@ARGV[0] eq 'MSDOS') {
		$os = 'MSDOS';
	} else {
		# in case the script is invoked on Windows, as from Makefile.bat,
		# one must pass an additionnal argument so that system calls
		# are issued as for MSDOS.
		print "only accept 'MSDOS' as optional argument. Exit!\n";
		exit;
	}
} else {
	$os = 'unix';
}

my $dbg = 0; # set to 1 to display messages

if ($os ne 'MSDOS') {
	@Cfiles = `ls *.c`;
	@Hfiles = `ls *.h`;
} else {
	# emulate the above in case we run the script on Windows
	@Cfiles = ();
	@Hfiles = ();
	@files_entries = `dir`; # dir *.c or dir *.h failed to return all files
	# and the results were different on DOS or Linux
	foreach $entry (@files_entries){
		if ($entry =~ /[\s\t]+([\w\d_\-]+\.c)/) {
			@Cfiles = (@Cfiles, $1);
		}
	}
	foreach $entry (@files_entries){
		if ($entry =~ /[\s\t]+([\w\d_\-]+\.h)/) {
			@Hfiles = (@Hfiles, $1);
		}
	}
	# at this stage @Cfiles and @Hfiles contain the list of C and H files
}

@Cfiles = (@Cfiles, @Hfiles);

# and what about the #define statements?
@maybeFortranCalls = ();
foreach $Cfile (@Cfiles){
    chomp $Cfile; # safer than chop: for DOS, no ending newline
	open(MYFILE,"<$Cfile");
	while(<MYFILE>){
		my $line = $_;

		# special case where the function is called though a C macro
		# of the sort: #define f f_
		if ($line =~ /^\#define[\s\t]*([\w\d_]+)[\s\t]+([\w\d_]+)/){
		    # print "# define $1 $2\n";
		    my $str1 = $1 . "_";
		    my $str2 = $2 ;
		    if ($str1 eq $str2) {
			#print "found a macro '\#define $1 $2' macro hiding '$1()' calls\n";
			#exit;
			# register the routine has being a Fortran call from C
			my $alreadyPresent = 0;
			foreach $maybeFortranCall (@maybeFortranCalls){
			    if  ($maybeFortranCall eq $1) { # $1 is without the ending "_"
				$alreadyPresent = 1;
			    }
			}
			if ($alreadyPresent == 0){
			    @maybeFortranCalls = ( @maybeFortranCalls, $1);
			}			
		    }
		}
	       

		if ($line =~ /([\w\d_]+)_\((.+)\)/){
			if (	($line =~ /^[\s\t]*\/\*/) || # only check /* on same line
				($line =~ /^[\s\t]*\/\//)
				) {
				# only a comment
			} else {
#				print "$Cfile: probable Fortran call to $1 with parameters $2\n";
				my $alreadyPresent = 0;
				foreach $maybeFortranCall (@maybeFortranCalls){
					if  ($maybeFortranCall eq $1) {
						$alreadyPresent = 1;
					}
				}
				if ($alreadyPresent == 0){
					@maybeFortranCalls = ( @maybeFortranCalls, $1);
				}
			}
		}
		
	}
	close MYFILE;
}

my @fortran_extensions =('fi','F','f90');
@fortranCalledFromC = ();
my $routineName = "undefined";  # the current FORTRAN subroutine being processed
my $CparameterString; # C-style parameter string (associative list)
my $CtypedParameterString; # for signatures; above is for calls

foreach $ext (@fortran_extensions){
	if ($os ne 'MSDOS') {
		@files = `ls *.$ext`;
	} else {
		# must emulate the above in case the script runs on Windows instead of Linux
		@files = ();
		@files_entries = `dir`; # dir *.$ext failed
		foreach $entry (@files_entries){
			chop $entry; # on DOS
			# we caught the problem: upon invoking dir,
			# one returns more than one file per line
			# well, actually not on DOS, only for Linux!
			if ($entry =~ /[\s\t]+([\w\d\-_]+\.$ext)/){
				@files = (@files, $1);
			}
		}
		# at this stage @files contains the list of all Fortran files
	}
    foreach $file (@files){
		chomp $file; # safer than chop: for DOS, no ending newline

    	open(MYFILE,"<$file");		
	my @args = (); # will be redefined each time we see a new subroutine with parameters

	# preprocessing: concatenate lines which are speparated by '& \n &'

	my $fileContents = "";
	my $simplifiedFile = "$file\_concatenated"; # to be deleted afterwards
	open MY_SIMPLIFIED_FILE, ">$simplifiedFile";
	while(<MYFILE>) { $fileContents .= $_; }
	$fileContents =~ s/&[\s\t]*\n[\s\t]*&//g;
	print MY_SIMPLIFIED_FILE $fileContents;
	close(MYFILE);
	close(MY_SIMPLIFIED_FILE);

	#reopen
	open MYFILE, "<$simplifiedFile";
	LINE_READ: while(<MYFILE>){
		my $line = $_;
		$retreiveArgType = 0; # a mode we enter after finding each subroutine
		
                # very special case: res_index in resindex.F relies on an include file resindex.fi 
		# for its parameters. Let's handle it manually:


		if ($line =~ /^[\s\t]*subroutine[\s\t]+res_index\(skew,mynorder,myn1,myn2,indexa,mynres\)/){
			if ($dbg==1){
		    		print "Specific treatment of res_index\n";
			}
		    # we know this routine is called from C
		    @fortranCalledFromC = (@fortranCalledFromC, "res_index");
		    $fortranFile{'res_index'}="resindex.F";
		    $CtypedParameterString{'res_index'} = "int* skew,int* mynorder,int* myn1,int* myn2,int indexa[4][1000],int* mynres";
		    $CparameterString{'res_index'} = "skew,mynorder,myn1,myn2,indexa,mynres";
		    $allTypesKnown{'res_index'} = 1;
		    $retreiveArgType = 0; # do not attempt to interpret following statements as param defs
		    next LINE_READ;
		}


		if ($line =~ /^[\s\t]*subroutine (.+)/) {
		#    if ($1 =~ /w_ptc_create_universe/) {print "UNIVERSE\n";exit;}
			@args=(); # default
			# entry should be unique	
		#	print "Fortran: found subroutine $1\n";
			$signature = $1;
			$signature =~ /(.+)\((.*)\)/;
			$routineName = $1;
			$fortranFile{$routineName}=$file;
			my $argString = $2;
			# clean-up spaces
			$argString =~ s/[\s\t]//g;

			my $calledFromC = 0;
			foreach $maybeFortranCall (@maybeFortranCalls){
				if ($maybeFortranCall eq $routineName) {
					$calledFromC = 1;
				}
			}

			if ($calledFromC == 1 ) {
				# should make sure the routine is not already here
				my $alreadyKnown = 0 ;
				foreach $entry (@fortranCalledFromC){
					if ($entry eq $routineName){
						if ( $dbg==1 ){ 
							print "routine $routineName already in list!\n";
						}
						$alreadyKnown = 1;
						# exit;
					}
				}
				# otherwise there will be overwrite (name clash?)

				if ($alreadyKnown==0){
					@fortranCalledFromC = (@fortranCalledFromC, $routineName);
				}			
				if ($argString eq ""){
					if ($dbg ==1 ) {
				    		print "routine $1 has NOARG\n";
					}
				    # routine has no arg
				    $allTypesKnown{$routineName}=1;
				    # now store the parameters C-string for this routine
				    $CparameterString{$routineName} = "";
				    $CtypedParameterString{$routineName} = "";
				    $retreiveArgTypeMode = 0;
				} else {

				    # now check the argument list for this subroutine
				    @args = split /,/, $argString;
				    if ($dbg==1) {
				    	print "FORTRAN routine $ routineName has args @args\n";
					}

				    if (scalar(@args)>0) {
					foreach $arg (@args) {
					    $type{$arg}='undefined';
					}
					$retreiveArgTypeMode = 1;
				    }
				    # now we'll try to retreive arguments one by one
				}
			}
			
			
		} # if /subroutine/
		else {
		if ($retreiveArgTypeMode == 1) { # if else should avoid entering here after above

			if ($line =~ /^[\s\t]*(real\(?[\w\d\(\)]*\)?|double precision|integer|character\*?\d*|complex|logical\(?[\w\d]*\)?)[\s\t]+(.+)/){
				my $aType = $1;
				my $argString = $2;
				
				# handle specific types such as character*16...
				if ($1 =~ /character\*?\d*/){ 
					$aType = 'character';
					if ($dbg == 1) { 
						print "CHARACTER in line '$line',with \$1=$aType, \$2=$argString=>" ;
					} 
				}

				if ($1 =~ /real\(kind\(1d0\)\)/){
					$aType = 'double precision';
				} else {
					if ($1 =~ /real\(?[\w\d\(\)]*\)?/) { $aType = 'real';}
				}
				if ($1 =~ /logical\(?[\w\d]*\)?/) { $aType="logical"}

				# in Fortran arrays appear as a(b) or a(b,c) where the type
				# of the array is retreived from a, hence we'll discard all parentheses
				$argString =~ s/\([\w_\d\s:\+\*,]+\)//g;
				# repeat above substitution to handle nested parentheses
				## $argString =~ s/\([\w_\d\s:\+\*,]+\)//g;
				$argString =~ s/[\s\t]//g; # remove blanks
				$argString =~ s/:://g; # sometimes we get type :: variable
				$argString =~ s/!(.+)//g; # remove Fortran-style comments
				# print "argstring is $argString\n";
				@someArgs = split /,/, $argString;
				
				if ( $dbg == 1 ) {
					print "COMPARE @args and @someArgs (from $argString)\n";
				}
 
				foreach $anArg (@someArgs){
					foreach $arg (@args){
						if ($anArg eq $arg){
							$type{$arg}= $aType;
							if ( $dbg ==1 ) {	
								print "$routineName: Type '$aType' found for arg '$arg'\n";
							}
							
						}
					}
				}
				# check wether we collected all args for the subroutine of interest
				$allTypesKnown{$routineName} = 1;
				foreach $arg (@args) {
					if ($type{$arg} ne 'undefined') {
						next;
					} else {
						$allTypesKnown{$routineName} = 0;
					}
				}
				if ($allTypesKnown{$routineName}==1) {
					# now store the parameters C-string for this routine
					my $paramStr = "";
					my $paramTypedStr = "";
					my $counter = 0;
					foreach $arg (@args){
						my $Ctype;
						if ($type{$arg} eq "integer") { $Ctype="int*"};
						if ($type{$arg} eq "real") { $Ctype="float*"};
						if ($type{$arg} eq "character") { $Ctype = "char*"};
						if ($type{$arg} eq "double precision") {$Ctype="double*"};
						if ($type{$arg} eq "complex") { $Ctype = "??"};
						if ($type{$arg} eq "logical") { $Ctype = "int*"};
						my $n = scalar(@args);
						if (($counter>0) && ($counter < $n)) {
							$paramStr .=",";
							$paramTypedStr .= ","; 
						}
						$counter++;
						$paramTypedStr = $paramTypedStr . "$Ctype $arg";
						$paramStr = $paramStr . "$arg";
						
					}
					$CparameterString{$routineName} = $paramStr;
					$CtypedParameterString{$routineName} = $paramTypedStr;
					$retreiveArgTypeMode = 0;
				}	# we know all the arg types for this FORTRAN routine
						
			}
		
		}	
		} # try to collect args
	}
	close MYFILE;
	
	if ($os ne 'MSDOS'){
		`rm $simplifiedFile`;	# delete temporary file (Unix or Linux)
	} else {
		`erase $simplifiedFile`; # delete temporary file (Windows)
	}
    }
}

my $nbRoutines = scalar(@fortranCalledFromC);

# now generate the wrapper files
my $wrapperFile = ">fortran_wrappers.c";
my $wrapperDefinesFile = ">fortran_wrappers.h";
my $fortranProtosFile = ">fortran_prototypes.h";
my $fortranWrappersProtosFile = ">fortran_wrappers_prototypes.h";
open WRAPPER_FILE, $wrapperFile;
open DEFINE_WRAPPER_FILE, $wrapperDefinesFile;
open FORTRAN_PROTOS_FILE, $fortranProtosFile;
open FORTRAN_WRAPPERS_PROTOS_FILE, $fortranWrappersProtosFile;

print WRAPPER_FILE "/* set of $nbRoutines wrappers to synchronize FORTRAN and C stdout buffers */\n";
print WRAPPER_FILE "/* when crossing the border upon calling FORTRAN from C. */\n\n";
print WRAPPER_FILE "#include <stdio.h>\n";
print WRAPPER_FILE "#include \"fortran_prototypes.h\"\n\n";
print WRAPPER_FILE "extern void call_fortran_flush_();\n\n";


print DEFINE_WRAPPER_FILE "\#ifndef _FORTRAN_WRAPPERS_H\n";
print DEFINE_WRAPPER_FILE "\#define _FORTRAN_WRAPPERS_H\n";
print DEFINE_WRAPPER_FILE "/* redirect FORTRAN calls to wrappers that synchronize FORTRAN and C stdout buffering */\n";
print DEFINE_WRAPPER_FILE "/* when crossing the border upon calling FORTRAN from C. */\n";
print DEFINE_WRAPPER_FILE "\#include \"fortran_wrappers_prototypes.h\"\n";

print FORTRAN_PROTOS_FILE "/* to avoid warnings of implicit declarations from fortran_wrappers.c */\n";
print FORTRAN_PROTOS_FILE "#ifndef _FORTRAN_PROTOTYPES_H\n";
print FORTRAN_PROTOS_FILE "#define _FORTRAN_PROTOTYPES_H\n";

print FORTRAN_WRAPPERS_PROTOS_FILE "/* to avoid warnings of implicit declaration from code calling the functions below */\n";
print FORTRAN_WRAPPERS_PROTOS_FILE "#ifndef _FORTRAN_WRAPPERS_PROTOTYPES_H\n";
print FORTRAN_WRAPPERS_PROTOS_FILE "#define _FORTRAN_WRAPPERS_PROTOTYPES_H\n";

foreach $subroutine (@fortranCalledFromC){
    if ($allTypesKnown{$subroutine}) {
	print DEFINE_WRAPPER_FILE "#define $subroutine\_ $subroutine\_wrapper\n";	
	my $filename = $fortranFile{$subroutine};                                            
	print WRAPPER_FILE "/* Wrap '$subroutine' defined in '$filename' */\n";
	print WRAPPER_FILE "void $subroutine\_wrapper(";
	print WRAPPER_FILE $CtypedParameterString{$subroutine};
	print WRAPPER_FILE "){\n";
	print WRAPPER_FILE "\tfflush(stdout);\n";
	print WRAPPER_FILE "\t$subroutine\_($CparameterString{$subroutine});\n";
	print WRAPPER_FILE "\tcall_fortran_flush_();\n";
	print WRAPPER_FILE "}\n";
       
	print FORTRAN_PROTOS_FILE "/* Wrap '$subroutine' defined in '$filename'*/\n";
	print FORTRAN_PROTOS_FILE "void $subroutine\_(";	
 	print FORTRAN_PROTOS_FILE $CtypedParameterString{$subroutine};
	print FORTRAN_PROTOS_FILE ");\n";
	
	print FORTRAN_WRAPPERS_PROTOS_FILE "/* Wrap '$subroutine' defined in '$filename'*/\n";	
	print FORTRAN_WRAPPERS_PROTOS_FILE "void $subroutine\_wrapper(";
	print FORTRAN_WRAPPERS_PROTOS_FILE $CtypedParameterString{$subroutine};
	print FORTRAN_WRAPPERS_PROTOS_FILE ");\n";
	
    } else {
	print WRAPPER_FILE "\n/* Skipped '$subroutine' defined in '$fortranFile{$subroutine}' due to incomplete type identification */\n\n";
    }
}
print DEFINE_WRAPPER_FILE "\#endif\n";
print FORTRAN_PROTOS_FILE "#endif\n";
print FORTRAN_WRAPPERS_PROTOS_FILE "#endif\n";

close WRAPPER_FILE;
close DEFINE_WRAPPER_FILE;
close FORTRAN_PROTOS_FILE;
close FORTRAN_WRAPPERS_PROTOS_FILE;
