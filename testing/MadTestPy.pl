#!/usr/bin/perl

# now supports a specific $debugMode to carry-out tests of a single target.
# sample invocation: MadTest.pl ./MadCvsExtract/madX debug=match
# In this case, HTML file becomes test_debug.html instead of test.htm.

# input: directory name in which madx and madxp are present
# output: a hierachy of directories to hold the tests' inputs and outputs
# at the same time an HTML document is created and moved to the web folder

# Still to do :
# (1) expand call-tree till the leafs
# (2) accomodate for different directory structures
# (3) allow for files loading from different directories in the call-tree
# (4) check potential troubles with namespace restricted to target name
# (5) ideally HTML formatting should be moved out of the code, relying on XML, XSLT and CSS instead

use MIME::Lite; # to send e-mail

open REPORT_FILE, ">MadTest_Report.txt";

$startTime = localtime;

my $path = $ENV{'PATH'};
my $newPath = $path . ":."; # can invoke local commands from MAD-X script
$ENV{'PATH'}=$newPath;



print REPORT_FILE "MadTest.pl report from $startTime\n";

$testReport = ""; # will be stored into an HTML document

$pwd = `pwd`;
chop $pwd;
$localRootDir = $pwd;

# for Makefile_develop and Makefile_nag, stderr is redirected
my $stderrFile = "stderr_redirected";

my $debugMode; 		# select a specific target and writes summary HTML file
			# to test_debug.htm instead of test.htm
my $debugTarget;	# meaningful iff $debugMode is set to 1

my $silentMode = 0;     # if set, no message will be sent. 0 by default
		
if ( $#ARGV < 0 ) {
    print "expect at least 1 argument: (1) MAD executable directory [(2) debug:<target>]  EXIT!\n" ;
    exit ;
} else {
    if (($#ARGV == 1) || ($#ARGV == 2)) {
	# second parameter for debug-mode
	my $option = @ARGV[1];
	if ($option =~ /^debug=(\w+)/){
	    $debugMode = 1;
	    $debugTarget = $1;
	    print "debug mode set for target '$debugTarget'\n";
	} elsif ($option =~ /^silent$/){
	    $silentMode = 1;
	} else {
	    print "the second optional argument '$option' should be of the form:";
	    print "debug=<target>\n";
	    print "EXIT!\n";
	    exit;
	}
    }
	
	$madDir = @ARGV[0];
	# check specified directory indeed contains executable
	$existsMadx = `ls $madDir/madx | wc -l`;
	if ($existsMadx == 0) {
		print REPORT_FILE "Madx missing in specified directory => exit!\n";
		exit;
	}
#	$existsMadxp = `ls $madDir/madxp | wc -l`;
#	if ($existsMadxp == 0){
#		print REPORT_FILE "Madxp missing in specified directory => exit!\n";
#		exit;
#	}
	# expand the full path
	$_ = $madDir;
	if (/^.\/([\w\d_\-\/.]+)/) {
		# local path specified
		$madDir = $localRootDir . "/" . $1;
	} else { 
		# full path already given 
	}

	# remove last "/" if any
	$_ = $madDir;
	if (/\/$/) { chop $madDir;  } # $ for end anchoring of the string


}


# checkout reference examples from the CVS
my $rmRes = `rm -rf ./madX-examples`;

if ($debugMode ==0){
	# default: we retrieve all the CVS
	print REPORT_FILE "extract the complete CVS repository containing examples\n";
	my $checkoutRes = `cvs -d :gserver:isscvs.cern.ch:/local/reps/madx-examples checkout madX-examples`;
	print REPORT_FILE $checkoutRes;
} else {
	print REPORT_FILE "extract only the CVS repository containing specific example";
	print REPORT_FILE " $debugTarget\n";
	# extract only the $debugTarget from the CVS
	my $checkoutRes = `cvs -d :gserver:isscvs.cern.ch:/local/reps/madx-examples checkout madX-examples/REF/$debugTarget`;
	print REPORT_FILE $checkoutRes;
}

# $samplesRootDir = '/afs/cern.ch/user/f/frs/public_html/mad-X_examples';

$samplesRootDir = "$localRootDir/madX-examples/REF";


$htmlRootDir = '/afs/cern.ch/user/n/nougaret/www/mad';

if ($debugMode == 0) {
	$htmlFile = "$htmlRootDir/test.htm"; # for the time being
} else {
	$htmlFile = "$htmlRootDir/test_debug.htm"; # for the time being
}
	
@targetDirs = `ls $samplesRootDir`;

# search all test examples' directories
# (only a subset of them will be processed by the automated test, as specified in a separate XML document)

# in future version, should grow the call-graph from the list of source mad files provided in the XML file

# not instead of being performed beforehand on all possible files, this dependency analysis should only be
# carried-out on the actual target directories...

# build the call-graph
TARGET_DIR: foreach $targetDir (@targetDirs) {

chop $targetDir;

if ($debugMode ==1) {
	if ($targetDir ne $debugTarget){
		next TARGET_DIR;
	}
}


print REPORT_FILE "target = '$targetDir'\n";

chdir("$samplesRootDir/$targetDir");
@subdirectories =`ls -d */`; # returns directories only (with '/' suffix)

@alldirectories = ('./',@subdirectories);
foreach $dir (@alldirectories) {
	chop $dir; # end-of-line
	chop $dir; # /
	$pwd = `pwd`; # for print below
	print REPORT_FILE "process directory $dir in $pwd\n";
	if ($dir eq "test") { next; } # Very specific case: twiss/test directory causes multiple key error in call-graph
	chdir("$samplesRootDir/$targetDir/$dir"); # go to subdir to open files with getListOfDependantFiles()
	@files = `ls`;
	foreach $file(@files) {
	chop $file;
	print REPORT_FILE "process $file"; 
	$dependentFileList = getListOfDependantFiles($file);
	# avoid name clash by prefixing
	$fileKey = $targetDir . "/" . $file;
	if ($dependencyList{$fileKey} ne "") {
		# actually no trouble if the entry is the same
		if ($dependencyList{$fileKey} eq $dependentFileList) {
		print REPORT_FILE "WARNING: multiple entry $fileKey in call-graph\n";
		$testReport .= "</p><font color=\"#FFCCBB\">WARNING: multiple entry $fileKey in call-graph</font></p>\n";

		} else {
		$testReport .= "</p><font color=\"#FFBBBB\">ERROR: multiple incompatible entry $fileKey in call-graph</font></p>\n";
		print REPORT_FILE "ERROR: multiple incompatible entry $fileKey in call-graph\n";
		print REPORT_FILE "previous entry = '$dependencyList{$fileKey}'\n";
		print REPORT_FILE "new entry = $dependentFileList\n";
		}
	} # let's make sure
	
	$dependencyList{$fileKey}=$dependentFileList; # comma-separated list of files that depend on key
	if ($dependencyList{$fileKey} ne "") {
		print REPORT_FILE " -> calls $dependencyList{$fileKey}\n";
	} else {
		print REPORT_FILE " -> calls no other file\n";
	}	    
	}
	chdir("$samplesRootDir/$targetDir/$dir");
}
chdir($localRootDir); # back to local root dir

}

@makefiles = ('Makefile_develop','Makefile_nag','Makefile'); # Makefile must be last
# because we need to have collected the results of the two previous makefiles in order
# to prepare the final web-page, which is only created for Makefile ...

# hierarchy to accomodate the local copy of the tests
`rm -rf $localRootDir/TESTING`; # cleanup
mkdir("$localRootDir/TESTING", 0777);

foreach $makefile (@makefiles) { # repeat for Makefile_develop, Makefile_nag, Makefile

# back to the initial directory
chdir $localRootDir;

my $regressionTest;
if ($makefile eq "Makefile") {
	$regressionTest = 1;
} else {
	$regressionTest = 0;
}


$localTestDir = "$localRootDir/TESTING/$makefile"; # separate work dirs for the 3 makefiles...
mkdir($localTestDir, 0777);

# following executables have been prepared by MadBuid.pl
$madxLink = $madDir . "/madx_$makefile";
$madxpLink = $madDir . "/madxp_$makefile";
print REPORT_FILE "madxLink is $madxLink\n";
print REPORT_FILE "madxpLink is $madxpLink\n";



@targets = `xsltproc --stringparam what list_targets ProcessScenario.xsl TestScenario.xml`; # all target functionalities

TARGET: foreach $target (@targets) {
chop $target;

$outcome{$target} = "success"; # by default, nobody will receive an e-mail about this target

if ($debugMode ==1){
	if ($target ne $debugTarget) {
		next TARGET;
	}
}

print REPORT_FILE "--- testing $target\n";

	if ($regressionTest) {
	$testReport .= "<table width=\"75%\" border=\"0\">\n";
	$testReport .= "<tr class='test_target'><td colspan=\"2\"><div align=\"center\"><strong>$target</strong></div></td></tr>\n";
	}

chdir($localRootDir); # top of the hierarchy

$targetDir = "$localTestDir/$target";
mkdir($targetDir, 0777) or die "fail to create directory $targetDir\n";

@tests = `xsltproc --stringparam what list_tests --stringparam target $target ProcessScenario.xsl TestScenario.xml`;


chdir("$localTestDir/$target") or die "fail to chdir to $localTestDir/$target\n"; # after processing stylessheets

# populate the local test directory with the input file as well as with all the files it depends from

$autonumber = 0;
foreach $test (@tests) {

	$autonumber++;
	chop $test;
	$command = $test;
	print REPORT_FILE "command='$command'\n";	
	
	# retreive the subdirectory relocation if any (specific subdirs are specified in the XML)
	$sourceSubDir = ""; # by default
	$_ = $command;
	if(/subdirectory=([\w_\-\.\d]+)/){
	$sourceSubDir = $1;
	}

	# retrieve the 'input file' name with a regular expression
	$_ = $command;
	/<[\s\t]+([\w\.\-_\d]+)[\s\t]+>/;
	$infilename = $1;
	print REPORT_FILE "the input file name is: $infilename\n";
	
	if ($sourceSubDir eq "" ) {
	$testCaseDir = "test_" . $autonumber;
	} else {
	$testCaseDir = "test_" . $autonumber . "_" . $sourceSubDir; 
	# keep the same name for the local test dir under 'target' as in the reference
	# and the files will need to be copied from location with the $testSubDir prefix later-on...
	}

	# before executing the command, make sure we remove
	# the optional subdirectory information which shows up after the comma
	# in case the source test directory contains a subdirectory structure...
	my $executableCommand;
	if ($sourceSubDir eq "") {
	$executableCommand = $command;
	} else {
	$_ = $command;
	s/,[\s\t]*[\w\d.\-_=]+//g;
	$executableCommand = $_;
	}

	if ($regressionTest) {
		my $key = $target . "_" . $testCaseDir;
		$testReport .= "<tr class='test_case'><td width=\"80%\">$testCaseDir: $executableCommand</td><td width=\"20%\"><table width=\"100%\" style=\"text-align: center\"><tr>$dev{$key} $nag{$key}</tr></table></td></tr>\n"; 
		# above sets column width for the whole table
	}


	mkdir("$localTestDir/$target/$testCaseDir",0777) or die "fail to create directory $testCaseDir\n";
	chdir("$localTestDir/$target/$testCaseDir") or die "fail to chdir to $localTestDir/$target/$testCaseDir\n";

	# copy the input file that corresponds to this specific test case
	if ($sourceSubDir eq "") {
		`cp $samplesRootDir/$target/$infilename .`;
	} else {
		`cp $samplesRootDir/$target/$sourceSubDir/$infilename .`;
	}
	
	$key = "$target/$infilename";

	# now copy additional input files for the test, according to the dependency information retreived above
	@inputs = split /,/, $dependencyList{"$target/$infilename"}; # prefixed with $target to avoid name clashes

	# now grow the @input list by expanding the dependency and by adding the root node

	# partial treatment: only look for files calling files under the same directory
	my $reccursionLevels =2;

	for ($i=0;$i<$reccursionLevels;$i++){
	foreach $input (@inputs){
		if ($dependencyList{"$target/$input"} ne "") 
		{
		my @secondLevelInputs = split /,/, $dependencyList{"$target/$input"};
		# or third etc... level according to reccursion level...
		foreach $secondLevelInput (@secondLevelInputs){
			if (/..\/([\w\d\-_\.\/]+)/) { 
			# file located in a directory with path starting above
			# currently handled differently by code appearing downthere
			# => do nothing for time-being
			# (incomplete: should also handle any tree structure)
			} else {
			print REPORT_FILE "Found $secondLevelInput called by $testCaseDir/$input and to be copied locally\n";
			# actually, should only add to the list if not already present
			my $existsInput = 0;
			INPUT_LOOP: foreach $existingInput (@inputs) {
				if ($secondLevelInput eq $existingInput) {
				$existsInput = 1;
				next INPUT_LOOP; 
				}
			}
			if ($existsInput ==0) {
				@inputs = (@inputs, $secondLevelInput);
			}
			} 
		}
		}
	} # grow list of inputs at one level
	}  # end growing the list of inputs that can be moved to the same input directory

	@inputs = ($infilename, @inputs); # add the root inputfile


	my @inputSubdirectories = (); # list of input subdirectories under the workdir in which MAD is invoked
	# these subdirectories will later be moved under "inputs" instead of "outputs". In most cases this
	# list is empty, but not for targets 'error' and 'twiss' for instance.


	# copyping inputs and dependent files locally
	foreach $input (@inputs) {
	print REPORT_FILE "Input is $input\n";
	$_ = $input;
	# SPECIFIC CASE: files that must be stored in a locally replicated hierarchy
	# considering the MAD call instructions mention a relative path...
	# note later-on should also accomodate for situations were the included files
	# are below, down the hierarchy... At this stage should re-implement path handling
	if (/\.\.\/([\w\d\-_\.\/]+)/) { # ../dir/file or ../file formats
		my $term = $1;
		if (/(\.\.\/[\w\d\-_]+)\/([\w\d-_\.]+)/) { # ../dir/file format only
		# file to be called is located up the hierarchy
		# in which case we need to reflect the directory tree structure
		# by creating directories if necessary
		$dependencyDir = $1;
		$dependencyFile = $2;
		} else {
		$dependencyDir = '..';
		$dependencyFile = $term;
		}
		# print REPORT_FILE "Found dependency file '$dependencyFile' under '$dependencyDir'\n";		
		$existsDir = `ls -d $dependencyDir | wc -l`;
		# print REPORT_FILE "existsDir now equals $existsDir\n";
		if ($existsDir == 1) {
		print REPORT_FILE "dependency directory '$dependencyDir' already exists\n";
		# simply copy the file
		if  ($sourceSubDir eq "") {	
			`cp $samplesRootDir/$target/$dependencyDir/$dependencyFile ./$dependencyDir`;
		} else {
			`cp $samplesRootDir/$target/$sourceSubDir/$dependencyDir/$dependencyFile ./$dependencyDir`; # should merge the two with '.' for $sourceSubDir
		}

		} else {
#		    print REPORT_FILE "dependency directory '$dependencyDir' will be created\n";
		# first create the directory, then copy the file
		mkdir($dependencyDir,0777);
		if ($sourceSubDir eq "") {
			`cp $samplesRootDir/$target/$dependencyDir/$dependencyFile ./$dependencyDir`;
		} else {
			`cp $samplesRootDir/$target/$sourceSubDir/$dependencyDir/$dependencyFile ./$dependencyDir`; # should merge the two with '.' for $sourceSubDir

		}
		}


		@secondLevelInputs = split /,/, $dependencyList{"$target/$dependencyFile"}; 
		# prefixed with $target to avoid name clashes (may be too coarse!)
		foreach $secondLevelInput (@secondLevelInputs) {
		print REPORT_FILE "Second-level copy of $secondLevelInput\n";
		$_ = $secondLevelInput;
		if (/\.\.\/([\w\d\-_\.\/]+)/) { # ../dir/file or ../file formats
			my $term;
			if (/(..\/[\w\d\-_]+)\/([\w\d-_\.]+)/) { # ../dir/file format only
			$secondDependencyDir = $1;
			$secondDependencyFile = $2;
			} else {
			$secondDependencyDir = '..';
			$secondDependencyFile = my $term;
			}

			$existsSecondDir = `ls -d $dependencyDir/$secondDependencyDir | wc -l`;
			if ($existsSecondDir) {
			print REPORT_FILE "dependency second directory '$dependencyDir/$secondDependencyDir' already exists\n";
			# simply copy the file
			if  ($sourceSubDir eq "") {
				`cp $samplesRootDir/$target/$dependencyDir/$secondDependencyDir/$secondDependencyFile ./$dependencyDir/$secondDependencyDir`;
			} else {
				# current focus
				`cp $samplesRootDir/$target/$sourceSubDir/$dependencyDir/$secondDependencyDir/$secondDependencyFile ./$dependencyDir/$secondDependencyDir`; # should merge the two with '.' for $sourceSubDir
				print REPORT_FILE "now copying second-level $samplesRootDir/$target/$sourceSubDir/$dependencyDir/$secondDependencyDir/$secondDependencyFile into ./$dependencyDir/$secondDependencyDir\n";
			}
			} else {
			print REPORT_FILE "Not ready yet to handle creation of secondary dir/n";
			$testReport .= "ERROR: ($makefile) Not ready yet to handle creation of second-level dir\n";
			}
		} else {
			# the second-level related file is located under the same directory
			if ($sourceSubDir eq "") {
			`cp $samplesRootDir/$target/$dependencyDir/secondLevelInput ./$dependencyDir`;
			} else {
			`cp $samplesRootDir/$target/$sourceSubDir/$dependencyDir/$secondLevelInput ./$dependencyDir`;
			}		    }

		} # second level
	
	} 
	else {

		# --- below to handle the specific case of twiss lhc.madx -> temp contents seem to be cleaned-up from .mad files...
		if ($input =~ /([\w\-_\d\.]+)\/([\w\d\-_\.]+)/) {

		# specific case of a file being called in a subdir of current directory
		my $calledFileSubDir = $1;
		my $calledFile = $2;
		# $calledFileSubDir = 'temp'; # FORCE FOR DBG
		# with the above, for lhc.out of twiss target, we indeed 
		# find the .mad files under 'TEMP', but if we keep the name they disappear
		# check if the subdir already exists - if not create it
		my $existsDir = `ls -d $calledFileSubDir | wc -l`;

		my $pwd = `pwd`; chop $pwd;
		if ($existsDir==1) {
		} else {
			# DBG: where we end-up for specific case of MB.12.mad
			mkdir($calledFileSubDir,0777);
			# add this dir to the list of input subdirectories that should be transferred from the workdir
			# in which the MAD command is executed into the "inputs" directory (otherwise the directory
			# would later go under "outputs" and undergo the side-by-side comparison...)
			@inputSubdirectories = ( @inputSubdirectories, $calledFileSubDir );

		}
		# and now do the copy
		# DBG: we end-up copying MB.12.mad into ./temp
		# print REPORT_FILE "DIRECTORY: now copying $samplesRootDir/$target/$input into ./$calledFileSubDir\n";
		if ($sourceSubDir eq ""){
			`cp $samplesRootDir/$target/$input ./$calledFileSubDir`;
		} else {
			`cp $samplesRootDir/$target/$sourceSubDir/$input ./$calledFileSubDir`;
		}
		

		#    @what = `ls ./$calledFileSubDir`;
		#   $howMany = scalar(@what);
		#    print REPORT_FILE "DIRECTORY $howMany contents of $pwd/$calledFileSubDir =\n";
		#    foreach $line (@what ) { print REPORT_FILE "$line\n"; }
		# --- above to handle the specific case of twiss lhc.madx


		} else {
		# file to be called in the same directory as the input file
		# print REPORT_FILE "for target '$target' and test input '$infilename', now copying additional '$input'\n";
		if ($sourceSubDir eq "") {
			`cp $samplesRootDir/$target/$input .`;
		} else {
			`cp $samplesRootDir/$target/$sourceSubDir/$input .`;
		}
		} # file to be called located in the same directory as the command's input file
	}

	}



	# check whether we should call madx or madxp
	my $madLink;
	$_ = $command;
	/.\/(madxp?)[\s\t]*/;
	if ($1 eq "madxp") { 
		$madLink = $madxpLink; 
		$madProgram = "madxp"; 
	} else {
		if ($1 eq "madx") { 
			$madLink = $madxLink; 
			$madProgram = "madx"; 
		} else {
			$madLink = "."; 
			$madProgram = "linkUnknown";
		}
	}

	`ln -s $madLink $madProgram`;

	# DBG: for twiss - test 5, MAD.12.mad etc in temp are replaced by MAD.12 etc...
	# hence temp is removed and rewritten whereas this does not seem to be the case
	# in the reference directory ... rights?
	
	if ($regressionTest ==0) { 
		# For Makefile_develop and Makefile_nag, we want to see stderr
		# => also redirect the standard error
		# by replacing:	'madx(p) < file.madx > file.out'
		# with:		'(madx(p) <file.madx > file.out) >& stderr_redirected'
		my $redirectedExecutableCommand = "($executableCommand) >& $stderrFile";
		$executableCommand = $redirectedExecutableCommand;
	}
	
	$testCaseStartTime = localtime;
	# run the test
	`$executableCommand`; 
	$testCaseEndTime = localtime;

	# retrieve the 'output file' name with a regular expression
	$_ = $command;
	/[\s\t]([\w._\d]*).out/;
	$outfilename = $1 . ".out";


	# list all by-product output files
	@allFilesNow = `ls`; # list of all input + output files after invoking 'mad' command 
	
	# remove the madx link from the list of files to be moved

	# grow the list of input files with the list of subdirectory so that they move together under 'inputs'
	foreach $dir (@inputSubdirectories){
	push (@inputs, $dir);
	}

	@outputs = ();
	foreach $file (@allFilesNow){
	chop $file;
	$isInput = 0;
	foreach $input (@inputs){
		if ($file eq $input) {
		$isInput =1;
		}
	}
	if (($isInput == 0) && ($file ne "madx") && ($file ne "madxp")) {        
		# ignore the madx/madxp entries which should stay on top
		push (@outputs, $file);
	}
	} ;

	#DBG
	my $dbg =1 ;
	if ($dbg==1){
	print REPORT_FILE "--- for $testCaseDir ---\n";
	print REPORT_FILE "list of input files: @inputs\n";
	print REPORT_FILE "list of output files: @outputs\n";
	print REPORT_FILE "------------------------\n";
	}

	# now move inputs and outputs into dedicated subdirectory
	$inputSubdir = "$localTestDir/$target/$testCaseDir/input";
	$outputSubdir = "$localTestDir/$target/$testCaseDir/output";
	mkdir($inputSubdir, 0777) or die "fail to create directory $inputSubdir\n";
	mkdir($outputSubdir, 0777) or die "fail to create directory $outputSubdir\n";	

	foreach $file (@inputs) { 
	# SPECIFIC CASE: files that must be stored in the locally stored hierarchy
	# with MAD call instructions referring to a relative path...
	$_ = $file;
	if (/..\//) { next; #skip. Later on should also deal with the case of directories under the test dir... 
		} else {
			`mv $file $inputSubdir/`; 			  
		}
	}
	foreach $file (@outputs) {`mv $file $outputSubdir/`; }


	chdir($outputSubdir) or die "fail to chdir to $outputSubdir\n";

	if ($regressionTest) { # i.e. makefile == Makefile_develop
		# now compare desired output and actual output

		# let's try to do a blind diff without trying to cure any 'standard' discrepancy
		$diffResFilename = "DIFFERENCES.txt";
		open (OUT,">$diffResFilename");
		foreach $file (@outputs) {

		# specific case: 'temp' entry refers to a temporary file name
		# (if other files, should be omitted in the same way...)
		if ($file eq "temp") { next; }

		# specific case: skip Postcript '.ps' and '.eps' files
		# (eventually Poscript files should be removed from the examples' reference CVS repo)
		if ($file =~ /\.ps$/) { next; }
		if ($file =~ /\.eps$/) { next; }

		# check there is always two files to be compared
		if ($sourceSubDir  eq "") {
			$fileCount = `ls $samplesRootDir/$target/$file | wc -l`;
			chop $fileCount;
		} else {
			$fileCount =  `ls $samplesRootDir/$target/$sourceSubDir/$file | wc -l`;
			chop $fileCount;
		}
		if ($fileCount == 0) {
			print OUT "# FAIL TO COMPARE $file: no such file for reference => FAILURE\n";
			$testReport .="<tr class='omit'><td width=\"70%\">$file</td><td width=\"30%\">no file for reference</td></tr>\n";
		} else {

			my $detailsLink;
			my $detailsHtmlFile;

			# specific case when the HTML file name is of the form XX.map or XX.map.htm 
			# webserver will fail to display the HTML although one can open it form the webfolder...
			# to overcome this limitation, we need to juggle with the HTML file name
			$_ = $file;
			s/\.map$/\.maAap/g;
			$f = $_;
		

			# handle both cases where there is a $sourceSubDir or not...
			$detailsLink = "./details/DiffResult_"."$target"."_"."$testCaseDir"."_"."$f.htm"; # weblink
			$detailsHtmlFile = "$htmlRootDir/details/DiffResult_"."$target"."_"."$testCaseDir"."_"."$f.htm"; # deliver		
		
			if ($sourceSubDir eq "") {
			$madDiffRes = 
				`$localRootDir/MadDiff.pl ./$file $samplesRootDir/$target/$file $detailsHtmlFile`;
			} else {
			$madDiffRes =
				`$localRootDir/MadDiff.pl ./$file $samplesRootDir/$target/$sourceSubDir/$file $detailsHtmlFile`;
			}

		
			$testReport .= "<tr class='$madDiffRes'><td width=\"70%\">$file</td><td width=\"30%\"><a href=\"$detailsLink\">$madDiffRes</a></td></tr>\n";
			print OUT "#COMPARING $file yields $madDiffRes\n";
			
			# summarize the information in $outcome{$target}
			if ($madDiffRes eq 'failure'){
				$outcome{$target}='failure';
			} else {
				if ($madDiffRes eq 'warning'){
					if ($outcome{$target} eq 'success'){
						$outcome{$target}='warning';
					}
					# otherwise, $outcome{$target} shall keep its value
					# whether it is 'success' or 'failure'
				}
			}

		}
		}
		close(OUT);	
	} else {
		# what we do to process the outcome of tests
		# compiled with Makefile_develop and Makefile_nag
		
		# info can be either (1) per-file : does it exist or not
		# or (2) per test-case i.e. [a] do we reach completion in the main output file
		# and if not, and [b] what is the stderr if any.
	
		# let's go for (2), i.e. per $testCaseDir
		
		# (1) try to open the output file, according to it show the 
		# cell as green / orange / red for completed / incomplete / absent
		my $existsOutput = 0;
		foreach $file (@outputs){
			if ($file eq $outfilename) {
				$existsOutput = 1;
			}
		}
		
		my $finishedNormally = 0;
		if ($existsOutput) {
			# now check whether the output completed normally
			open OUTFILE, "<$outputSubdir/$outfilename"; # output just moved
			while (<OUTFILE>){
				if (/\+[\s\t]+MAD\-X[\s\t]+([\d\.]+)[\s\t]+finished normally[\s\t]+\+/){
					$finishedNormally = 1;
				}
			}
			close OUTFILE; 
		}
		
		my $existsStderr = 0;
		@stderrFiles = `ls $outputSubdir/$stderrFile`;
		if (scalar(@stderrFiles) >0) {
			# check if not empty
			my $linesCount = `wc -l $stderrFile`;
			if ($linesCount >0) {
				$existsStderr=1;
			}
		}
		
		my $outputStatus; # used by CSS to apply font-color on stderr web page
		
		if ($existsOutput ==0 ){
			$outputStatus = 'failure';
		} else {
			if ($finishedNormally) {
				if ($existsStderr) {
					$outputStatus = 'warning';
				} else {
					$outputStatus = 'success'; # no link to stderr
				}
			} else {
				$outputStatus = 'failure';
			}
		}
		
		# update $outcome{$target}
		# summarize the information in $outcome{$target}
		if ($outputStatus eq 'failure'){
			$outcome{$target}='failure';
		} else {
			if ($outputStatus eq 'warning'){
				if ($outcome{$target} eq 'success'){
					$outcome{$target}='warning';
				}
				# otherwise, $outcome{$target} shall keep its value
				# whether it is 'success' or 'failure'
			}
		}		
		
		
		
				
		# (2) if there is an stderr file, then create an HTML file
		# and draw a link from the main web page

		# HTML files and links to make stderr viewable from the web
		my $m;
		if ($makefile eq "Makefile_develop") {
			$m="dev";
		}
		if ($makefile eq "Makefile_nag"){
			$m="nag";
		}
			
		# weblink:
		my $errorLink = "<a href=\"./details/Error_".$m."_".$target."_".$testCaseDir.".htm"."\">$m</a>";
		# deliver: 
		my $errorHtmlFile = "$htmlRootDir/details/Error_".$m."_".$target."_".$testCaseDir.".htm"; 
		`touch $errorHtmlFile`; # for the time-being
		
		
		errorWebPage("$outputSubdir/$stderrFile",
			$errorHtmlFile,
			$outputStatus,
			$testCaseStartTime,
			$testCaseEndTime);
		
		my $key = $target . "_" . $testCaseDir;
		if ($makefile eq "Makefile_develop") {
			$dev{$key}="<td class=\"$outputStatus\">$errorLink</td>";
		}
		if ($makefile eq "Makefile_nag"){
			$nag{$key}="<td class=\"$outputStatus\">$errorLink</td>";
		}


			
	}
} # tests

	if ($regressionTest) {
	$testReport .= "</table>\n"; # end of table per target
	}
} # targets
} # for each Makefile_develop, Makefile_nag, Makefile ...

$endTime = localtime;

# we know that at this stage $makefile equals "Makefile" anyway
# i.e. that we are finally performing regression testing

# create web page
my $html = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">';
$html .= "<html>\n";
$html .= "<head>\n";
$html .= "<title>MAD test result</title>\n";
$html .= "<link rel=stylesheet href=\"./MadTestWebStyle.css\" type=\"text/css\">"; 
# CSS stylesheet
$html .= "</head>\n";
$html .= "<!-- automatically generated by the MAD testing script -->\n";
$html .= "<body>\n";
$html .= "<p>Test started $startTime, ended $endTime</p>\n";
$html .= "\n<p>To update example X: \n";
$html .= "'cvs -d:gserver:isscvs.cern.ch/local/reps/madx-examples checkout madX-examples/REF/X' ";
$html .= "then modify input files, re-run and commit outputs to the CVS.</p>\n";
$html .= $testReport;
$html .= "</body>\n";
$html .= "</html>\n";
open(OUTHTML, ">$htmlFile");
print OUTHTML $html;
close OUTHTML;

# back to the initial directory
chdir $localRootDir;

# then send e-mails to the people responsible for the various tests' targets
if (($debugMode ==0) && ($silentMode==0) ) { # for time being, do not notify anybody
	foreach $target (@targets){
		# do we need to notify anyone for this target?
		if (($outcome{$target} eq 'failure')||($outcome{$target} eq 'warning')){
			# retreive the responsible person for each test target
			$responsible{$target} = `xsltproc --stringparam what responsible --stringparam target $target ProcessScenario.xsl TestScenario.xml`;	
			# retreive e-mail address of the responsible person
			my $emailRecipient = $responsible{$target};
			my $emailSubject = "MAD Unsucessful test for '$target'";
			my $emailContent .= "The test that started on $startTime and completed on $endTime resulted in a $outcome{$target}.\n";
			$emailContent .= "See detailed report on:\n";
			$emailContent .= "http://test-mad-automation.web.cern.ch/test-mad-automation\n";
			$emailContent .= "\nThis e-mail has been sent to you because you are registered as the responsible person for the '$target' package\n";
			my $msg = MIME::Lite->new(
				From	   =>	'Jean-Luc.Nougaret@cern.ch',
				'Reply-To' =>	'mad-automation-admin@cern.ch',
				To	   =>	$emailRecipient,
				Cc	   =>	'mad-automation-admin@cern.ch',
				Subject	   =>	$emailSubject,
				Data	   =>	$emailContent
			);
			$msg->send;
			print REPORT_FILE "sent e-mail to $emailRecipient\n";
		} else {
			# for this targets, all tests were sucessful => no need to notify
		}
	}
}

# for debug, inform me that the test prog ran till completion
my $lastMsg = MIME::Lite->new(
			      From	   =>	'Jean-Luc.Nougaret@cern.ch',
			      To	   =>	'Jean-Luc.Nougaret@cern.ch',
			      Subject	   =>	'Last message from MadTestPy.pl',
			      Data	   =>	'MadTestPy.pl ran succesfully till completion (not killed)'
			      );
$lastMsg->send;


print REPORT_FILE "script terminated\n";
close REPORT_FILE;

sub getListOfDependantFiles {
my $parentFilename = $_[0];
my $dependentFileList = "";
open(IN,$parentFilename);

LINE: while(<IN>){
	# take into account MAD various ways of reading a file, namely using either "call,file" or "readmytable,file"
	my $fileRetreival = 0 ;
	my $childCount = 0;
	my @childs;
	# MAD syntax too permissive: hard to grep commands
	# not left-anchored as can be an action after an 'if'
	# but make sure it has not been commented (negated !) as is frequently the case...

	if (/^[\s\t]*!/) { next LINE; } # this is a comment (!) - do not bother

	if (/[\s\t]*[Rr][Ee][Aa][Dd][Mm][Yy][Tt][Aa][Bb][Ll][Ee],?[\s\t]*[Ff][Ii][Ll][Ee][\s\t]*=[\s\t]*[\"\']?([\w\._\-\d\/]+)[\"\']?[\s\t]*/){
	@childs[$childCount++] = $1;
	$fileRetreival = 1;

	}

	if (/[\s\t]*[Cc][Aa][Ll][Ll],?[\s\t]*[Ff][Ii][Ll][Ee][\s\t]*=[\s\t]*[\"\']?([\w\._\-\d\/]+)[\"\']?[\s\t]*;/) {
	@childs[$childCount++] = $1;
	$fileRetreival = 1;
	}
	
	# another way to read a table
	if (/[\s\t]*[Rr][Ee][Aa][Dd][Tt][Aa][Bb][Ll][Ee],?[\s\t]*[Ff][Ii][Ll][Ee][\s\t]*=[\s\t]*[\"\']?([\w\._\-\d\/]+)[\"\']?[\s\t]*;/){
	@childs[$childCount++] = $1;
	$fileRetreival = 1;	    
	}

	# another - rare - instruction that calls a script available in the input directory
	# along the other input files. Such script must be copied locally along the other
	# files.
	# this is for instance the case for the read.magnet.errors perl-script
	# in twiss/test_5/ or foot
	if (/[Ss][Yy][Ss][Tt][Ee][Mm][\s\t]*,?[\s\t]*[\"\']([\w\._\-\d]+)[\s\t]*<?[\s\t]*([\w\.\_\-\d\/]*)[\s\t]*>?(.*)[\"\']/){
		my $cmd=$1;
		my $arg=$2;
		
		if (($cmd =~/mkdir/)||($cmd =~/ls/)||($cmd =~ /cat/)||($cmd =~/rm/)||($cmd =~ /grep/) ||($cmd =~ /echo/) || ($cmd =~ /ln/) || ($cmd =~ /cp/)){
			# command invocation, not implying a file call
		
		} else {
		
		    # ... many other cases like 'mkdir' to be handled before we conclude
		    # the command corresponds to a file...

		    if (($cmd ne "perl") && ($cmd ne "gnuplot") && ($cmd ne "python")){
			@childs[$childCount++] = $cmd; # the command ...
		    }
		    if ($arg ne "") {
			@childs[$childCount++] = $arg; # as well as its argument if any
			# print "found command '$cmd', with input argument '$arg'\n";
		    }

		    $fileRetreival = 1;
		}
	    }
	# NOTE: the above is FRAGILE. In many cases, one might expect the System call
	# to start with a command not corresponding to a local file in the input dir.
	
	
	# another - rare -  instruction that calls a file from another
	if (/[\s\t]*sxfread[\s\t]*,?[\s\t]*file[\s\t]*=[\s\t]*[\"\']?([\w\._\-\d\/]+)[\"\']?[\s\t]*;/) {
	@childs[$childCount++] = $1;
	$fileRetreival = 1;
	}


	if ($fileRetreival == 1) {
	#    $child = $1;
	# before adding this child, make sure it's not already part of the childs the parent depends from
	foreach $child (@childs) {
		if ($dependentFileList =~ /$child,/) {
			# print REPORT_FILE "'$child' already belonging to dependentFileList '$dependentFileList' => omit insertion\n";
		} else {
			# print REPORT_FILE "add child '$child' to list '$dependentFileList'\n";
			$dependentFileList = $dependentFileList . $child . ",";
		}
	}
	}
}
close(IN);
$_ = $dependentFileList ; # output arg
}

sub errorWebPage {
	my $errorFile = $_[0];
	my $htmlFile = $_[1];
	my $outputStatus = $_[2];
	my $startTime = $_[3];
	my $endTime = $_[4];
	my $contents ="";
	
	$contents .= "<table width=\"75%\" border=\"0\">\n";
	$contents .= "<p>Test started $startTime, ended $endTime.<p>\n";	
	$contents .= "<tr class=$outputStatus><td width=\"80%\">Contents of stderr</td><td width=\"20%\">$outputStatus</td></tr>\n";

	open(INERROR, "<$errorFile");
	while (<INERROR>){
		chop;
		$contents .= "<tr><td colspan=\"2\">$_</td><tr>\n";
	}
	$contents .= "</table>\n";
	close(INERROR);
	
	# create web page
	my $html = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">';
	$html .= "<html>\n";
	$html .= "<head>\n";
	$html .= "<title>MAD stderr file</title>\n";
	$html .= "<link rel=stylesheet href=\"../MadTestWebStyle.css\" type=\"text/css\">"; 
	# CSS stylesheet
	$html .= "</head>\n";
	$html .= "<!-- automatically generated by the MAD testing script -->\n";
	$html .= "<body>\n";
	$html .= $contents;
	$html .= "</body>\n";
	$html .= "</html>\n";
	open(OUTHTML, ">$htmlFile");
	print OUTHTML $html;
	close OUTHTML;
}

