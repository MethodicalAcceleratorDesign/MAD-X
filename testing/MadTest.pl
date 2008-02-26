#!/usr/bin/perl

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

$startTime = localtime;

$testReport = ""; # will be stored into an HTML document

$pwd = `pwd`;
chop $pwd;
$localRootDir = $pwd;

# checkout reference examples from the CVS
my $rmRes = `rm -rf ./madX-examples`;
my $checkoutRes = `cvs -d :kserver:isscvs.cern.ch:/local/reps/madx-examples checkout madX-examples`;

# $samplesRootDir = '/afs/cern.ch/user/f/frs/public_html/mad-X_examples';

$samplesRootDir = "$localRootDir/madX-examples/REF";


$htmlRootDir = '/afs/cern.ch/user/n/nougaret/www/mad';
$htmlFile = "$htmlRootDir/test.htm"; # for the time being



if ( $#ARGV != 0 ) {
    print "expect 1 argument: (1) MAD executable directory! EXIT!\n" ;
    exit ;
} else {
    $madDir = @ARGV[0];
    # check specified directory indeed contains executable
    $existsMadx = `ls $madDir/madx | wc -l`;
    if ($existsMadx == 0) {
	print "Madx missing in specified directory => exit!\n";
	exit;
    }
    $existsMadxp = `ls $madDir/madxp | wc -l`;
    if ($existsMadxp == 0){
	print "Madxp missing in specified directory => exit!\n";
	exit;
    }
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

    $madxLink = $madDir . "/madx";
    $madxpLink = $madDir . "/madxp";

    print "madxLink is $madxLink\n";
    print "madxpLink is $madxpLink\n";

}



@targetDirs = `ls $samplesRootDir`;


# search all test examples' directories
# (only a subset of them will be processed by the automated test, as specified in a separate XML document)

# in future version, should grow the call-graph from the list of source mad files provided in the XML file

# not instead of being performed beforehand on all possible files, this dependency analysis should only be
# carried-out on the actual target directories...

foreach $targetDir (@targetDirs) {

    chop $targetDir;

    # DBG
   # if (($targetDir ne "ptc_accel")&&($targetDir ne "ptc_madx_interface")) {next;} # only one target
    #if ($targetDir ne "twiss"){next;}

    print "target = '$targetDir'\n";

    chdir("$samplesRootDir/$targetDir");
    @subdirectories =`ls -d */`; # returns directories only (with '/' suffix)

    @alldirectories = ('./',@subdirectories);
    foreach $dir (@alldirectories) {
	chop $dir; # end-of-line
	chop $dir; # /
	$pwd = `pwd`; # for print below
	print "process directory $dir in $pwd\n";
	if ($dir eq "test") { next; } # Very specific case: twiss/test directory causes multiple key error in call-graph
	chdir("$samplesRootDir/$targetDir/$dir"); # go to subdir to open files with getListOfDependantFiles()
	@files = `ls`;
	foreach $file(@files) {
	    chop $file;
	    print "process $file"; 
	    $dependentFileList = getListOfDependantFiles($file);
	    # avoid name clash by prefixing
	    $fileKey = $targetDir . "/" . $file;
	    if ($dependencyList{$fileKey} ne "") {
		# actually no trouble if the entry is the same
		if ($dependencyList{$fileKey} eq $dependentFileList) {
		    print "WARNING: multiple entry $fileKey in call-graph\n";
		    $testReport .= "</p><font color=\"#FFCCBB\">WARNING: multiple entry $fileKey in call-graph</font></p>\n";

		} else {
		    $testReport .= "</p><font color=\"#FFBBBB\">ERROR: multiple incompatible entry $fileKey in call-graph</font></p>\n";
		    print "ERROR: multiple incompatible entry $fileKey in call-graph\n";
		    print "previous entry = '$dependencyList{$fileKey}'\n";
		    print "new entry = $dependentFileList\n";
		}
	    } # let's make sure
	    
	    $dependencyList{$fileKey}=$dependentFileList; # comma-separated list of files that depend on key
	    if ($dependencyList{$fileKey} ne "") {
		print " -> calls $dependencyList{$fileKey}\n";
	    } else {
		print " -> calls no other file\n";
	    }	    
	}
	chdir("$samplesRootDir/$targetDir/$dir");
    }
    chdir($localRootDir); # back to local root dir

}


$localTestDir = "$localRootDir/TESTING";
`rm -rf $localRootDir/TESTING`; # cleanup
mkdir($localTestDir, 0777);


@targets = `xsltproc --stringparam what list_targets ProcessScenario.xsl TestScenario.xml`; # all target functionalities

foreach $target (@targets) {
    chop $target;
    # DBG
   # if (($target ne "ptc_accel")&&($target ne "ptc_madx_interface")) {next; } # only one target
   # if ($target ne "twiss") {next;}

    print "--- testing $target\n";

    $testReport .= "<table width=\"75%\" border=\"0\">\n";
    $testReport .= "<tr class='test_target'><td colspan=\"2\"><div align=\"center\"><strong>$target</strong></div></td></tr>\n";

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
	print "command='$command'\n";	
	
	# retreive the subdirectory relocation if any (specific subdirs are specified in the XML)
	$sourceSubDir = ""; # by default
	$_ = $command;
	if(/subdirectory=([\w_\-.\d]+)/){
	    $sourceSubDir = $1;
	}

	# retrieve the 'input file' name with a regular expression
	$_ = $command;
	/<[\s\t]+([\w._\d]+)[\s\t]+>/;
	$infilename = $1;
	print "the input file name is: $infilename\n";
	
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

	$testReport .= "<tr class='test_case'><td width=\"70%\">$testCaseDir: $executableCommand</td><td width=\"30%\"></td></tr>\n"; 
        # above sets column width for the whole table


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
			    print "Found $secondLevelInput called by $testCaseDir/$input and to be copied locally\n";
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
	    print "Input is $input\n";
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
		# print "Found dependency file '$dependencyFile' under '$dependencyDir'\n";		
		$existsDir = `ls -d $dependencyDir | wc -l`;
		# print "existsDir now equals $existsDir\n";
		if ($existsDir == 1) {
		    print "dependency directory '$dependencyDir' already exists\n";
		    # simply copy the file
		    if  ($sourceSubDir eq "") {
			`cp $samplesRootDir/$target/$dependencyDir/$dependencyFile ./$dependencyDir`;
		    } else {
			`cp $samplesRootDir/$target/$sourceSubDir/$dependencyDir/$dependencyFile ./$dependencyDir`; # should merge the two with '.' for $sourceSubDir
		    }

		} else {
#		    print "dependency directory '$dependencyDir' will be created\n";
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
		    print "Second-level copy of $secondLevelInput\n";
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
			    print "dependency second directory '$dependencyDir/$secondDependencyDir' already exists\n";
			    # simply copy the file
			    if  ($sourceSubDir eq "") {
				`cp $samplesRootDir/$target/$dependencyDir/$secondDependencyDir/$secondDependencyFile ./$dependencyDir/$secondDependencyDir`;
			    } else {
				# current focus
				`cp $samplesRootDir/$target/$sourceSubDir/$dependencyDir/$secondDependencyDir/$secondDependencyFile ./$dependencyDir/$secondDependencyDir`; # should merge the two with '.' for $sourceSubDir
				print "now copying second-level $samplesRootDir/$target/$sourceSubDir/$dependencyDir/$secondDependencyDir/$secondDependencyFile into ./$dependencyDir/$secondDependencyDir\n";
			    }
			} else {
			    print "Not ready yet to handle creation of secondary dir/n";
			    $testReport .= "ERROR: Not ready yet to handle creation of second-level dir\n";
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
			mkdir($calledFileSubDir,0777);
			# add this dir to the list of input subdirectories that should be transferred from the workdir
			# in which the MAD command is executed into the "inputs" directory (otherwise the directory
			# would later go under "outputs" and undergo the side-by-side comparison...)
			@inputSubdirectories = ( @inputSubdirectories, $calledFileSubDir );

		    }
		    # and now do the copy
		    # print "DIRECTORY: now copying $samplesRootDir/$target/$input into ./$calledFileSubDir\n";
		    if ($sourceSubDir eq ""){
			`cp $samplesRootDir/$target/$input ./$calledFileSubDir`;
		    } else {
			`cp $samplesRootDir/$target/$sourceSubDir/$input ./$calledFileSubDir`;
		    }
		    

		#    @what = `ls ./$calledFileSubDir`;
		#   $howMany = scalar(@what);
		#    print "DIRECTORY $howMany contents of $pwd/$calledFileSubDir =\n";
		#    foreach $line (@what ) { print "$line\n"; }
		# --- above to handle the specific case of twiss lhc.madx


		} else {
		    # file to be called in the same directory as the input file
		    # print "for target '$target' and test input '$infilename', now copying additional '$input'\n";
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
	if ($1 eq "madxp") { $madLink = $madxpLink; $madProgram = "madxp"; } else {
	    if ($1 eq "madx") { $madLink = $madxLink; $madProgram = "madx"; } else {
		$madLink = "."; $madProgram = "linkUnknown";
	    }
	}

	`ln -s $madLink $madProgram`;

	`$executableCommand`;


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
	    print "--- for $testCaseDir ---\n";
	    print "list of input files: @inputs\n";
	    print "list of output files: @outputs\n";
	    print "------------------------\n";
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

	# now compare desired output and actual output

	# let's try to do a blind diff without trying to cure any 'standard' discrepancy
	$diffResFilename = "DIFFERENCES.txt";
	open (OUT,">$diffResFilename");
	foreach $file (@outputs) {

	    # specific case: 'temp' entry refers to a temporary file name
	    # (if other files, should be omitted in the same way...)
	    if ($file eq "temp") { next; }


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
		$testReport .="<tr class='omit'><td>$file</td><td>no file for reference</td></tr>\n";
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

		
		$testReport .= "<tr class='$madDiffRes'><td>$file</td><td><a href=\"$detailsLink\">$madDiffRes</a></td></tr>\n";
		print OUT "#COMPARING $file yields $madDiffRes\n";

	    }
	}
	close(OUT);	
      
    } # tests
    $testReport .= "</table>\n"; # end of table per target

} # targets

$endTime = localtime;

# create web page
my $html = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">';
$html .= "<html>\n";
$html .= "<head>\n";
$html .= "<title>MAD test result</title>\n";
$html .= "<link rel=stylesheet href=\"./MadTestWebStyle.css\" type=\"text/css\">"; # CSS stylesheet
$html .= "</head>\n";
$html .= "<!-- automatically generated by the MAD testing script -->\n";
$html .= "<body>\n";
$html .= "<p>Test started $startTime, ended $endTime</p>\n";
$html .= $testReport;
$html .= "</body>\n";
$html .= "</html>\n";
open(OUTHTML, ">$htmlFile");
print OUTHTML $html;
close OUTHTML;


# then send an e-mail to the MAD team

$msg = MIME::Lite->new(
		       From     => 'MAD.test.program@cern.ch',
		       To       => 'Jean-Luc.Nougaret@cern.ch',
		       Subject  => "Automated MAD Testing",
		       Data     => "This is an automated message.\nSee report on:\nhttp://test-mad-automation.web.cern.ch/test-mad-automation\n"
		       );
$msg->send;

print "script terminated\n";

sub getListOfDependantFiles {
    my $parentFilename = $_[0];
    my $dependentFileList = "";
    open(IN,$parentFilename);

    LINE: while(<IN>){
	# take into account MAD various ways of reading a file, namely using either "call,file" or "readmytable,file"
	my $fileRetreival = 0 ;
	my $child;
	# MAD syntax too permissive: hard to grep commands
	# not left-anchored as can be an action after an 'if'
	# but make sure it has not been commented (negated !) as is frequently the case...

	if (/^[\s\t]*!/) { next LINE; } # this is a comment (!) - do not bother

	if (/[\s\t]*[Rr][Ee][Aa][Dd][Mm][Yy][Tt][Aa][Bb][Ll][Ee],?[\s\t]*[Ff][Ii][Ll][Ee][\s\t]*=[\s\t]*[\"\']?([\w\._\-\d\/]+)[\"\']?[\s\t]*/){
	    $child = $1;
	    $fileRetreival = 1;

	}

	if (/[\s\t]*[Cc][Aa][Ll][Ll],?[\s\t]*[Ff][Ii][Ll][Ee][\s\t]*=[\s\t]*[\"\']?([\w\._\-\d\/]+)[\"\']?[\s\t]*;/) {
	    $child = $1;
	    $fileRetreival = 1;
	}
	
	# another - rare -  instruction that calls a file from another
	if (/[\s\t]*sxfread[\s\t]*,?[\s\t]*file[\s\t]*=[\s\t]*[\"\']?([\w\._\-\d\/]+)[\"\']?[\s\t]*;/) {
	    $child = $1;
	    $fileRetreival = 1;
	}



	if ($fileRetreival == 1) {
	#    $child = $1;
	    # before adding this child, make sure it's not already part of the childs the parent depends from
	    $_ = $dependentFileList;
	    if (/$child,/) {
		# print "'$child' already belonging to dependentFileList '$dependentFileList' => omit insertion\n";
	    } else {
		# print "add child '$child' to list '$dependentFileList'\n";
		$dependentFileList = $dependentFileList . $child . ",";
	    }
	}
    }
    close(IN);
    $_ = $dependentFileList ; # output arg
}


