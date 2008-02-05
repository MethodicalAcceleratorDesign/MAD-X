#!/usr/bin/perl

# input: directory name in which madx and madxp are present
# output: a hierachy of directories to hold the tests' inputs and outputs
# at the same time an HTML document is created and moved to the web folder

# Still to do :
# (1) expand call-tree till the leafs
# (2) accomodate for different directory structures
# (3) Calling madfile with symbolic name?
# (4) allow for files loading from different directories in the call-tree
# (5) check potential troubles with namespace restricted to target name
# (6) ideally HTML formatting should be moved out of the code, relying on XML, XSLT and CSS instead

use MIME::Lite; # to send e-mail

$startTime = localtime;

$testReport = ""; # will be stored into an HTML document

$samplesRootDir = '/afs/cern.ch/user/f/frs/public_html/mad-X_examples';

$htmlFile = '/afs/cern.ch/user/n/nougaret/www/mad/test.htm'; # for the time being

$pwd = `pwd`;
chop $pwd;
$localRootDir = $pwd;

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
    print "--- testing $target\n";

    $testReport .= "<table width=\"75%\" border=\"0\">\n";
    $testReport .= "<tr class='test_target'><td colspan=\"2\"><div align=\"center\"><strong>Testing $target</strong></div></td></tr>\n";

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

	$testReport .= "<tr class='test_case'><td width=\"70%\">$testCaseDir</td><td width=\"30%\"></td></tr>\n"; 
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
	@inputs = ($infilename, @inputs); # add the root inputfile
	

	# copyping inputs and dependent files locally
	foreach $input (@inputs) {
	    print "Input is $input\n";
	    $_ = $input;
	    # SPECIFIC CASE: files that must be stored in a locally replicated hierarchy
	    # considering the MAD call instructions mention a relative path...
	    if (/(..\/[\w\d\-_]+)\/([\w\d-_.]+)/) { # ../dir/file format only
		# file to be called is located up the hierarchy
		# in which case we need to reflect the directory tree structure
		# by creating directories if necessary
		$dependencyDir = $1;
		$dependencyFile = $2;
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

		# check whether this file in turns depends on another file
		@secondLevelInputs = split /,/, $dependencyList{"$target/$dependencyFile"}; 
		# prefixed with $target to avoid name clashes (may be too coarse!)
		foreach $secondLevelInput (@secondLevelInputs) {
		    print "Second-level copy of $secondLevelInput\n";
		    $_ = $secondLevelInput;
		    if (/(..\/[\w\d\-_]+)\/([\w\d-_.]+)/) { # ../dir/file format only
			$secondDependencyDir = $1;
			$secondDependencyFile = $2;
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
                # file to be called in the same directory as the input file
		# print "for target '$target' and test input '$infilename', now copying additional '$input'\n";
		if ($sourceSubDir eq "") {
		    `cp $samplesRootDir/$target/$input .`;
		} else {
		    `cp $samplesRootDir/$target/$sourceSubDir/$input .`;
		}
	    }

	}

	# before executing the command, make sure we remove
	# the optional subdirectory information which shows up after the comma
	# in case the source test directory contains a subdirectory structure...
	if ($sourceSubDir eq "") {
	    $executableCommand = $command;
	} else {
	    $_ = $command;
	    s/,[\s\t]*[\w\d.\-_=]+//g;
	    $executableCommand = $_;
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
	@allFilesNow = `ls -I mad*`; # list of all input + output files after invoking 'mad' command 
        # ignore the madx/madxp entries which should stay on top
	
	# remove the madx link from the list of files to be moved

	@outputs = ();
	foreach $file (@allFilesNow){
	    chop $file;
	    $isInput = 0;
	    foreach $input (@inputs) {
		if ($file eq $input) {
		    $isInput = 1;
		} # otherwise the file has been produced upon invoking MAD must go into the ouput subdir
	    }
	    if ($isInput == 0) {
		push (@outputs, $file);
	    }
	} ;

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
		$testReport .="<tr class='omit'><td>$file</td><td>NO FILE FOR REFERENCE</td></tr>\n";
	    } else {
		if ($sourceSubDir eq "") {
		    $diffRes = `diff ./$file $samplesRootDir/$target/$file | wc -l`;
		} else {
		    $diffRes =  `diff ./$file $samplesRootDir/$target/$sourceSubDir/$file | wc -l`;
		}
		#print OUT $diffRes;
		chop $diffRes;
		
		if ($diffRes == 0 ) {
		    print OUT "# COMPARING $file yields $diffRes different lines => SUCCESS\n";
		    $testReport .= "<tr class='success'><td>$file</td><td>SUCCESS</td></tr>\n";
		} else {

		    # attempt juggling with the files' contents to check whether the match is actually
		    # better than it looks
		    # if this is an .out output file, comparison should omit the header and footer
		    if ($file =~ /[\w_\-.\d]+.out/ ){
			$secondDiffRes = specificOutputComparison("./$file","$sourceSubDir");
			print OUT "# COMPARING $file yields $diffRes different lines";
			print OUT " => cleaning header/footer yields $secondDiffRes differences =>";
			if ($secondDiffRes ==0) {
			    print OUT "WARNING\n";
			    $testReport .="<tr class='warning'><td>$file</td><td>WARNING</td></tr>\n";
			} else {
			    print OUT "FAILURE\n";
			    $testReport .="<tr class='failure'><td>$file</td><td>FAILURE</td></tr>\n";
			}
		    } else {
			# this is a 'by-product' output file
			# otherwise by checking the numerical values' precision is preserved
			$secondDiffRes = roundedOutputComparison("./$file","$sourceSubDir");
			print OUT "# COMPARING $file yields $diffRes different lines";
			print OUT " => tampering with contents yields $secondDiffRes differences =>";
			if ($secondDiffRes ==0) {
			    print OUT "WARNING\n";
			    $testReport .="<tr class='warning'><td>$file</td><td>WARNING</td></tr>\n";
			} else {
			    print OUT "FAILURE\n";
			    $testReport .="<tr class='failure'><td>$file</td><td>FAILURE</td></tr>\n";
			}			
		    }
		}
	    }
	}
	close(OUT);	
      
    } # tests
    $targetReport .= "</table>\n"; # end of table per target

} # targets

$endTime = localtime;

# create web page
$html = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">';
$html .= "<html>\n";
$html .= "<head>\n";
$html .= "<title>MAD test result</title>\n";
$html .= "<link rel=stylesheet href=\"./MadTestWebStyle.css\" type=\"text/css\">"; # CSS stylesheet
$html .= "</head>\n";
$html .= "<!-- automatically generated by the MAD testing script -->\n";
$html .= "<body>\n";
$html .= "<p>Test started $startTime, ended $endTime</p>";
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
    # print "process $parentFilename";
    while(<IN>){
	# take into account MAD various ways of reading a file, namely using either "call,file" or "readmytable,file"
	my $fileRetreival = 0 ;
	# MAD syntax too permissive: hard to grep commands
	if (/[Rr][Ee][Aa][Dd][Mm][Yy][Tt][Aa][Bb][Ll][Ee],?[\s\t]*[Ff][Ii][Ll][Ee][\s\t]*=[\s\t]*[\"\']?([\w._\-\d\/]+)[\"\']?[\s\t]*;/){
	    my $child;
	    $child = $1;
	    $fileRetreival = 1;

	}
	if (/[Cc][Aa][Ll][Ll],?[\s\t]*[Ff][Ii][Ll][Ee][\s\t]*=[\s\t]*[\"\']?([\w._\-\d\/]+)[\"\']?[\s\t]*;/) {
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


sub specificOutputComparison {
    # specific matching of the .out files
    my $outfilename = $_[0];
    my $sourceSubDir = $_[1];
     
    # Specific pre-processing of the outfiles which feature headers and tailers that may change
    # from one run to the next.
    my $pwd = `pwd`;
    chop $pwd;
    # copy the original reference output file locally since we are going to clean-up from
    # its more or less standard headers and footers...
    my $referenceOutfilename = $outfilename . '_ref';
    
    if ($sourceSubDir eq "") {
	`cp $samplesRootDir/$target/$outfilename ./$referenceOutfilename`;
    } else {
	`cp $samplesRootDir/$target/$sourceSubDir/$outfilename ./$referenceOutfilename`;
    }


    print "now about to compare '$outfilename' and '$referenceOutfilename' in directory '$pwd'";
    
    # before performing the comparison, discard heading comments and tail of the output files
    # (such cooking required until we manage to have more structured i/o file formats)
    
    my $headerEndMark = "Execution Time Stamp";
    my $footerStartMark = "terminated - total number of elements";
    
    print "now about to clean files up\n";
    @filesToBeCleaned = ($outfilename,$referenceOutfilename);
    
    foreach $inputFileName (@filesToBeCleaned) {
	my $outputFileName = $inputFileName . "_clean";
	
	open(INFILE,$inputFileName);
	open(OUTFILE,">$outputFileName");
	
	my $copy = 0 ;
	my $neverWritten = 1; # to avoid case where no header end-mark is found, yielding to empty files

	while(<INFILE>){
	    $currentLine = $_;
	    if (/$footerStartMark/) { $copy = 0;}
	    if ($copy == 1) { 
		print OUTFILE $currentLine;
		$neverWritten = 0 ;
			  }
	    if (/$headerEndMark/){ $copy = 1; }
	}
	
	close(INFILE);
	close (OUTFILE);
    }
    print "files have been cleaned-up\n";
    
    my $differencesLines;

    if ($neverWritten == 1) {
	$differencesLines = "EMPTY_FILE";
    } else {

	my $cleanedOutfilename = $outfilename . "_clean";
	my $cleanedReferenceOutfilename = $referenceOutfilename . "_clean";
	
	$differencesLines = `diff $cleanedOutfilename $cleanedReferenceOutfilename | wc -l `;
	chop $differencesLines;
    }
    $_ = $differencesLines ;
};


sub roundedOutputComparison {
    my $leastSignificantDigits = 5 ; # for time-being consider 4 least significant digits may differ in exponential notation
    my $outfilename = $_[0];
    my $sourceSubDir = $_[1];   
    # very approximate test
    # replace all the least significant digits by 'round' mark
    # and then carry-out the diff once again
    # alternate test, check that the ratio is equal to one at the specified precision
     # Specific pre-processing of the outfiles which feature headers and tailers that may change
    # from one run to the next.

    my $pwd = `pwd`;
    chop $pwd;
    # copy the original reference output file locally since we are going to clean-up from
    # its more or less standard headers and footers...
    my $referenceOutfilename = join("",$outfilename,'_ref');

    if ($sourceSubDir eq "") {
	`cp $samplesRootDir/$target/$outfilename ./$referenceOutfilename`;
    } else {
	`cp $samplesRootDir/$target/$sourceSubDir/$outfilename ./$referenceOutfilename`;
    }
    print "now about to compare '$outfilename' and '$referenceOutfilename' in directory '$pwd'";
    
    print "now about to clean roundup files up\n";
    my @filesToBeRoundedUp = ($outfilename,$referenceOutfilename);
    
    foreach $inputFileName (@filesToBeRoundedUp) {
	my $outputFileName = $inputFileName . "_rounded";
	
	open(INFILE,$inputFileName);
	open(OUTFILE,">$outputFileName");
	
	my $copy = 0 ;
	
	while(<INFILE>){
	    goto notest;
	    if(/(([+-]?\d.\d+)(\d{$leastSignificantDigits})([eE][+-]?\d+))/) {
		print "pattern found in $1 with ";
		print "MSD $2, and ";
		print "LSD $3, and ";
		print "exponent $4\n"; exit;
	    }
	  notest:
	    # in the future the 'round' marker should be adapted in length so as to match
	    # the precision expressed in terms of least significant bits...
	    s/(([+-]?\d.\d+)(\d{$leastSignificantDigits})([eE][+-]?\d+))/\2round\4/g; # \1 stands for $1 in a regexp

	    # a set of comments which depend on the time at which the test is launched
	    s/@[\s\t]*ORIGIN[\s\t]+[\w\d\-\/_.%\s\t\"]*/@ REPLACED COMMENT THAT CHANGES FROM ONE MAD VERSION TO THE NEXT\n/g;
	    s/@[\s\t]*DATE[\s\t]+[\w\d\-\/_.%\s\t\"]*/@ REPLACED COMMENT THAT CHANGES FROM ONE RUN TO THE NEXT\n/g;
	    s/@[\s\t]*TIME[\s\t]+[\w\d\-\/_.%\s\t\"]*/@ REPLACED COMMENT THAT CHANGES FROM ONE RUN TO THE NEXT\n/g;

	    print OUTFILE $_;
	}
	
	close(INFILE);
	close (OUTFILE);
    }
    print "files have been rounded-up\n";
    
    my $roundedOutfilename = $outfilename . "_rounded";
    my $roundedReferenceOutfilename = $referenceOutfilename . "_rounded";
    
    my $differencesLines = `diff $roundedOutfilename $roundedReferenceOutfilename | wc -l `;
    chop $differencesLines;

    $_ = $differencesLines ; 

};
