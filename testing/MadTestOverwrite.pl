#!/usr/bin/perl
# A script to invoke madx/madxp on a specific target, and to store the modified outputs into
# the madx-examples CVS repository, to alleviate the task of the people in charge of overwriting
# test cases whenever a new release brought unsignificant changes...


$target = $ARGV[0];
$root = "/afs/cern.ch/user/n/nougaret/scratch0/mad-automation";

@targets = `xsltproc --stringparam what list_targets $root/ProcessScenario.xsl $root/TestScenario.xml`;

$localDirectory = `pwd`; chop $localDirectory;

my $found = 0;
foreach $targ (@targets){
    chop $targ; # remove new line char
    if ($target eq $targ){
	$found = 1;
	print "about to retreive information about target '$target'\n";

	# do the CVS extraction locally
	my $localCvsExtractDir = "./localCvsExtractDir_$target";
	mkdir $localCvsExtractDir;
	chdir "$localCvsExtractDir";
	# do the CVS extraction
	my $cvsRootDir = ':kserver:isscvs.cern.ch:/local/reps/madx-examples';
	`cvs -d $cvsRootDir checkout madX-examples/REF/$target`;
	my $targetDir = "$localDirectory/$localCvsExtractDir/madX-examples/REF/$target";
	chdir $targetDir;
	$now=`pwd`; print "now in $now\n";

	my @tests = `xsltproc --stringparam target $target --stringparam what list_tests $root/ProcessScenario.xsl $root/TestScenario.xml`;
	foreach $test (@tests){
	    chop $test;
	    print "$test\n";
	    if ($test =~ /\.\/madx(p?)[\s\t]*<[\s\t]*([\w\d\.\_\-]+)[\s\t]*>[\s\t]*([\w\d\.\_\-]+)[\s\t]*,[\s\t]*subdirectory=([\w\d\-\_\.]+)/){
		# chop $test;
		if ($1 eq 'p'){
		    $program = 'madxp';
		} else {
		    $program = 'madx';
		}
		$input = $2;
		$output = $3;
		$subdir = $4;
		print "program is '$program' input is '$input', output is '$output', dir is '$subdir'\n";
		
		if ($subdir ne ""){
		    chdir $subdir; # going into the sub-directory of the test
		}
		
		# path for madx/madxp should link to the latest production version
		my $where = `pwd`; print "invoked from $where\n";
		my $programRoot = "/user/nougaret/MAD-X/madX"; 
		my $command = "$programRoot/$program < $input > $output";


		`$command`; # execute the command
		print "note: MAD commands will be picked up in $programRoot\n";

		# now look at the list of files that need CVS commit
		@locallyModified = `cvs status | grep Locally`;
		foreach $mod (@locallyModified){
		    chop $mod;
		    if ($mod =~ /^File:[\s\t]*([\w\-\_\.]+)[\s\t]*Status: Locally Modified[\s\t]*$/){
			my $file = $1;
			# need to commit this file
			print "now invoke 'cvs commit -m \"automated commit\"' $file\n";
			`cvs commit -m \"automated commit - overwrote all outputs by MadTestOverwrite.pl\" $file`;
		    }
		}


		if ($subdir ne ""){
		    chdir $targetDir; # going-back one level up
		}


	    } else {
		print "no match of the MAD command line!\n";
		exit 1;
	    } # failed to match MAD command line
	} # for each test
	chdir $targetDir; # back to where we were...
	# clean-up the location of the extract (rmtree...)
    } # target found (only once)

    
}

if ($found eq 0) {
    print "non-existing target '$target'\n";
    exit 1;
}


# retrieve the MAD command with input / output arguments

#`xsltproc --stringparam what retreiveXXX ProcessScenario.xsl TestScenario.xml`;

