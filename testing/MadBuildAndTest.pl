#!/usr/bin/perl

# retreive directory hosting the Mad program from the command line
$thisProgramName = $0;
$_ = $thisProgramName;
/^([\w\d\-_.\/]+)\/MadBuildAndTest.pl/;

if ($1 ne "") {
	$hostDirectory = $1;
} else {
	$hostDirectory = ".";
}
chdir($hostDirectory);

my $runTest = `./MadTrigTest.pl`;

# at this stage ./MadCvsExtract/madx dir created locally
# (currently overwritten by MadBuild.pl)

if ($runTest eq "run-test") {

	# irrespective of which directory MadTrigTest.pl ended in,
	# MadBuild.pl starts from the location set by calling MadBuildAndTest.pl program

	`./MadBuild.pl`; # issues e-mail upon completion
	
	# irrespective of which directory MadBuild.pl ended in,
	# MadTest.pl starts from the location set by calling MadBuildAndTest.pl program
	
	`./MadTest.pl ./MadCvsExtract/madX`; # issues e-mail upon completion
}
