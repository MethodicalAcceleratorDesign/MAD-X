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

# should also retreive the latest release tag
my @return = `./MadTrigTest.pl`; # extracts MAD's CVS repository with no tag

my $runTest;
my $releaseTag;
foreach $retString (@return){
	my $recognized = 0;
	if ($retString =~ /^trigger[\s\t]*=[\s\t]*([\w\-]+)$/){
		$runTest = $1;
		$recognized = 1;
	}
	if ($retString =~ /^releaseTag[\s\t]*=[\s\t]*([\w\-\d_]+)$/){
		$releaseTag = $1;
		$recognized = 1;
	}
	if ($recognized == 0 ) {
		print "Unknown string $retString returned by ./MadTrigTest.pl\n";
	}
}

# print "\$runTest is '$runTest', and \$releaseTag is '$releaseTag'\n";

# at this stage ./MadCvsExtract/madx dir created locally
# (currently overwritten by MadBuild.pl)

if ($runTest eq "run-test") {

	# the build-and-test procedure always applies to the last release
	# even if CVS commits took place afterwards. 

	# irrespective of which directory MadTrigTest.pl ended in,
	# MadBuild.pl starts from the location set by calling MadBuildAndTest.pl program
	
	`./MadBuild.pl $releaseTag`; # issues e-mail upon completion
	
	# irrespective of which directory MadBuild.pl ended in,
	# MadTest.pl starts from the location set by calling MadBuildAndTest.pl program
	
	`./MadTest.pl ./MadCvsExtract/madX`; # issues e-mail upon completion
}
