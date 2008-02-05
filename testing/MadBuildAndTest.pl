#!/usr/bin/perl

# retreive directory hosting the Mad program from the command line
$thisProgramName = $0;
$_ = $thisProgramName;
/^([\w\d\-_.\/]+)\/MadBuildAndTest.pl/;

$hostDirectory = $1;
chdir($hostDirectory);

`./MadBuild.pl`;
print "build completed\n";
# at this stage ./MadCvsExtract/madx dir created locally
`./MadTest.pl ./MadCvsExtract/madX`;
print "test completed\n";
