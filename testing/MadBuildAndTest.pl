#!/usr/bin/perl
chdir("/user/nougaret/MAD/bin"); # temporary
`./MadBuild.pl`;
print "build completed\n";
# at this stage ./MadCvsExtract/madx dir created locally
`./MadTest.pl ./MadCvsExtract/madX`;
print "test completed\n";
