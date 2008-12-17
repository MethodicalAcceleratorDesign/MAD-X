#!/usr/bin/perl

# trigger build and test if the following condition is satisfied:
# (1) find-out if the latest tag of the form 'madX_3_04_22' 
#     is located below the last 'madX_test' tag, which means
#     some release occured after the last test of the code.
# (2) tag the CVS with 'madX_test'.

# ... then issue a trigger to the do following:
# (3) do the build when extracting for the last release tag, ignoring
#     all changes that occured in between.
# (4) run the test suite

use MIME::Lite; # to send e-mail

use File::Path; # to remove directory trees


my $trigger;


@extractedPackages = ('madX');

$pwd = `pwd`;
chop $pwd;
$extractDir = $pwd . "/MadCvsExtract_trigger" ;
rmtree($extractDir);
mkdir($extractDir, 0777);
chdir($extractDir);


$cvsDir = ":gserver:isscvs.cern.ch:/local/reps/madx" ;


# Do we need to do the check-out or do we rely on FesaBuild.pl instead?
foreach(@extractedPackages) {
    my $pack = $_;
    # print "Extract package $pack from CVS\n";
    `cvs -d $cvsDir checkout $pack`;
}

# find-out the latest release
chdir('./madX');
my $representative = 'madxd.h';
my @log = `cvs log $representative`;

@releases = ();
@tests = ();
foreach $line (@log){
    if ($line =~/^[\s\t]*madX-(\d+)_(\d+)_(\d+)[\s\t]*:[\s\t]*([\d\.]+)[\s\t]*$/){
    	my $release = "$1_$2_$3";
	$release_revision{$release}=$4;
    	@releases = (@releases, $release);
	# print "found release $1_$2_$3, with revision $4\n";
    }
    if ($line =~/^[\s\t]*test-(\d+)_(\d+)_(\d+)[\s\t]*:[\s\t]*([\d\.]+)[\s\t]*$/){
    	my $test = "$1_$2_$3";
	$test_revision{$test}=$4;
    	@tests = (@tests, $test);
	# print "found test $1_$2_$3, with revision $4\n";
    }

}

my @sortedReleases = sort byDecreasingReleaseNumber @releases;
my @sortedTests = sort byDecreasingReleaseNumber @tests;

my $lastRelease = @sortedReleases[0];
my $beforeLastRelease = @sortedReleases[1];
my $lastTest = @sortedTests[0];

# decide whether a new release took place, in which case we shall
# trigger the build and test procedure.

if ($lastTest eq $lastRelease) {
	# there's no need to run the test again.
	$trigger = 'do-nothing';
} else {
	# also account for the very first time
	my $newTest = $lastRelease;
	my $newTestTag = "test-" . $newTest;
	# tag the CVS repository
	# ...
	`cvs tag $newTestTag $representative`;
	
	# now find-out all the work from the contributors
	# that went into the CVS in between the two last releases

	@authors = (); # global variable modified by recordWork
	recordWork("madX-$beforeLastRelease", "madX-$lastRelease");

	# sort @authors by alphabetical order
	@authors = sort @authors;

	my $workReport = "";
	$workReport .= "The MAD Build & Test script has detected a tentative candidate release for $lastRelease.\n\n";
	$workReport .= "Since last candidate release, the following changes have been made:\n";
	$workReport .= "\t-Lines-of-code added/deleted between $beforeLastRelease and $lastRelease:\n";
	foreach $auth (@authors){
		$workReport .= "\t\t$auth: +$linesAdded{$auth} -$linesDeleted{$auth}\n";
	}
	$workReport .= "\nSee detailed work report on:\n";
	$workReport .= "http://test-mad-automation.web.cern.ch/test-mad-automation/workReport.html\n";
	$workReport .= "\nFrom now on, the test procedure will start and may take several days.\n";
	$workReport .= "\nAt the end of the test procedure, module keepers will be informed ";
	$workReport .= "in case of discrepancy between the test's outcome and the reference. \n";

	# ... and send a summary to the list of watchers by e-mail

	$msg = MIME::Lite->new(
		       From       => 'Jean-Luc.Nougaret@cern.ch',
		       'Reply-To' => 'mad-automation-admin@cern.ch',
		       To         => 'mad-module-keepers@cern.ch',
		       Subject    => "MAD $lastRelease release candidate ready for testing.",
		       Data       => $workReport
		       );
	$msg->send;
	$trigger = 'run-test';
	# will cause the MadTest.pl script to run.

}
chdir($pwd); # back to the top menu
rmtree($extractDir);

# returned strings for MadBuildAndTest.pl
print "releaseTag=madX-$lastRelease\n"; # \n-terminated
print "trigger=$trigger\n"; # \n-terminated
exit 0;

sub byDecreasingReleaseNumber {
	# sort always assumes $a is to be compared with $b, such as '$f{$a} <=> f{$b}'
	# here, all numbers are of the form x+_x+_x+
	my @aNumbers = split /_/, $a;
	my @bNumbers = split /_/, $b;
	my $counter = 0;
	my $result = 0; # default
	foreach $aNumber (@aNumbers){
		if ($aNumber<@bNumbers[$counter]){
			$result = 1;
			return $result;
		}
		if ($aNumber>@bNumbers[$counter]){
			$result = -1;
			return $result;
		}
		$counter++;
	}
	return $result;
}

sub recordWork {
	$rel1 = $_[0]; # first release
	$rel2 = $_[1]; # second release
	@files = `ls *.*`;
	@makefiles = `ls Make*`;
	@files = (@files, @makefiles);
	# print "rel1 is $rel1, rel2 is $rel2\n";
	foreach $file (@files){
		chop $file;
		my @logs = `cvs log $file`;
		# first look for the revisions associated to the first and second releases
		foreach $log (@logs){

			my $pattern1 = "^[\\s\\t]*$rel1\[\\s\\t]*:[\\s\\t]*([\\d\\.]+)[\\s\\t]*\$";
			# print "pattern is $pattern\n";
			if ($log =~/$pattern1/){
				$rev1 = $1;
				# print "$file: $rel1 -> revision 1 is $rev1\n";
			}
			my $pattern2 = "^[\\s\\t]*$rel2\[\\s\\t]*:[\\s\\t]*([\\d\\.]+)[\\s\\t]*\$";
			if ($log =~ /$pattern2/){
				$rev2 = $1;
				# print "$file: $rel2 -> revision 2 is $rev2\n";
			}
		}
		
		# now look for the authors having modified the code between release 1 and 2
		# (including the latter).
		# to this end, check all revisions located between $rev1 and $rev2 for given file
		
		my @numbers1 = split /\./, $rev1;
		my @numbers2 = split /\./, $rev2;
		
		# re-scan log for this file in order to find-out revision work
		# in a CVS log info looks like the following:
		# revision: 1.2.3
		# date: ... author: ...
		my $retreiveDetailedWorkInfoMode = 0 ; # some state-variable

		LOG_SCAN: foreach $log (@logs){
			if ($retreiveDetailedWorkInfoMode == 1){
				#print "Mode: retreive detailed work info\n";
				#print "Process line $log\n";
				if ($log =~ /author[\s\t]*:[\s\t]*(\w+);/){
					my $candidate = $1;
					# print "found work-unit carried-out by $candidate\n";
					# make sure author is not already in the list
					my $alreadyRecorded = 0;
					foreach $auth (@authors){
						if ($auth eq $candidate) {
							$alreadyRecorded = 1;
						}
					}

					$log =~ /;[\s\t]+lines:[\s\t]+\+(\d+)[\s\t]+\-(\d+)[\s\t]*;/;
					# print "work is +$1 -$2 LOC\n";
					if ($alreadyRecorded == 0 ) {
						@authors = ( @authors, $candidate );
						$linesAdded{$candidate} = $1;
						$linesDeleted{$candidate} = $2;  
					} else {
						$linesAdded{$candidate} += $1;
						$linesDeleted{$candidate} += $2;
					
					}
				}
				$retreiveDetailedWorkInfoMode = 0;
			}
			if ($log =~ /^[\s\t]*revision[\s\t]+([\d\.]+)[\s\t]*$/){
				# same number of dots?
				my $rev = $1;
				my @numbers = split /\./, $rev; 
				# print "about to compare $rev1 < $rev < $rev2 ?\n";
				if (scalar(@numbers) == scalar(@numbers2)) {
					# now check that $rev1 < $rev <= $rev2
					# and record the author & work in this case
					# print "compare $rev1 < $rev < $rev2 ?\n";
					my $counter = 0;
					foreach $number (@numbers){
						if ($number<@numbers1[$counter]) { next LOG_SCAN; }
						if ($number>@numbers2[$counter]) { next LOG_SCAN; }
						$counter++;
					}
					
					# special case: $rev = $ rev1 => of no interest
					if ($rev eq $rev1) { next LOG_SCAN; }
						
					# if we reach this point, the comparison was succesful
					# =>advance to the next log line contains the author
					$retreiveDetailedWorkInfoMode = 1;
					# print "$file: contribution found $rev1 < $rev <= $rev2\n";
					
				} else {
					next LOG_SCAN;
				}
				
			}
		}
		
		
		
	}
}
