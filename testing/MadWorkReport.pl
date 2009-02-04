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

$website = "/afs/cern.ch/user/n/nougaret/www/mad/";

@extractedPackages = ('madX');

$pwd = `pwd`;
chop $pwd;
$extractDir = $pwd . "/MadCvsExtract_testWorkReport" ;
rmtree($extractDir);
mkdir($extractDir, 0777);
chdir($extractDir);


$cvsDir = ":gserver:isscvs.cern.ch:/local/reps/madx" ;

my $startTime = localtime;

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

	
# now find-out all the work from the contributors
# that went into the CVS in between the two last releases

@authors = (); # global variable modified by recordWork
recordWork("madX-$beforeLastRelease", "madX-$lastRelease");


chdir($pwd); # back to the top menu
rmtree($extractDir);

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

sub recordWork { # should also pass the output file name
    $rel1 = $_[0]; # first release
    $rel2 = $_[1]; # second release

    $rootDir = "$website"; 
    chop $rootDir; # don't need the final "/"
    $workReportFile = "workReport.html";
    open WORKFILE, ">$rootDir/$workReportFile";
    $workReportContent = "";

    @files = `ls *.*`;
    @makefiles = `ls Make*`;
    @files = (@files, @makefiles);
    # print "rel1 is $rel1, rel2 is $rel2\n";
    FILE_LOOP: foreach $file (@files){
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
	
	if ($rev1 eq $rev2){
	    # no work went into latest release, may skip
	    next FILE_LOOP;
	}

	# now look for the authors having modified the code between release 1 and 2
	# (including the latter).
	# to this end, check all revisions located between $rev1 and $rev2 for given file
	
	# re-scan log for this file in order to find-out revision work
	# in a CVS log info looks like the following:
	# revision: 1.2.3
	# date: ... author: ...
	my $parserState = "skipHeader" ; # some state-variable
	my $nextParserState = ""; # default		
	my @logsOfInterest = `cvs log -r$rev1:$rev2 $file`;

#	print "cvs log -r$rev1:$rev2 $file\n";

	foreach $log (@logsOfInterest){
  
#	    print "Current Parser State = $parserState\n";
    
	    if ($parserState eq "skipHeader"){
		if ($log =~ /^description:[\s\t]*/){
		    $nextParserState = "skipLine";
		} else {
		    $nextParserState = "skipHeader";
		}
	    }

	    if ($parserState eq "skipLine"){
		if ($log =~ /^\-\-\-(\-)+[\s\t]*/){
		    $nextParserState = "readRev";
		} else {
		    if ($log =~ /^\=\=\=(\=)+[\s\t]*/){
			$nextParserState = "EOF";
		    } else {
			print "line '-----' expected, but not found!\n";
			return; # leave function
			# should clean-up the HTML file containing the work-report first
		    }
		}
	    }
	    
	    
	    if ($parserState eq "readRevisionInfo"){
		#print "Mode: retreive detailed work info\n";
		#print "Process line $log\n";
		if ($log =~ /author[\s\t]*:[\s\t]*(\w+);/){
		    $author = $1;
		    # print "found work-unit carried-out by $candidate\n";
		    # make sure author is not already in the list
		    my $alreadyRecorded = 0;
		    foreach $auth (@authors){
			if ($auth eq $author) {
			    $alreadyRecorded = 1;
			}
		    }
		    
		    $log =~ /;[\s\t]+lines:[\s\t]+\+(\d+)[\s\t]+\-(\d+)[\s\t]*;/;
		    # print "work is +$1 -$2 LOC\n";
		    if ($alreadyRecorded == 0 ) {
			@authors = ( @authors, $author );
			$workByAuthor{$author}="";
			$linesAdded{$author} = $1;
			$linesDeleted{$author} = $2;  
		    } else {
			$linesAdded{$author} += $1;
			$linesDeleted{$author} += $2;
			
		    }
		}
		$nextParserState = "readComment" ;
	    }



	    if ($parserState eq "skipRevisionInfo"){
		# in case the work has been registered for the previous release, we want to discard it here
		# (a parallel sequential machine dedicated to rejecting work from previous release)
		$nextParserState = "skipComment";
	    }
	    if ($parserState eq "skipComment"){ # also, just in case the work relates to the previous release
		# (a parallel sequential machine dedicated to rejecting work from previous release)
		if ($log =~ /\-\-\-\-(\-)+/){
		    # actually we find-out we should actually be in "skipLine" state
		    $nextParserState = "readRev";
		} else {
		    if ($log =~ /(========)+/){
			$nextParserState = "EOF";
		    } else {
			$nextParserState = "skipComment"; # until encounter a '---' or '===' marker
		    }
		}
	    }


	    
	    if ($parserState eq "readComment"){
		

		if ($log =~ /\-\-\-\-(\-)+/){
		    # actually we find-out we should actually be in "skipLine" state
		    # mark the end of the comment with a specific beacon "EndOfComment"
		    $workByAuthor{$author} .= "EndOfComment\n";
		    $nextParserState = "readRev";
		} else {
		    if ($log =~ /(========)+/){
			# same as above
			$workByAuthor{$author} .= "EndOfComment\n";
			$nextParserState = "EOF";
		    } else {
			# should concatenate the comment's string here...
			$workByAuthor{$author} .= "$log";
			$nextParserState = "readComment"; # until encounter a '---' or '===' marker
		    }
		}

#		$nextParserState = "skipLine"; # unless the comment spreads over successive lines!
		

	    }

	    if ($parserState eq "readRev"){
		if ($log =~ /^[\s\t]*revision[\s\t]+([\d\.]+)[\s\t]*$/){
		    # same number of dots?
		    my $rev = $1;
		    
		    if ($rev eq $rev1){
			# print "revision $rev of $file belongs to last release => skip!\n";
			# don't want to register work that went in to the previous release: skip
			$nextParserState = "skipRevisionInfo";
		    } else {
			
			# if we reach this point, the comparison was succesful
			# =>advance to the next log line contains the author
			
			$nextParserState = "readRevisionInfo";
			# print "$file: contribution found $rev1 < $rev <= $rev2\n";

		    }
		    
			
		} else {
		    print "Encountered error while parsing file\n";
		    exit;
		} # if ($log
	    } # if (parserState


	    $parserState = $nextParserState;
	    
	} # foreach $log  line
	

    } # for each file

    # sort list of authors in alphabetical order
    @authors = sort (@authors);

    foreach $author (@authors){
	$workReportContent .= "<tr class=\"work-author\"><td>$author</td></tr>\n";
	@workComments = split /EndOfComment\n/, $workByAuthor{$author}; # "EndOfComment" is our specific beacon

	@cleanedWorkComments =();
	# remove multiple entries, e.g. occuring when several files are commited at the same time
	SCAN_COMMENTS: foreach $comment (@workComments){
	    foreach $registeredComment (@cleanedWorkComments){
		if ($comment eq $registeredComment){
		    next SCAN_COMMENTS;
		}
	    }
	    @cleanedWorkComments = ( @cleanedWorkComments, $comment);
	}

	
	# new HTML table line for each comment
	$oddOrEven = "even";
	foreach $comment (@cleanedWorkComments){

	    $workReportContent .= "<tr class=\"work-log-$oddOrEven\"><td>$comment</td></tr>\n";

	    if ($oddOrEven eq "even"){
		$oddOrEven = "odd";
	    } else {
		$oddOrEven = "even";
	    }

	}
    }

    my $endTime = localtime; # note that startTime defined outside this subroutine

    # create web page
    my $html = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">';
    $html .= "<html>\n";
    $html .= "<head>\n";
    $html .= "<title>Work log between releases $rel1 and $rel2 </title>\n";
    $html .= "<link rel=stylesheet href=\"./MadTestWebStyle.css\" type=\"text/css\">"; # CSS stylesheet at same level
    $html .= "</head>\n";
    $html .= "<!-- automatically generated by the MAD test script -->\n";
    $html .= "<body>\n";
    $html .= "<p>Work log between releases $rel1 and $rel2 </p>\n";   
    $html .= "<p>Log report started $startTime, ended $endTime</p>\n";
    $html .= "<table width=\"75%\" border=\"0\">\n";
    $html .= $workReportContent;
    $html .= "</table>\n";
    $html .= "</body>\n";
    $html .= "</html>\n";

    
    # also copy the work report into the MAD-X News html page
    my $newsFile = "/afs/cern.ch/eng/sl/MAD-X/pro/docu/Introduction/news.html";
    open NEWSFILE, "<$newsFile";
    my $news = '';
    while(<NEWSFILE>){
	my $line = $_;
	if ($line =~ /[\s\t]*<!\-\- do not remove \- START WORK LOG \- do not remove \-\->[\s\t]*/) { 
	    # insert contents of the work report here
	    $news .= "<!-- do not remove - START WORK LOG - do not remove -->\n";
	    $news .= "<p>Work log between releases $rel1 and $rel2 </p>\n";   
	    $news .= "<p>Log report started $startTime, ended $endTime</p>\n";
	    $news .= "<table width=\"75%\" border=\"0\">\n";
	    $news .= $workReportContent;
	    $news .= "</table>\n";
	} else{
	    $news .= $line;
	}
    }
    close NEWSFILE;
    open NEWSFILE, ">$newsFile";
    print NEWSFILE $news;
    close NEWSFILE;

    print WORKFILE $html;
    close WORKFILE;

} # sub
