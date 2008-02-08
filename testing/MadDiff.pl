#!/usr/bin/perl

# performs a diff of two files and displays the output in HTML

# inputs: the two files to be compared, and the name of the HTML file to be created (full path)
# output: prints SUCCESS, WARNING or ERROR of the comparison to stdout

if ( $#ARGV != 2 ) {
    print "expect 3 arguments: (1) first file (2) second file and (3) target htmlfile! EXIT!\n" ;
    exit ;
} else {
    $leftFilename = @ARGV[0];
    $rightFilename = @ARGV[1];
    $htmlFile = @ARGV[2];

    # debug
    print "MadDiff invoked for (1) '$leftFilename', (2) '$rightFilename' and (3)'$htmlFile'\n";
}

my $startTime = localtime;

my $diffReport = "";

$diffReport .= "<table>\n";
$diffReport .= "<tr><td width=\"50%\">$leftFilename</td><td width=\"50%\">$rightFilename</td></tr>\n";

my $diffResult = `diff --side-by-side --width 180 $leftFilename $rightFilename`; # should not exceed webpage max
my @lines = split /\n/, $diffResult;


foreach $line (@lines) {

    # first replace "<" and ">" with the "&lt;" "&gt;" equivalent in HTML otherwise heading for trouble in display
    $_ = $line;
    s/</&lt;/g;
    s/>/&gt;/g;
    $line = $_;

    # separate the left and right parts of the line

    my @parts = split /[|]/, $line; # look for parts that differ

    my $nbParts = scalar(@parts);


    if ($nbParts == 2) {
	# before concluding the difference is a failure, check it is not only a warning
	my @knownPatterns = (
			     # patterns showing-up in .out files
			     '^[\s\t]*\+[\s\t]+MAD-X[\s\t]+\d+\.\d+\.\d+[\s\t]+\+[\s\t]*$',
			     '^[\s\t]*\+[\s\t]+Code Modification Date:[\s\t]+\d+\.\d+\.\d+[\s\t]+\+[\s\t]*$',
			     '^[\s\t]*\+[\s\t]+Execution Time Stamp:[\s\t]+\d+\.\d+\.\d+[\s\t]+\d+\.\d+\.\d+[\s\t]+\+[\s\t]*$',
			     '^[\s\t]*\+[\s\t]+MAD-X[\s\t]+\d+\.\d+\.\d+[\s\t]+finished normally[\s\t]+\+[\s\t]*$',
			     # other patterns that often come together
			     '^[\s\t]*@[\s\t]+ORIGIN[\s\t]+%',
			     '^[\s\t]*@[\s\t]+DATE[\s\t]+%',
			     '^[\s\t]*@[\s\t]+TIME[\s\t]+%'
			  );

	my $matchKnownPattern = 0; # default
	foreach $pattern (@knownPatterns) {
	    if ( $parts[0] =~ /$pattern/) {
		if ($parts[1] =~ /$pattern/) {
		    $matchKnownPattern = 1;
		}
	    }
	}

	if ($matchKnownPattern == 1) {
	    $diffReport .= "<tr class=\"different-warning\"><td>$parts[0]</td><td>$parts[1]</td></tr>\n";
	} else {
	    # before concluding the difference is a failure, check for numerical rounding errors
	    $diffReport .= "<tr class=\"different-failure\"><td>$parts[0]</td><td>$parts[1]</td></tr>\n";
	}
    } else {
	

	if ($nbParts == 1){

	    # try to see whether left and right parts of the line are identical

	    # force split in the middle
	    $length = length($parts[0]);
	    # theoretically should be odd, accounting for the gutter...
	    my $leftPart;
	    my $rightPart;
	    if ($length%2==1) { 
		# odd
		# to be checked
		$leftPart = substr $line, 0, ($length-1)/2;
		$gutter = $line[($length-1)/2];
		$rightPart = substr $line, ($length-1)/2+1, ($length-1)/2; 
	    } else {
		# even
		# to be checked
		# NO GUTTER !!
		$leftPart = substr $line, 0, ($length/2);
		$rightPart = substr $line, ($length/2)+1, ($length/2);
	    }
	    
	    # rework alignments

	    $_ = $leftPart;
	    s/^[\s\t]+//g; # no tab nor space at the beginning
	    s/[\t\s]+/ /g; # replace multiple tabs and space as one space
	    s/[\t\s]+$//g; # no tab nor space at the end
	    $alignedLeftPart = $_;

	    $_ = $rightPart;
	    s/^[\s\t]+//g; # no tab nor space at the beginning
	    s/[\t\s]+/ /g;
	    s/[\t\s]+$//g; # no tab nor space at the end
	    $alignedRightPart = $_;


	    # make sure the left and right parts are really identical, otherwise MadDiff got lost
	    if ($leftPart eq $rightPart) {
		$diffReport .= "<tr class=\"identical\"><td>$leftPart</td><td>$rightPart</td></tr>\n";
	    } else {
		if ($alignedLeftPart eq $alignedRightPart) {
		    $diffReport .= "<tr class=\"almost-identical\"><td>$leftPart</td><td>$rightPart</td></tr>\n";
		}
		else {
		    # before concluding we're lost check whether the part only belongs to left or right
		    # diff --side-by-side issues a '<' for parts belonging to the first file only
		    # ! remember '<' has already been replaced by '&lt;' earlier on for HTML display...
		    if ($line =~ /&lt;$/) {
			# remove the '<' char at the end of the line
			$_ = $line;
			s/[\s\t]+&lt;//g;
			$cleanLine = $_;
			$diffReport .= "<tr class=\"only-left\"><td>$cleanLine</td><td></td></tr>\n";
		    }
		    # also try to the same for part only belonging to the right part...
		    else {
			# diff --side-by-side issues a '<' for parts belonging to the first file only
			# ! remember '>' has already been replaced by '&gt;' earlier on for HTML display...
			
			if ($line =~ /^[\s\t]*&gt;/) {
			   
			    # remove the '>' char at the beginning of the line
			    $_ = $line;
			    
			    # curious: following fails when trying to add spaces or tabs before ampersand
			    s/&gt;//g; # there seems to be a trouble around this line
				$cleanLine = $_;
			       
				$diffReport .= "<tr class=\"only-right\"><td></td><td>$cleanLine</td></tr>\n";
			}	
			else {
			    # finally ...
			    # keep the line as it is - just got lost!
			    $diffReport .= "<tr class=\"got-lost\"><td colspan=\"2\">$line</td></tr>\n";
			}
		    }
		}
	    }
	} # $nbParts = 1, i.e. could not split as two left and right parts separated by "|" (from diff)
    }


    # if 1 part, there must be a whitespace in the middle, meaning identical
    
    # if more than two parts there is a problem


}

my $endTime = localtime;

$diffReport .= "</table>\n";

# create web page
my $html = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">';
$html .= "<html>\n";
$html .= "<head>\n";
$html .= "<title>Diff of $leftFilename and $rightFileName</title>\n";
$html .= "<link rel=stylesheet href=\"../MadTestWebStyle.css\" type=\"text/css\">"; # CSS stylesheet one level up
$html .= "</head>\n";
$html .= "<!-- automatically generated by the MAD test script -->\n";
$html .= "<body>\n";
$html .= "<p>Test started $startTime, ended $endTime</p>\n";
$html .= $diffReport;
$html .= "</body>\n";
$html .= "</html>\n";
open(OUTHTML, ">$htmlFile");
print OUTHTML $html;
close OUTHTML;
