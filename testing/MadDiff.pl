#!/usr/bin/perl

# performs a diff of two files and displays the output in HTML

# inputs: the two files to be compared, and the name of the HTML file to be created (full path)
# output: prints 'success','quasi-success', 'warning' or 'failure' of the comparison to stdout

my $tolerance = 0.001; # for time-being, hard-code maximum incertitude as 0.1%, i.e. 0.001
my $absoluteTolerance = 10e-15;
my $maxWidth = 360; # different lines may appear as identical after truncation
# for instance numerical discrepancies may not show up until width is increased
# to say 360 or so. In many cases 180 or so would be sufficient.

if ( $#ARGV != 2 ) {
    print "expect 3 arguments: (1) first file (2) second file and (3) target htmlfile! EXIT!\n" ;
    exit ;
} else {
    $leftFilename = @ARGV[0];
    $rightFilename = @ARGV[1];
    $htmlFile = @ARGV[2];

}


my $retStatus = 'undefined'; # default

my $startTime = localtime;

my $diffReport = "";

# we know $rightFilename can be rather long - of the form .../.../.../REF/...
# let's only display from ref
$shortenRightFilename = $rightFilename;
$shortenRightFilename =~ s/(.+)\/REF\/(.+)/\(cvs\)\/REF\/\2/;

$diffReport .= "<table>\n";
$diffReport .= "<tr><td width=\"50%\">$leftFilename</td><td width=\"50%\">$shortenRightFilename</td></tr>\n";

my $diffResult = `diff --side-by-side --width $maxWidth $leftFilename $rightFilename`; # should not exceed webpage max
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

	my $matchKnownPattern = 0; # default

	# --- specific case: 'Number of warning[s]:(\d)' match in case actual numbers on both sides are equal
	$warningPattern = '^[\s\t]*Number of warning[s]?[\s\t]*:[\s\t]*(\d+)[\s\t]*$';
	if ($parts[0] =~ /$warningPattern/){
	    my $leftWarningNumber = $1;
	    if ($parts[1] =~ /$warningPattern/){
		my $rightWarningNumber = $1;
		if ($leftWarningNumber == $rightWarningNumber){
		    $matchKnownPattern = 1 ; # will display and return status as a quasi-success, 
		    # i.e. as no warning nor failure
		} else { goto matchKnownPattern; }
	    } else { goto matchKnownPattern; }
	} else { goto matchKnownPattern; }
	# --- 

      matchKnownPattern:
 

	# before concluding the difference is a failure, check it is not only a warning
	# to be checked: looks like patterns do not works a strings in some cases
	my @knownPatterns = (
			     # patterns showing-up in .out files
			     '^[\s\t]*\+[\s\t]+MAD-X[\s\t]+\d+\.\d+\.\d+[\s\t]+\+[\s\t]*$',
			     '^[\s\t]*\+[\s\t]+Code Modification Date:[\s\t]+\d+\.\d+\.\d+[\s\t]+\+[\s\t]*$',
			     '^[\s\t]*\+[\s\t]+Execution Time Stamp:[\s\t]+\d+\.\d+\.\d+[\s\t]+\d+\.\d+\.\d+[\s\t]+\+[\s\t]*$',
			     '^[\s\t]*\+[\s\t]+MAD-X[\s\t]+\d+\.\d+\.\d+[\s\t]+finished normally[\s\t]+\+[\s\t]*$',
			     '^[\s\t]*\+[\s\t]*\w+ Version[\-\w:!\s]*\+[\s\t]*$',
			     # patterns showing-up in postcript output
			     '%%CreationDate:[\s\t]*\d\d\/\d\d\/\d\d',
			     '\(MAD-X[\s\t]+\d+\.\d+.\d+[\s\t]+\d\d\/\d\d\/\d\d[\s\t]+\d\d\.\d\d\.\d\d\)',
			     # other patterns that often come together
			     '^[\s\t]*@[\s\t]+ORIGIN[\s\t]+%',
			     '^[\s\t]*@[\s\t]+DATE[\s\t]+%',
			     '^[\s\t]*@[\s\t]+TIME[\s\t]+%',
			     # other patterns
			     '^[\s\t]*sec.s since start: \d+[\s\t]+since last call: \d+',
			     # patterns for gnuplot
			     '^[\s\t]*set title "no-title[\s\t]+MAD-X[\s\t]+\d+\.\d+\.\d+[\s\t]+\d\d\/\d\d\/\d\d[\s\t]+\d\d\.\d\d\.\d\d"',
			     # patterns specific to ptc_normal
			     '^[\s\t]*Reading Data took[\s\t]+[\d\.]+[\s\t]+second\(s\) of Execution Time',
			     '^[\s\t]*Program[\s\t]+&lt;[\w\d]+&gt;[\s\t]+used[\s\t]+[\d\.]+[\s\t]+second\(s\) of Execution Time', 
                             # above, <> have been replaced by &lt; &gt; respectively
			     '^[\s\t]*Program[\s\t]+&lt;[\w\d]+&gt;[\s\t]+used a total of[\s\t]+[\d\.]+[\s\t]+second\(s\) of Execution Time' 
			     # above, <> have been replaced by &lt; &gt; respectively
			  );


	foreach $pattern (@knownPatterns) {
	    if ( $parts[0] =~ /$pattern/) {
		if ($parts[1] =~ /$pattern/) {
		    $matchKnownPattern = 1;
		}
	    }
	}

	if ($matchKnownPattern == 1) {
	    $diffReport .= "<tr class=\"different-warning\"><td>$parts[0]</td><td>$parts[1]</td></tr>\n";
	    if (($retStatus eq "undefined")||($retStatus eq "success")) { 
		$retStatus = "quasi-success";
	    } # otherwise keep its worst value
	} else {

	    # should repeat the test for identical left and right here

	    my $leftPart = $parts[0];
	    $leftPart =~ s/^[\s\t]+//g; # no tab nor space at the beginning
	    $leftPart =~ s/[\t\s]+/ /g; # replace multiple tabs and space as one space
	    $leftPart =~ s/[\t\s]+$//g; # no tab nor space at the end
	    $leftPart =~ s/[\n\r\f]//g; # no newline, no carriage return, no form-feed (maybe useless)
	    

	    my $rightPart = $parts[1];
	    $rightPart =~ s/^[\s\t]+//g; # no tab nor space at the beginning
	    $rightPart =~ s/[\t\s]+/ /g;
	    $rightPart =~ s/[\t\s]+$//g; # no tab nor space at the end
	    $rightPart =~ s/[\n\r\f]//g; # no newline, no carriage return, no form-feed (maybe useless)

	    
	    if ( $leftPart eq $rightPart ) {
	   	$diffReport .= "<tr class=\"quasi-success\"><td>$parts[0]</td><td>$parts[1]</td></tr>\n";
	    	if (($retStatus eq "undefined")||($retStatus eq "success")) { 
		$retStatus = "quasi-success";
	   	 } # otherwise keep its worst value
	    } else {
	    
		    # before concluding the difference is a failure, check for numerical rounding errors
		    $numericalMatch = 1; # default
#		    @leftChunks = split /[\s\t:=]+/, $parts[0];
#		    @rightChunks = split /[\s\t:=]+/, $parts[1];

		    @leftChunks = split /[\s\t:=]+/, $leftPart;
		    @rightChunks = split /[\s\t:=]+/, $rightPart;
	    
		    if (scalar(@leftChunks) == scalar(@rightChunks)) {

	    
			CHUNK: for ($i=0; $i<scalar(@leftChunks); $i++) {
			    $leftChunk = $leftChunks[$i];
			    $rightChunk = $rightChunks[$i];

			    if ($leftChunk eq $rightChunk) { next CHUNK; } # substrings match, a fortiori numbers
			    if (($leftChunk =~ /^[\+\-]?0$/) && ($rightChunk =~ /^[\+\-]?0$/ )){ next CHUNK; }
			    # one is '0' while the other '-0', which may happen with Fortran signed zeroes
			    else {
				# numerical test
				# tolerance expressed in %. 0 % means zero tolerance down to double float precision

				# for some reason, looks like we have to double the '\' when the pattern is in a string
				# and also have to add "\" before the right-anchor "$"
				# ... to be clarified...
				# try to replace " by '
				$numPattern = '^[\+\-]?\d+\.?\d*[eE]?[\+\-]?\d*$'; # account for various formats found with Mad
				if (($leftChunk =~ /$numPattern/) && ($rightChunk =~ /$numPattern/)) {
				    $leftValue = $leftChunk;
				    $rightValue = $rightChunk;

				    if ($tolerance==0) {
					if ($leftValue == $rightValue ) { $numericalMatch = 1; 	} 
					else {
					    $numericalMatch = 0 ;
					    # force to leave the for-loop
					    $i = scalar(@leftChunks);
					}
				    }
				    else {
					# define relative tolerance as the maximum incertitude between the two values
					$mean = ($leftValue+$rightValue)/2.0;
					if ($mean != 0.0){
					    if ( (abs($leftValue-$rightValue)/abs($mean)) < $tolerance ) {
						$numericalMatch = 1;
					    } else {
						# still need to compare with absolute tolerance
						if ((abs($leftValue)<$absoluteTolerance) && (abs($rightValue)<$absoluteTolerance)) {
						    $numericalMatch = 1;
						}
						else {
						    $numericalMatch = 0 ;
						    # force to leave the for-loop
						    $i = scalar(@leftChunks);
						} 
					    }
					} else { # mean == 0.0 - special case
					    if (abs($leftValue) > ($tolerance/2.0)){
						# still need to compare with absolute tolerance
						if (abs($leftValue)<$absoluteTolerance){
						    $numericalMatch = 1;
						} else {
						    $numericalMatch = 0;
						    # force to leave the for-loop
						    $i = scalar(@leftChunks);
						}
					    } else {
						$numericalMatch = 1;
					    }
					}
				    
				    }

			} else { 
			    $numericalMatch = 0;
			    # force to leave the for-loop
			    $i = scalar(@leftChunks);
			}
			
		    }
		}

	    } else { $numericalMatch = 0; } # don't bother looking for the details

	    if ($numericalMatch == 1) {
		$diffReport .= "<tr class=\"numerical-match\"><td>$parts[0]</td><td>$parts[1]</td></tr>\n";
		if ($retStatus ne "failure") { $retStatus = "warning"; } # otherwise keep worst value
	    } else {
		$diffReport .= "<tr class=\"different-failure\"><td>$parts[0]</td><td>$parts[1]</td></tr>\n";
		$retStatus = "failure"; # anyway
	    }
	    
	    }
	    
	    
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
		$gutter = substr $line, ($length-1)/2, 1;
		$rightPart = substr $line, ($length-1)/2+1, ($length-1)/2; 
	    } else {
		# even
		# to be checked
		# NO GUTTER !!
		$leftPart = substr $line, 0, ($length/2);
		$rightPart = substr $line, ($length/2), ($length/2);
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
		if ($retStatus eq "undefined") {
		    $retStatus = "success";
		} # otherwise keep its previous value
	    } else {
		if ($alignedLeftPart eq $alignedRightPart) {
		    $diffReport .= "<tr class=\"almost-identical\"><td>$leftPart</td><td>$rightPart</td></tr>\n";
		    if ($retStatus eq "undefined"){
			$retStatus = "success";
		    }# otherwise keep its worst value
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
			$retStatus = "failure"; # whatever its value so far
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
			    $retStatus = "failure"; # whatever its value so far...
			}	
			else {
				# one last chance: 'left part|' or '|right part'
				if ($line =~ /^\|/) {
					$txt = $line;
					$txt =~ s/^\|//g; # to be checked
					$diffReport .= "<tr class=\"only-right\"><td></td><td>$txt</td></tr>\n";
			    		$retStatus = "failure"; # whatever previous val	
				} else {
					if ($line =~ /\|$/){
						$txt = $line;
						$txt =~ s/\|$//g; # to be checked
						$diffReport .= "<tr class=\"only-left\"><td>$txt</td><td></td></tr>\n";
			    			$retStatus = "failure"; # whatever previous val					
					} else {
		
						# finally ...
			    			# keep the line as it is - just got lost!
			    			$diffReport .= "<tr class=\"got-lost\"><td colspan=\"2\">$line</td></tr>\n";
			    			$retStatus = "failure"; # what ever its value so far...			
					}
				}
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
$html .="<p>Legend:</p>\n";
$html .="<table>\n";
$html .="<tr class=\"identical\"><td>Lines match</td></tr>\n";
my $tolPercent = $tolerance * 100.0;
my $tolAbsolute = $absoluteTolerance;
$html .="<tr class=\"numerical-match\"><td>Numerical data match within $tolPercent% relative tolerance, or $tolAbsolute absolute tolerance</td></tr>\n";
$html .="<tr class=\"different-warning\"><td>Lines differ expectedly</td></tr>\n";
$html .="<tr class=\"different-failure\"><td>Lines differ unexpectedly</td></tr>\n";
$html .="<tr class=\"only-left\"><td>Part present in one file, absent of the other</td></tr>\n";
$html .="<tr class=\"got-lost\"><td>File differencing program got lost</td></tr>\n";
$html .="</table>\n";
$html .= $diffReport;
$html .= "</body>\n";
$html .= "</html>\n";
open(OUTHTML, ">$htmlFile");
print OUTHTML $html;
close OUTHTML;


# return to STDOUT

print $retStatus;
