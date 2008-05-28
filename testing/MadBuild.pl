#!/usr/bin/perl

# acrontab command: 25 * * * * lxplus (/afs/cern.ch/user/n/nougaret/myMAD/madX/AUTOMATION/MadBuild.pl) > /afs/cern.ch/user/n/nougaret/MADbuilt.out 2>&1


use MIME::Lite; # to send e-mail

if ( $#ARGV != 0 ) {
    print "expect 1 argument (1) releaseTag! EXIT!\n" ;
    exit ;
} else {
    $releaseTag = @ARGV[0];
}

open REPORT_FILE, ">MadBuild_Report.txt";

$startTime = localtime;

#$htmlFile = '/afs/cern.ch/user/n/nougaret/www/mad/build.htm'; # for the time being

$website = "/afs/cern.ch/user/n/nougaret/www/mad/";

# first set environment variables for lf95 and NAG compilers
# this is necessary for the acron job which does not perform a login
# that would set the variables transparently.

my $path = $ENV{'PATH'}; # should handle case of empty string
my $newPath = $path . ":/afs/cern.ch/sw/fortran/nag/f95.5.361/bin:/afs/cern.ch/sw/fortran/lahey/lf9562/bin";
$ENV{'PATH'}=$newPath;

my $ldLibPath = $ENV{'LD_LIBRARY_PATH'}; # should handle case of empty string
my $newLdLibPath = $ldLlibPath . ":/afs/cern.ch/sw/fortran/lahey/lf9562/lib";
$ENV{'LD_LIBRARY_PATH'}=$newLdLibPath;

# flexlm license manager for NAG compiler
$ENV{'NAGF95_ROOT'}="/afs/cern.ch/sw/fortran/nag/f95.5.361";
$ENV{'LM_LICENSE_FILE'}="/afs/cern.ch/sw/fortran/nag/f95.5.361/license.dat";

@extractedPackages = ('madX');


$pwd = `pwd`;
chop $pwd;
$extractDir = join("", $pwd, "/MadCvsExtract") ;


mkdir($extractDir, 0777);
chdir($extractDir);


$cvsDir = ":kserver:isscvs.cern.ch:/local/reps/madx" ;


foreach(@extractedPackages) {
    my $pack = $_;
    print REPORT_FILE "Extract package $pack from CVS\n";
    `cvs -d $cvsDir checkout -r $releaseTag $pack`; # now extract for specified release tag
}

# build
my $buildReport = "";

chdir('./madX');

$buildReport = "<table width=\"75%\" border=\"0\">\n";

@makefiles = ("Makefile_develop","Makefile_nag","Makefile"); # last one overwrites madx, madxp

@targets = ("madx","madxp");

foreach $makefile (@makefiles){

    foreach $target (@targets){

	my $detailedBuildReport = "";
	$detailedBuildReport = "<table width=\"75%\" border=\"0\">\n";
	my $startTime  = localtime;

	`make clean`;
	`rm $target`;
	`rm $target\_$makefile`; # later-on, $target copied into $target\_$makefile ...

	my $warnings = 0; # default

	my $makeResult = `make -f $makefile $target 2>&1`; 
	
	# if we wanted to redirect stdout and stderr into separate files...
        # system("make -f $makefile $target 1>build_result_stdout 2>build_result_stderr");

	my $warnings = lookForWarnings($makeResult);
	# colorize warnings
	$makeResult =~ s/([wW])arning(s?)/<font class=\"warning-font\">\1arning\2<\/font>/g;
	# undo the above in case the warning appears as '0 warnings'
	$makeResult =~ s/([\s\t]+0[\s\t]+)<font class=\"warning\-font\">([Ww])arning(s?)<\/font>/\1\2arning\3/g;
	# make sure line feeds display correctly in HTML
	$makeResult =~ s/\n/<\/td><\/tr>\n<tr><td>/g;
	$makeResult = "<tr><td colspan=\"2\">" . $makeResult . "</td></tr>\n\n";

	
	my $nbOfTargets = `ls $target | wc -w`;
	if ($nbOfTargets == 1) {
		# keep all the targets as madx_Makefile_develop,
		# madx_Makefile_nag and madx_Makefile, as needed by MadTest.pl
		`cp $target $target\_$makefile`;
		if ($warnings==0) {
			$compilationOutcome{$target} = 'success';
		} else {
			$compilationOutcome{$target} = 'warning';
		}
	}
	else { 
	    $compilationOutcome{$target} = 'failure';
	}
	$detailedBuildReport .= "<tr class =\"$compilationOutcome{$target}\"><td colspan=\"2\">$target</td><td>$compilationOutcome{$target}</td></tr>\n";
	$detailedBuildReport .= "<tr><td>$makeResult</td><tr>\n";

	my $outcome = $compilationOutcome{$target};
	my $endTime = localtime;
	$detailedBuildReport .= "</table>\n";
	my $htmlFile = "build_" . $makefile . "_" . $target . ".htm";
	$buildReport .= "<tr class=\"$outcome\"><td>make -f $makefile $target</td><td><a href=\"$htmlFile\">$outcome</a></td></tr>\n";
	
	createWebPage($htmlFile,$detailedBuildReport,$startTime,$endTime);

    } # for each target
	
} # for each makefile

`make clean`; # final make clean, mostly to remove all the .o files

$buildReport .= "</table>\n";

$endTime = localtime;

createWebPage("build.htm",$buildReport, $startTime, $endTime ); # main page



# then send an e-mail
  $msg = MIME::Lite->new(
			 From       => 'Jean-Luc.Nougaret@cern.ch', # if "" instead of '', '@'->'\@'
			 'Reply-To' => 'mad-automation-admin@cern.ch',
			 To         => 'mad-automation-admin@cern.ch',
			 Subject    => "Automated MAD Build $compilationOutcome{'madx'} for madx, $compilationOutcome{'madxp'} for madxp",
			 Data       => "This is an automated e-mail. Check report on\nhttp://test-mad-automation.web.cern.ch/test-mad-automation/build.htm"
			);
  $msg->send;

 close REPORT_FILE;

sub lookForWarnings {

	my $makeOutput = $_[0]; # global variable
	my $warnings = 0;
	my @lines = split /\n/, $makeOutput;
	
	foreach $line (@lines) {
		while ( $line =~ /([\s\t]*([\w\d_\-:]+)[\s\t]+[Ww]arning(s?))/g ){
			# /g modifier to make it a progressive match
			my $prefix = $2;
			if ($prefix ne "0") { $warnings=1;}
		}
	}
	
	return $warnings;
}
  
sub createWebPage {
    my $file = $website . $_[0];
    my $buildReport = $_[1];
    my $startTime = $_[2];
    my $endTime = $_[3];
    $html = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">';
    $html .= "<html>\n";
    $html .= "<head>\n";
    $html .= "<title>MAD build result</title>\n";
    $html .= "<link rel=stylesheet href=\"./MadTestWebStyle.css\" type=\"text/css\">"; # CSS stylesheet
    $html .= "</head>\n";
    $html .= "<!-- automatically generated by the MAD build script -->\n";
    $html .= "<body>\n";
    $html .= "<p>Build started $startTime, ended $endTime</p>\n";
    $html .= $buildReport;
    $html .= "</body>\n";
    $html .= "</html>\n";   
    open(OUT_HTML,">$file");
    print OUT_HTML $html;
    close OUT_HTML;
}
