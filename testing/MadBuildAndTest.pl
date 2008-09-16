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

open REPORT_FILE, ">MadBuildAndTest_Report.txt";
my $now = localtime;
print REPORT_FILE "MadBuildAndTest.pl report from $now\n";
	
my $releaseTag;
my $runTest;
my $debugMode;
	
my $argsNumber = $#ARGV+1;
	
if ($argsNumber>0){
	$debugMode = 1;
} else {
	$debugMode = 0;
}
	


if ($debugMode==0){
	# should also retreive the latest release tag
	print REPORT_FILE "entering MadTrigTest.pl\n";
	my @return = `./MadTrigTest.pl`; # extracts MAD's CVS repository with no tag
	print REPORT_FILE "MadTrigTest.pl completed\n";

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
			print REPORT_FILE "Unknown string $retString returned by ./MadTrigTest.pl\n";
		}
	}
} else {
	# expect "debug=madX-x_yy_zz" where madX-x_yy_zz is the release tag
	$firstArg = $ARGV[0];
	if ($firstArg =~ /debug=(madX\-\d_\d\d_\d\d)/){
		$releaseTag = $1;
		$runTest = 'run-test';
	} else {
		print "Expect either no argument (automatic trigger)or argument of form:\n";
		print "\tdebug=madX-x_yy_zz\n";
		# no need to kill the child process: it was not yet forked
		# print "Kill the child process\n";
		#kill 9, $child_pid; # kill the child process
		print REPORT_FILE "Exit!\n";
		exit;
	}
}

my $child_pid;

if ($runTest eq "run-test") {
	$child_pid = fork();
} else {
	print REPORT_FILE "No new release detected => no need to run the test-suite\n";
	close REPORT_FILE;
	exit;
}


if (not defined $child_pid){
	print "no system resources to fork process\n";
	exit;
}

if ($child_pid==0){
    # this child process
    # refresh the AFS token every 6 hours. Otherwise the token
    # would expire after 25 hours.
    # (note this trick works for up to 10 days according to IT support)
    my $start = localtime;
    open TICKETS_HISTORY, ">MadBuildAndTest_Tickets_History.txt"; # first time opening
    print TICKETS_HISTORY "Tracking AFS/Kerberos tickets refreshing since $start\n";
    close TICKETS_HISTORY; # will be successively opened and closed to force flushing

  INFINITE_LOOP:
    open TICKETS_HISTORY, ">>MadBuildAndTest_Tickets_History.txt"; # append   
    my @tokens = `/usr/bin/tokens`;	# provide path to cope with reduced acron environment
    my @klist = `/usr/sue/bin/klist`;	# provide path to cope with reduced acron environment
    my $now = localtime;
    print TICKETS_HISTORY "\n\nAFS and Kerberos tickets on $now\n";
    print TICKETS_HISTORY "======================= running tokens\n";
    print TICKETS_HISTORY @tokens;
    print TICKETS_HISTORY "======================= running klist\n";
    print TICKETS_HISTORY @klist;
    print TICKETS_HISTORY "now trying to invoke 'kinit -R' and 'aklog'\n";
    close TICKETS_HISTORY; # open / close in place of flushing
    sleep 21600; # 6 hours
    `/usr/sue/bin/kinit -R`;	# provide path to cope with reduced acron environment
    `/usr/bin/aklog`;		# provide path to cope with reduced acron environment
    # could also set enviromment variables within perl,
    # as done within MadBuild.pl for flexlm for instance...
    
    # check if the child process' parent is dead, it should also kill itself
    $parent_pid = getppid(); # get parent process' pid
    $cnt = kill 0, $parent_pid;
    if ($cnt == 0){
	exit;
    }   
    goto INFINITE_LOOP;
}

if ($child_pid) {
	# non-zero pid means we are in the parent process, which received the child's pid

	# print REPORT_FILE "\$runTest is '$runTest', and \$releaseTag is '$releaseTag'\n";

	# at this stage ./MadCvsExtract/madx dir created locally
	# (currently overwritten by MadBuild.pl)
    
    
        # first create the work-log report
        `./MadWorkReport.pl`; # creates and deletes a CVS extract locally

	# the build-and-test procedure always applies to the last release
	# even if CVS commits took place afterwards. 

	# irrespective of which directory MadTrigTest.pl ended in,
	# MadBuild.pl starts from the location set by calling MadBuildAndTest.pl program
	
	print REPORT_FILE "entering MadBuild.pl\n";	
	`./MadBuild.pl $releaseTag`; # issues e-mail upon completion
	print REPORT_FILE "MadBuild.pl completed\n";
	
	# irrespective of which directory MadBuild.pl ended in,
	# MadTest.pl starts from the location set by calling MadBuildAndTest.pl program
	
	print REPORT_FILE "entering MadTest.pl\n";
	`./MadTest.pl ./MadCvsExtract/madX`; # issues e-mail upon completion
	print REPORT_FILE "MadTest.pl completed\n";

	close REPORT_FILE;
	
	# kill the child process in charge of refreshing the AFS token every 6 hours	
	kill 9, $child_pid; # kill the child process

} # parent-process
