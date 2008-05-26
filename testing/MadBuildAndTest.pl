#!/usr/bin/perl

my $child_pid = fork();

if (not defined $child_pid){
	print "no system resources to fork process\n";
	exit;
}

if ($child_pid==0){
	# this child process
	# refresh the AFS token every 6 hours. Otherwise the token
	# would expire after 25 hours.
	# (note this trick works for up to 10 days according to IT support)
	INFINITE_LOOP:
	sleep 21600; # 6 hours
	`kinit -R`;
	`aklog`;
	
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

	my $releaseTag;
	my $runTest;
	my $debugMode;
	
	my $argsNumber = $#ARGV+1;
	
	if ($argsNumber>0){
		$debugMode = 1;
	} else {
		$debugMode = 0;
	}
	
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

	if ($debugMode==0){
		# should also retreive the latest release tag
		print "entering MadTrigTest.pl\n";
		my @return = `./MadTrigTest.pl`; # extracts MAD's CVS repository with no tag
		print "MadTrigTest.pl completed\n";

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
	} else {
		# expect "debug=madX-x_yy_zz" where madX-x_yy_zz is the release tag
		$firstArg = $ARGV[0];
		if ($firstArg =~ /debug=(madX\-\d_\d\d_\d\d)/){
			$releaseTag = $1;
			$runTest = 'run-test';
		} else {
			print "Expect either no argument (automatic trigger)or argument of form:\n";
			print "\tdebug=madX-x_yy_zz\n";
			print "Kill the child process\n";
			kill 9, $child_pid; # kill the child process
			print "Exit!\n";
			exit;
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
	
		print "entering MadBuild.pl\n";	
		`./MadBuild.pl $releaseTag`; # issues e-mail upon completion
		print "MadBuild.pl completed\n";
	
		# irrespective of which directory MadBuild.pl ended in,
		# MadTest.pl starts from the location set by calling MadBuildAndTest.pl program
	
		print "entering MadTest.pl\n";
		`./MadTest.pl ./MadCvsExtract/madX`; # issues e-mail upon completion
		print "MadTest.pl completed\n";
	} else {
		print "No new release detected => no need to run the test-suite\n";
	}

	# kill the child process in charge of refreshing the AFS token every 6 hours	
	kill 9, $child_pid; # kill the child process

} # parent-process