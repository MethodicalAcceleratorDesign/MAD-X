#!/usr/bin/perl

# assume a Perl script in charge of compiling MAD on the Windows platform
# is waiting on an agreed-upon socket => the Server must be started before
# the Client, if not the program should return an error otherwise it would
# wait for ever!

# accept an argument manual-trigger or automatic-trigger, the former meaning we
# immediately trigger the compilation, the later meaning the program will listen
# to port 7075 for a request emanating from the automatic build and test procedure.

#$debug = 'no'; # global variable also used by notify() subroutine and others
# this program only for tests => set $debug='yes' from 10 June 2009
$debug = 'yes';


if ($#ARGV!=0){
    die "expect one argument: 'now' or 'wait-for-trigger'\n";
} else {
    if ($ARGV[0] eq 'now'){
	$mode = 'now';
    } else {
	if ($ARGV[0] eq 'wait-for-trigger') {
	    $mode = 'wait-for-trigger';
	} else {
	    die "incorrect argument: should be either 'now' or 'wait-for-trigger'\n";
	}
    }
    
}


# KILL ITSELF IN CASE ALREADY RUNNING

my @check = `ps -aef | grep MadWindowsCompileClient.pl`;
foreach $line (@check){
    chop $line;
    $line =~ /^\w+[\s\t]+(\d+)/; # format of output: username pid parent ...
    my $pid = $1; # pid of the process being considered
    my $myPid = $$; # $$ is the pid of this process (special perl variable!#@?)
   
    if ($pid ne $myPid){
	# check this is indeed a perl command
	if ($line =~ /\/usr\/bin\/perl[\s\t]/) {
	    # `kill -9 $pid`; # kill process, unless this is this process's pid
	    # no: instead should kill itself
	    if ($mode eq 'now') {
		my $warning = "MadWindowsCompileClient.pl already running => abort new process.\n";
		print $warning;
		notify($warning);
	    } # otherwise don't print the message otherwise the cron job
	    # will send an e-mail
	    exit;
	} else {
	    # skip
	}
    } # later-on should send a signal that this process would either accept
    # to kill itself, or reject in case it is already engaged in a compilation
    # on the Windows platform => before dying, close the listening socket port
}



# fork process to spawn a branch that will periodically refresh the AFS kerberos tokens
my $child_pid = fork();

if (not defined $child_pid){
    notify("no system resources to fork process => exit");
    exit;
}

if ($child_pid==0){
    # this is the child process
    # refresh the AFS token every 6 hours. Otherwise the token
    # would expire after 25 hours.
    # (note that this trick works for up to 10 days according to IT support)
    my $start = localtime;
    # the recipe for refresing AFS/Kerberos tokens through a child process is known to work
    # up to ten days. As a consequence, this process should die gracefully before ten days
    # of running and wait for the acrontab to relaunch it with a time of grain of 5 min
    # Hence the NFS client that triggs the remote compilation should have a time-out of the
    # same order.
    my $counter = 0;
    while(1){
	$counter ++; # increment every 6 hours
	if ($counter > (4*5)) { # life-expectancy set to 5 days (could as well try with 9 days)
	    # time to die gracefully and let the acrontab restart a new process with
	    # fresh AFS and Kerberos tokens. In the meantime, the trigger should be able
	    # to wait for the client socket port to reappear...
	    # should also kill the parent process...
	    my $parent_pid = getppid();
	    kill 9, $parent_pid;
	    exit;
	}
	my $now = localtime;
	sleep 21600; # 6 hours
	`/usr/sue/bin/kinit -R`;
	`/usr/bin/aklog`;
	# check if the child process' parent is dead. If so, should kill itself
	my $parent_pid = getppid(); # get parent's pid
	$cnt = kill 0, $parent_pid;
	if ($cnt == 0){
	    exit;
	}
    }

}

# else ...

if ($child_pid){
    # non-zero pid means we are in the parent process, which received the child's pid



    # output of this program on stdout: SUCCESS or FAILURE:<message>

    # this program must run on a machine such as abcopl1 on which NFS mounted partitions
    # can be viewed through Samba by the remote Windows machine.

    my $windowsHost = 'abpc10788';
    # $windowsHost = 'abcopl1'; # 29 september 2009 - for test purposes

    $executablesAfsWebFolder = "/afs/cern.ch/user/n/nougaret/www/mad/windows-binaries"; # global

    $madForWindowsSambaFolder = "/user/nougaret/MAD-X-WINDOWS/madX"; 
    # problem: won't be seen on pcslux99!!! => cannot automate fully !!!
    # => for the time-being this process will need to be launched manually.

    # where binaries are delivered on the web for subsequent retreival by users

    $madWindowsCompilationDir = $madForWindowsSambaFolder; # global used by other routines
    $madWindowsDeliveryDir = "/afs/cern.ch/user/n/nougaret/www/mad/windows-binaries"; # global
    # also used by other routines
    @windowsTargets = ('madx.exe','madxp.exe','mpars.exe'); # Windows/DOS deliverables
    # above is global as used by other routines as well

    use IO::Socket::INET;
    use MIME::Lite;
    use Sys::Hostname;


    # two local and remote ports for communication with the Windows host in the two directions
    $socketPortWindows = 7070; # agreed-up with client (>1024 for non-root)
    $socketPortLinux = 7071; # could be the same as above
    # one local port to listen for requests emanating from the automated build and test process
    my $listeningPort = 7075;


    my $thisLinuxHost = hostname;


    # before asking the Windows host to trigger the compilation, we must first make sure that the
    # Samba folder MAD-X-WINDOWS/madX contains the latest CVS (more precisely latest released tagged
    # version - for the time being, we'll simply pick-up the latest contents of the repository)
    # ... well the `cvs update` is actually carried-out in sub updateMadForWindowsSambaFolder() !!!

    # wait to be waken-up (or start right now if mode is 'now')

    INFINITE_LOOP: while(1){

	my $invokeCompilation = 0; # default = do nothing

	if ($mode eq 'wait-for-trigger'){
	    my $listeningSocket = new IO::Socket::INET(
						       LocalHost => $thisLinuxHost,
						       LocalPort => $listeningPort,
						       Proto => 'tcp',
						       Listen => 1,
						       Reuse => 1,
						       ) or die "Can't bind : $@\n";

	    unless ($listeningSocket) {
		notify("failed to open TCP socket $listeningPort on $thisLinuxHost to receive command => will die\n");
		die 'failed to open TCP socket $listeningPort on $thisLinuxHost';
	    }

	    my $receive = $listeningSocket->accept();
	    while (<$receive>){
		if (/^Trigger Windows compilation$/){
		    # sent triggering signal via socket to Windows compilation server
		    $invokeCompilation = 1;
		    last;
		}
	    }
	    close $listeningSocket;
	} # else mode is 'now' and we should trigger the remote compilation right now
	
	
	if ($mode eq 'now'){
	    $invokeCompilation = 1 ;
	}


	if ($invokeCompilation == 1){
	    notify("MadWindowsCompileClient.pl will now forward the compilation request to the Windows host machine.");

	    updateMadForWindowsSambaFolder();



	    # $thisLinuxHost = 'abcopl1';
	    # print "the Linux box is '$thisLinuxHost'\n";

	    my $sock = new IO::Socket::INET ( 
					      PeerAddr => $windowsHost,
					      PeerPort => $socketPortWindows,
					      Proto => 'tcp'
					      ); 

	    unless ($sock) {
		notify("Could not create socket $socketPortWindows to connect to $windowsHost => will die\n");
		die "Could not create socket: $!\n" unless $sock; 
	    }

	    print "will now send message to port $socketPortWindows of $windowsHost\n";
	
	    print $sock "$thisLinuxHost asks: Compile MAD for Windows!\n";
	
	    $startTime = localtime; # global
	    $endTime; # global, will be set later-on
	
	    close($sock);
	
	    # now wait for the message signalling that the compilation completed



	    my $clientSock = new IO::Socket::INET(
						  LocalHost => $thisLinuxHost,
						  LocalPort => $socketPortLinux,
						  Proto => 'tcp',
						  Listen => 1,
						  Reuse => 1
						  );

	    die "Could not create client socket: $!\n" unless $clientSock;
	    
	    print "$thisLinuxHost accepts messages sent through socket $socketPortLinux\n";
	    my $newClientSock = $clientSock->accept();


	    
	    while (<$newClientSock>){
		print $_;
		if (/Compilation completed/){
		    $endTime = localtime;
		    print "OK: the compilation completed on Windows side\n";
		    checkWindowsCompilationOutcome();
		    print "=> installed the executables in the AFS web folder\n";
		    last INFINITE_LOOP; # leave the while loop
		}

		# should leave loop (timeout) in case there's no reply by the Windows-side server,
		# in which case, the executables will need to be delivered manually

	    }
	
	    close ($clientSock);
	    
	} # if $invokeCompilation == 1

	if ($mode eq 'now') {
	    # should compile only once and then leave the infinite loop to complete the program
	    last;
	} else {
	    print "now wait for wake-up by next compilation-triggering signal\n";
	}

	# debug:
	my $whereAmI = `pwd`;
	notify("at the end of the compilation, the Linux box client is in '$whereAmI'\n");

    } # while(1): wait forever to be woken-up by automated build-and-test program

    # do we really need to kill the child process?
    # in principle not, but the child would commit suicide only 6 hours later due to loop duration
    kill 9, $child_pid;

} # this is the parent process (not the child forked process refreshing AFS/Kerberos tokens)


sub checkWindowsCompilationOutcome {
    my $initialDir = `pwd`;
    # check the delivery directory contents
    foreach $target (@windowsTargets){
	# check that the executable has been created within the last hour
	my $ls = `ls -l $madWindowsCompilationDir/$target`;

	# debug
	notify("for target '$target', we see : '$ls'");

	# pick the date and time at which the executables have been created
	$ls =~ /(\w{3})[\s\t]+(\d{1,2})[\s\t]+(\d+:\d+)[\s\t]/ ;

	my $month = $1;
	my $day = $2;
	my $time = $3;

	my $now = localtime;
	print "now=$now\n";

	$now =~ /^\w{3}[\s\t]+(\w{3})[\s\t]+(\d{1,2})[\s\t]+(\d+:\d+:)\d+/ ;
	# forget about the year...

	my $monthNow = $1;
	my $dayNow = $2;
	my $time = $3;

	# debug
	notify("monthNow is '$monthNow', month is '$month', dayNow is '$dayNow', day is '$day'");
#	if (0){ # for the time being, always deliver the executables, without checking anything
	if (($monthNow != $month)||($dayNow != $day)){
	    print "Mistmatch of day and month => executables were not created\n";
	} else {
	    # now check that compilation occured within on hour from now
	    
	    # now install the executables in the AFS web folder
	    my $source = "$madWindowsCompilationDir/$target";
	    my $destination = "$madWindowsDeliveryDir/$target";
	    my $result = `cp $source $destination`;
	    
	    # debug:
	    notify("just copied '$source' into '$destination' => outcome = '$result'");

	}


    } # for each $target (@windowsTargets)

    # if everything ok...




    # now notify that the Windows executables are ready
    my $grepVersion = `grep myversion $madWindowsCompilationDir/madxd.h`; # hard-coded !?

    # debug:
    notify("now grep my version in '$madWindowsCompilationDir/madx.h'");

    $grepVersion =~ /MAD-X (\d+\.\d+\.\d+)/;
    $madVersion = $1; # global, also used in subroutine 'deliverHtmlPage';

    deliverHtmlPage();
    

    if ($debug eq 'no') {
	my $msg = MIME::Lite->new(
				  From => 'Jean-Luc.Nougaret@cern.ch',
				  To => 'mad-windows-watchers@cern.ch',
#				  To => 'Jean-Luc.Nougaret@cern.ch',
				  Subject => 'MAD-X for Windows updated',
				  Data => "Dear colleagues,\n\nPlease take note that MAD-X version $madVersion is now available on Windows.\n\nThe new releases are available for download on the new Web page:\nhttps://test-mad-automation.web.cern.ch/test-mad-automation/windows-binaries/executables.htm\n\nRegards,\nJean-Luc"
				  );
	$msg->send;
    } else {
	notify("MAD-X for Windows has been updated");
    }
    
    chdir($initialDir);

} # subroutine checkWindowsCompilationOutcome


sub deliverHtmlPage {
    
    # at this stage, the Windows binaries have been delivered to the
    # AFS web folder already

    my $htmlFile = "$executablesAfsWebFolder/executables.htm"; 
    # only for Windows?

    my $contents =''; # blank at first

    # grep size of the binaries located in the AFS web folder
    my @binaries = `ls -l $executablesAfsWebFolder/*.exe`;

#    my $nBinaries = scalar(@binaries);
#    notify("in '$executablesAfsWebFolder', 'found $nBinaries'");
	
    $contents .= "Version $madVersion compiled with Lahey Fortran and Microsoft Visual C++:\n";
    $contents .= "<table width=\"75%\" border=\"0\">\n";
    my $oddOrEven = 'even'; # to colorize successive lines differently
    foreach $binary (@binaries){
	chop $binary; # end of line
#	notify("line:$binary");
#	notify("in '$executablesAfsWebFolder', 'found $binary'");
# -rw-r--r--  1 nougaret pz  658664 Oct  1 12:06 /afs/cern.ch/user/n/nougaret/www/mad/windows-binaries/mpars.exe
	$binary =~ /(\d+)[\s\t]+(\w{3})[\s\t]+(\d{1,2})[\s\t]+(\d+:\d+)[\s\t]+[^\s]+\/(\w+\.exe)$/;
	my $size = $1;
	my $megabytes = $size / 1000000;
	my $month = $2;
	my $day = $3;
	my $time = $4;
	my $executable = $5;
#	notify("size='$size',exec='$executable',descr='$description{$executable}'");
	$description{'madx.exe'} = "standard version";
	$description{'madxp.exe'} = "version including PTC";
	$description{'mpars.exe'} = "\"parser-only\" version";
	if ($oddOrEven eq 'odd'){
	    $oddOrEven = 'even';
	} else {
	    $oddOrEven = 'odd';
	}
	$contents .= "<tr class=\"$oddOrEven\"><td>Download</td><td><a href=\"./$executable\">$executable</a></td><td>($megabytes Megabytes)</td><td>for the $description{$executable}.</td></tr>\n";
    }
    $contents .= "</table>\n";

    # create web page in the correct AFS web folder location
    my $html = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">';
    $html .= "<html>\n";
    $html .= "<head>\n";
    $html .= "<title>MAD-X downloadable executables</title>\n";
    $html .= "<link rel=stylesheet href=\"../MadTestWebStyle.css\" type=\"text/css\">"; # CSS stylesheet one level up
    $html .= "</head>\n";
    $html .= "<!-- generated by Windows compilation script -->\n";
    $html .= "<body>\n";
    $html .= "<p>Windows compilation started $startTime, ended $endTime</p>\n";
    $html .= $contents;
    $html .= "</body>\n";
    $html .= "</html>\n";
    open(OUTHTML, ">$htmlFile");
    print OUTHTML $html;
    close OUTHTML;

    # debug
    notify("created file '$htmlFile'");
    
    # now move HTML file into the AFS target web folder
    
}


sub updateMadForWindowsSambaFolder{
    my $localDir = `pwd`;
    chdir($madForWindowsSambaFolder);
    # ideally we should do a complete clean-up here.
    print "invoke CVS update in $madForWindowsSambaFolder. Ideally should do a complete clean-up before\n";
    `cvs update`;
    $cvsStatus = `cvs status`;
    if ($debug eq 'yes'){
	notify("outcome of `cvs update`: $cvsStatus");
    }
    chdir ($localDir); # back to where we were before entering the sub

}


sub notify{
    if ($debug eq 'yes') {
	my $message = $_[0];
	my $msg = MIME::Lite->new(
				  From => 'MAD-X Windows compilation robot',
				  ReplyTo => 'Jean-Luc.Nougaret@cern.ch',
				  To => 'Jean-Luc.Nougaret@cern.ch',
				  Subject => 'automatic notification',
				  Data => $message
				  );
	$msg->send;    
    } # else do nothing

}
