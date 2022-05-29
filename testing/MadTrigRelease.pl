#!/usr/bin/perl

# trigger build and test if the following condition is satisfied:
# (1) find-out if the latest tag of the form 'madX_3_04_22_prod' 
#     is located below the last 'madX_prod' tag, which means
#     some release occured after the last release of the code.
# (2) tag the CVS with 'madX_prod'.
# (3) send work report to hep-project-madx@cern.ch
# (4) trigger Windows compilation and send mails to madx-windows-watchers@cern.ch

# this process is supposed to complete in less than a day, hence it is unnecessary
# to bother about refreshing the Kerberos/AFS tokens as we do in MadBuildAndTest.pl for instance...


# if already running, this script should die or kill the existing instance...

use MIME::Lite; # to send e-mail

use File::Path; # to remove directory trees

$debug = 'no'; # this one is default and is required in automatic mode
#$debug = 'yes'; # this one for manual tests to avoid sending e-mails to the community, and avoid tagging the CVS

@extractedPackages = ('madX');

$rootDir = '/afs/cern.ch/user/n/nougaret/scratch0/mad-automation';

$extractDir = $rootDir . "/MadCvsExtract_prod_assert" ;
#rmtree($extractDir); # should check it exists first
mkdir($extractDir, 0777);
chdir($extractDir);

$cvsDir = ":gserver:isscvs.cern.ch:/local/reps/madx" ;

# Redo a specific checkout
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
@prods = ();
foreach $line (@log){
    if ($line =~/^[\s\t]*madX-(\d+)_(\d+)_(\d+)_prod[\s\t]*:[\s\t]*([\d\.]+)[\s\t]*$/){
    	my $release = "$1_$2_$3";
	$release_revision{$release}=$4;
    	@releases = (@releases, $release);
	 notify("found release $1_$2_$3_prod, with revision $4\n");
    }
    if ($line =~/^[\s\t]*prod-(\d+)_(\d+)_(\d+)[\s\t]*:[\s\t]*([\d\.]+)[\s\t]*$/){
    	my $prod = "$1_$2_$3";
	$prod_revision{$prod}=$4;
    	@prods = (@prods, $prod);
	# print "found prod $1_$2_$3, with revision $4\n";
    }

}

my @sortedReleases = sort byDecreasingReleaseNumber @releases;
my @sortedProds = sort byDecreasingReleaseNumber @prods;

my $lastRelease = @sortedReleases[0];
my $beforeLastRelease = @sortedReleases[1];
my $lastProd = @sortedProds[0];

# decide whether a new release took place

if (($lastProd eq $lastRelease) && ($debug ne 'yes')) {
    # there's no need to release to production.
    chdir($rootDir);
    rmtree($extractDir);
    notify("MadTrigRelease.pl completed (NO RELEASE).");
    exit 0; 
} else {
    # also account for the very first time
    my $newProd = $lastRelease;
    my $newProdTag = "prod-" . $newProd;
    # tag the CVS repository
    # ...
    if ($debug ne "yes"){
	`cvs tag $newProdTag $representative`;
    }
    
    # generate HTML page for the work report between the two production releases
    my $beforeLastTag = "madX-$beforeLastRelease"."_prod";
    my $lastTag = "madX-$lastRelease"."_prod";
    `$rootDir/MadWorkBetweenReleases.pl $beforeLastTag $lastTag`;
    

    # now find-out all the work from the contributors
    # that went into the CVS in between the two last releases
    
    @authors = (); # global variable modified by recordWork
    recordWork($beforeLastTag, $lastTag); # between production releases
    
    # sort @authors by alphabetical order
    @authors = sort @authors;
    
    my $workReport = "";
    $workReport .= "MAD-X $lastRelease has been released.\n\n";
    $workReport .= "Since last release, the following changes have been made:\n";
    $workReport .= "\t-Lines-of-code added/deleted between $beforeLastRelease and $lastRelease:\n";
    foreach $auth (@authors){
	$workReport .= "\t\t$auth: +$linesAdded{$auth} -$linesDeleted{$auth}\n";
    }
    $workReport .= "\nSee detailed work log on:\n";
    $workReport .= "http://test-mad-automation.web.cern.ch/test-mad-automation/workReport.html\n";
    
    # ... and send a summary to the list of watchers by e-mail
    if ($debug ne 'yes'){
	$msg = MIME::Lite->new(
			       From       => 'Jean-Luc.Nougaret@cern.ch',
			       'Reply-To' => 'mad-automation-admin@cern.ch',
			       To         => 'hep-project-madx@cern.ch',
			       Subject    => "MAD $lastRelease released.",
			       Data       => $workReport
			       );
	$msg->send;
    }
    #
    # now trigger Windows compilation (note: this one runs as acrontab process with access to both AFS and NFS)
    #
    use IO::Socket::INET;
    use MIME::Lite;
    use Sys::Hostname;

    my $windowsHost = 'abpc10788';
    my $thisLinuxHost = hostname;
    my $socketPortWindows = 7070;
    my $socketPortLinux = 7071; # could be the same as above (different machine)

    # following should be globals as they are used from subroutines...
    $executablesAfsWebFolder = "/afs/cern.ch/user/n/nougaret/www/mad/windows-binaries"; # global
    $madForWindowsSambaFolder = "/user/nougaret/MAD-X-WINDOWS/madX"; 
    # where binaries are delivered on the web for subsequent retreival by users
    $madWindowsCompilationDir = $madForWindowsSambaFolder; # global used by other routines
    $madWindowsDeliveryDir = "/afs/cern.ch/user/n/nougaret/www/mad/windows-binaries"; # global
    # also used by other routines
    # @windowsTargets = ('madx.exe','madxp.exe','mpars.exe'); # Windows/DOS deliverables
    @windowsTargets = ('madx.exe'); # Windows/DOS deliverables
    # above is global as used by other routines as well

    notify("MadWindowsCompileClient.pl will now forward the compilation request to the Windows host machine.");

    updateMadForWindowsSambaFolder();

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
	    
    INFINITE_LOOP: while (<$newClientSock>){
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
    chdir($rootDir); # back to the top menu
    rmtree($extractDir);
    notify("MadTrigRelease.pl completed (RELEASED).");
    if ($debug ne 'yes'){
	$msg = MIME::Lite->new(
			       From       => 'Jean-Luc.Nougaret@cern.ch',
			       'Reply-To' => 'mad-automation-admin@cern.ch',
			       To         => 'mad-windows-watchers@cern.ch',
			       Subject    => "MAD-X for Windows updated",
			       Data       => "Dear colleagues,\n\nPlease take note that MAD-X version $madVersion is now available on Windows.\n\nThe new releases are available for download on the new Web page:\nhttps://test-mad-automation.web.cern.ch/test-mad-automation/windows-binaries/executables.htm\n\nRegards\nJean-Luc"
			       );
	$msg->send;
    }
    
    exit 0; 
}


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
#				  To => 'mad-windows-watchers@cern.ch',
				  To => 'Jean-Luc.Nougaret@cern.ch',
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
    # my @binaries = `ls -l $executablesAfsWebFolder/*.exe`;
    my @binaries = `ls -l $executablesAfsWebFolder/madx.exe`; # from March 26th 2009, only one executable

#    my $nBinaries = scalar(@binaries);
#    notify("in '$executablesAfsWebFolder', 'found $nBinaries'");
	
    $contents .= "<p>Version $madVersion compiled with Lahey Fortran and Microsoft Visual C++:</p>\n";
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
#	$description{'madxp.exe'} = "version including PTC";
#	$description{'mpars.exe'} = "\"parser-only\" version";
	if ($oddOrEven eq 'odd'){
	    $oddOrEven = 'even';
	} else {
	    $oddOrEven = 'odd';
	}
	$contents .= "<tr class=\"$oddOrEven\"><td>Download</td><td><a href=\"./$executable\">$executable</a></td><td>($megabytes Megabytes)</td><td>for the $description{$executable}.</td></tr>\n";
    }
    $contents .= "</table>\n";
    $contents .= "<p>Version 3.04.53 accepting sequences files with BV flag, as until March 2009:</p>\n";
    $contents .= "<table width=\"75%\" border=\"0\">\n";
    $oldExecutable = "";
    $megabytes = 123456789; # should put the actual value here
    $contents .= "<tr class=\"even\"><td>Download</td><td><a href=\"./madx-old.exe\">madx-old.exe</a></td><td>(2.6132 Megabytes)</td><td>for the archived version, without PTC.</td></tr>\n";
    $contents .= "<tr class=\"odd\"><td>Download</td><td><a href=\"./madxp-old.exe\">madxp-old.exe</a></td><td>(6.7554 Megabytes)</td><td>for the archived version, including PTC.</td></tr>\n";
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
    if ($debug eq 'yes'){
	my $message = $_[0];
	my $msg = MIME::Lite->new(
				  From => 'MAD-X Windows compilation robot',
				  ReplyTo => 'Jean-Luc.Nougaret@cern.ch',
				  To => 'Jean-Luc.Nougaret@cern.ch',
				  Subject => 'automatic notification',
				  Data => $message
				  );
	$msg->send;    
    }
}

