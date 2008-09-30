#!/usr/bin/perl

# assume a Perl script in charge of compiling MAD on the Windows platform
# is waiting on an agreed-upon socket => the Server must be started before
# the Client, if not the program should return an error otherwise it would
# wait for ever!

# output of this program on stdout: SUCCESS or FAILURE:<message>

my $windowsHost = 'abpc10788';
# $windowsHost = 'abcopl1'; # 29 september 2009 - for test purposes

$executablesAfsWebFolder = "/afs/cern.ch/user/n/nougaret/www/mad/windows-binaries"; # global

$madForWindowsSambaFolder = "/user/nougaret/MAD-X-WINDOWS/madX"; 
# problem: won't be seen on pcslux99!!! => cannot automate fully !!!
# => for the time-being this process will need to be launched manually.

# where binaries are delivered on the web for subsequent retreival by users

my $madWindowsCompilationDir = $madForWindowsSambaFolder;
my $madWindowsDeliveryDir = "/afs/cern.ch/user/n/nougaret/www/mad/windows-binaries";
my @windowsTargets = ('madx.exe','madxp.exe','mpars.exe'); # Windows/DOS deliverables

use IO::Socket::INET;
use MIME::Lite;
use Sys::Hostname;


$socketPortWindows = 7070; # agreed-up with client (>1024 for non-root)
$socketPortLinux = 7071; # could be the same as above

my $thisLinuxHost = hostname;


# before asking the Windows host to trigger the compilation, we must first make sure that the
# Samba folder MAD-X-WINDOWS/madX contains the latest CVS (more precisely the latest released tagged
# version - for the time being, we'll simply pick-up the latest contents of the repository)

updateMadForWindowsSambaFolder();



# $thisLinuxHost = 'abcopl1';
# print "the Linux box is '$thisLinuxHost'\n";

my $sock = new IO::Socket::INET ( 
				  PeerAddr => $windowsHost,
				  PeerPort => $socketPortWindows,
				  Proto => 'tcp'
				  ); 

die "Could not create socket: $!\n" unless $sock; 

print "will now send message to port $socketPortWindows of $windowsHost\n";

print $sock "$thisLinuxHost asks: Compile MAD for Windows!\n";

my $startTime = localtime;
my $endTime; # will be set later-on

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
	last; # leave the while loop
    }

    # should leave this loop (timeout) in case there's no reply by the Windows-side server,
    # in which case, the executables will need to be delivered manually

}

close ($clientSock);


sub checkWindowsCompilationOutcome {
    my $initialDir = `pwd`;
    # check the delivery directory contents
    foreach $target (@windowsTargets){
	# check that the executable has been created within the last hour
	my $ls = `ls -l $madWindowsCompilationDir/$target`;
	# pick the date and time at which the executables have been created
	$ls =~ /(\w{3})[\s\t](\d{1,2})[\s\t](\d+:\d+)[\s\t]/ ;

	my $month = $1;
	my $day = $2;
	my $time = $3;

	my $now = localtime;
	print "now=$now\n";

	$now =~ /^\w{3}[\s\t](\w{3})[\s\t](\d{1,2})[\s\t](\d+:\d+:)\d+/ ;
	# forget about the year...

	my $monthNow = $1;
	my $dayNow = $2;
	my $time = $3;

#	if (0){ # for the time being, always deliver the executables, without checking anything
	if (($monthNow != $month)||($dayNow != $day)){
	    print "Mistmatch of day and month => executables were not created\n";
	} else {
	    # now check that compilation occured within on hour from now
	    
	    # now install the executables in the AFS web folder
	    my $source = "$madWindowsCompilationDir/$target";
	    my $destination = "$madWindowsDeliveryDir/$target";
	    `cp $source $destination`;
	}


    } # for each $target (@windowsTargets)

    # if everything ok...




    # now notify that the Windows executables are ready
    my $grepVersion = `grep myversion $madWindowsCompilationDir/madxd.h`; # hard-coded !?
    $grepVersion =~ /MAD-X (\d+\.\d+\.\d+)/;
    $madVersion = $1; # global, also used in subroutine 'deliverHtmlPage';

    deliverHtmlPage();
    

    my $msg = MIME::Lite->new(
			      From => 'Jean-Luc.Nougaret@cern.ch',
			      ReplyTo => 'Jean-Luc.Nougaret@cern.ch',
			      To => 'Jean-Luc.Nougaret@cern.ch',
			      Subject => 'MAD-X for Windows updated',
			      Data => "Dear colleagues,\n\nPlease take note that MAD-X version $madVersion is now available on Windows.\n\nThe new releases are available for download on the new Web page:\nhttps://test-mad-automation.web.cern.ch/test-mad-automation/windows-binaries/executables.htm\n\nRegards,\nJean-Luc"
			      );
    $msg->send;
    
    chdir($initialDir);

} # subroutine


sub deliverHtmlPage {
    
    # at this stage, the Windows binaries have been delivered to the
    # AFS web folder already

    my $htmlFile = "$executablesAfsWebFolder/executables.htm"; 
    # only for Windows?

    my $contents =''; # blank at first

    # grep size of the binaries located in the AFS web folder
    my @binaries = `ls -l $executablesAfsWebFolder/*.exe`;
    $contents .= "Version $madVersion compiled with Lahey Fortran and Microsoft Visual C++:\n";
    $contents .= "<table>\n";
    foreach $binary (@binaries){
	chop $binary; # end of line
	$binary =~ /(\d+)\s(\w{3})\s(\d{1,2})\s(\d+:\d+)\s[^\s]+\/(\w+\.exe)$/;
	my $size = $1;
	my $megabytes = $size / 1000000;
	my $month = $2;
	my $day = $3;
	my $time = $4;
	my $executable = $5;
	$description{'madx.exe'} = "standard version";
	$description{'madxp.exe'} = "version including PTC";
	$description{'mpars.exe'} = "\"parser-only\" version";
	$contents .= "<tr><td>Download</td><td><a href=\"./$executable\">$executable</a></td><td>($megabytes Megabytes)</td><td>for the $description{$executable}.</td></tr>\n";
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


    # now move HTML file into the AFS target web folder
    
}


sub updateMadForWindowsSambaFolder{
    my $localDir = `pwd`;
    chdir($madForWindowsSambaFolder);
    # ideally we should do a complete clean-up here.
    print "invoke CVS update in $madForWindowsSambaFolder. Ideally should do a complete clean-up before\n";
    `cvs update`;
    chdir ($localDir); # back to where we were before entering the sub

}
