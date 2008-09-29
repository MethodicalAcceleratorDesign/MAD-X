#!/usr/bin/perl

# assume a Perl script in charge of compiling MAD on the Windows platform
# is waiting on an agreed-upon socket => the Server must be started before
# the Client, if not the program should return an error otherwise it would
# wait for ever!

# output of this program on stdout: SUCCESS or FAILURE:<message>

my $windowsHost = 'abpc10788';
# $windowsHost = 'abcopl1'; # 29 september 2009 - for test purposes

my $madWindowsCompilationDir = "/user/nougaret/MAD-X-WINDOWS/madX";
my $madWindowsDeliveryDir = "/afs/cern.ch/user/n/nougaret/www/mad/windows-binaries";
my @windowsTargets = ('madx.exe','madxp.exe','mpars.exe'); # Windows/DOS deliverables

use IO::Socket::INET;
use MIME::Lite;
use Sys::Hostname;


$socketPortWindows = 7070; # agreed-up with client (<1024 for non-root)
$socketPortLinux = 7071; # could be the same as above

my $thisLinuxHost = hostname;



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

print "$thisLinuxHost accepts messages sent throught socket $socketPortLinux\n";
my $newClientSock = $clientSock->accept();

while (<$newClientSock>){
    print $_;
    if (/Compilation completed/){
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

	if (0){ # for the time being, always deliver the executables, without checking anything
#	if (($monthNow != $month)||($dayNow != $day)){
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
    my $msg = MIME::Lite->new(
			      From => 'Jean-Luc.Nougaret@cern.ch',
			      ReplyTo => 'Jean-Luc.Nougaret@cern.ch',
			      To => 'Jean-Luc.Nougaret@cern.ch',
			      Subject => 'MAD-X for Windows updated',
			      Data => "Dear colleagues,\n\nPlease take note that MAD-X version xxx is now available on Windows.\n\nThe new releases are available for download on the usual Web page:\nhttp://cern.ch/project-madwindows/MAD-X/\n\nRegards,\nJean-Luc"
			      );
    $msg->send;
    
    chdir($initialDir);

} # subroutine
