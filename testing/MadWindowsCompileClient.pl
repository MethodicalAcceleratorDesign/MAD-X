#!/usr/bin/perl

# assume a Perl script in charge of compiling MAD on the Windows platform
# is waiting on an agreed-upon socket => the Server must be started before
# the Client, if not the program should return an error otherwise it would
# wait for ever!

# output of this program on stdout: SUCCESS or FAILURE:<message>

my $windowsHost = 'abpc10788';


use IO::Socket;

use Sys::Hostname;


$socketPortWindows = '800'; # agreed-up with client (<1024 for non-root)
$socketPortLinux = '801'; # could be the same as above

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
					 localHost => $thisLinuxHost,
					 localPort => $socketPortLinux,
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
	print "=> should now install the executables in the AFS web folder\n";
	last; # leave the while loop
    }
}

close ($clientSock);

sub checkWindowsCompilationOutcome {
    my $initialDir = `pwd`;
    # check the delivery directory contents
    chdir($initialDir);

}
