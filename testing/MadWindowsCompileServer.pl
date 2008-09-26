#!/usr/bin/perl

# wait until woken-up by reception of a "Compile MAD for Windows!\n" message
# on the receiving socket

# this script should be scheduled on each startup of the Windows box, with the
# Start->All Programs->Accessories->System Tools->Scheduled Tasks.
# 
# the script is installed under /afs/user/n/nougaret/scractch0/mad-automation
# together with the other scripts, BUT IF I DO SO IT CAN'T BE VIEWED THROUGH SAMBA
# HENCE NEED TO FIND SOMETHING ELSE...


use IO::Socket;
use Cwd; # to get current directory on Windows (would work on Linux as well)

my $os = 'windows';

my $socketPortWindows = '800'; # agreed-up with client
my $socketPortLinux = '801'; # could be the same as above

my $windowsHost = 'abpc10788';

my $sock = new IO::Socket::INET ( 
				  LocalHost => $windowsHost,
				  LocalPort => $socketPortWindows,
				  Proto => 'tcp',
				  Listen => 1,
				  Reuse => 1
				  ); 

die "Could not create client socket: $!\n" unless $sock; 


print "now accepting messages sent throught socket $socketPortWindows\n";

my $receiving = $sock->accept();

my $clientHost;

while (<$receiving>){
    print $_;
    if (/([\w\d\-\_]+) asks: Compile MAD for Windows!/){
	
	# following is specific to the DOS system
	if ($os eq 'windows'){
	    my $dir = getcwd; # get current working directory
	    print "we are now in $dir\n";
	    chdir("Z:\\MAD-X-WINDOWS\\madX");
            $dir = getcwd;
	    print "we are now in $dir\n";
	    # first make sure we are in sync with the CVS repository
	    my @needsUpdate = `cvs status | grep Need`;
	    print "list of files requiring refreshment from the CVS:\n";
	    foreach $file (@needsUpdate){
		print "$file\n";
	    }
	    # if is cumbersome on Windows, we may assume the CVS synch is already achieved
	    # on the Linux side.


	    print "now trigger compilation of madx.exe, madxp.exe and pars.exe\n";
	    # invoke compilation of madx.exe, madxp.exe and pars.exe
	    system("Makefile.bat");
	    # moving these files to the proper AFS web folder is left to the calling
	    # Linux-side Perl script...
	    # (an e-mail will be sent to observers of the Windows distribution of MAD-X)

	    # chdir("Z:\\MAD-X\\madX\");
	}


	$clientHost = $1; 
	print "clientHost is '$clientHost'\n";

       	my $replySock = new IO::Socket::INET(
					      PeerAddr => $clientHost,
					      PeerPort => $socketPortLinux,
					      Proto => 'tcp'
					     );

	die "Could not create socket $socketPortLinux to reply: $!\n" unless $replySock;
	print $replySock "Compilation completed\n";
	print "sent acknowledgement to $clientHost\n";	
	close($replySock);

	# once compliation completed, send an acknowledgement message to the
	# the Linux box that triggered the compilation in the first place

	#last;

    }
} # should stay here for ever...




close($sock);
