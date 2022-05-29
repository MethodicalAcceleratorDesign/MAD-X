#!/usr/bin/perl

# this triggers remote invocation of the MAD-X compilation on Windows,
# through an intermediate machine on which the NFS directory seen through
# Samba is mounted.


use IO::Socket::INET;

my $intermediateNfsLinuxHost = 'cs-ccr-dev1'; # a.k.a abcopl1
my $intermediateSocketPort = 7075;

# this value should be such as to allow the acrontab to recreate the
# MadWindowsCompileClient.pl process in case it just died after 10 days
# of so or running, after which the AFS and Kerberos tokens can no
# longer be refreshed and need to be recreated with a brand new process

my $socket;
my $minutes = 5;
TWO_TRIALS: for ($i=1; $i<$minutes; $i++) {
    $socket = new IO::Socket::INET(
				      PeerAddr => $intermediateNfsLinuxHost,
				      PeerPort => $intermediateSocketPort,
				      Proto => 'tcp',
				      );
    if ($socket) { last TWO_TRIALS; } else {
	sleep 60;
    }# exit the loop
}



# die "Could not create socket: $!\n" unless $socket;

if ($socket){
    print $socket "Trigger Windows compilation\n";
    print "success\n"; # output in case of success

    close($socket);
} else {
    # output in case of failure
    print "failure: failed to create socket $intermediateSocketPort\n";
}

# we don't wait the reply, as the intermediate Linux perl script will handle it
