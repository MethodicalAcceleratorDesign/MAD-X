#!/usr/bin/perl

# this triggers remote invocation of the MAD-X compilation on Windows,
# through an intermediate machine on which the NFS directory seen through
# Samba is mounted.


use IO::Socket::INET;

my $intermediateNfsLinuxHost = 'cs-ccr-dev1'; # a.k.a abcopl1
my $intermediateSocketPort = 7075;

my $socket = new IO::Socket::INET(
				  PeerAddr => $intermediateNfsLinuxHost,
				  PeerPort => $intermediateSocketPort,
				  Proto => 'tcp'
				  );

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
