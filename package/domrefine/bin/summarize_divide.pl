#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: cat CLUSTER_LIST | $PROGRAM DIR
";

use DomRefine::Read;

### Settings ###
my %OPT;
getopts('', \%OPT);

if (@ARGV != 1) {
    print STDERR $USAGE;
    exit 1;
}
my ($DIR) = $ARGV[0];

### Main ###
-t and die $USAGE;
while (my $c = <STDIN>) {
    chomp($c);
    if (-s "$DIR/$c.log") {
	my ($i, $result);
	open(LOG, "$DIR/$c.log") || die;
	while (my $line = <LOG>) {
	    chomp($line);
	    if ($line =~ /^i=(\d+)\S+(.*)/) {
		($i, $result) = ($1, $2);
		last;
	    }
	}
	close(LOG);
	if (defined $i) {
	    print join("\t", $c,"i$i", $result), "\n";
	} else {
	    print "\n";
	}
    }
}
