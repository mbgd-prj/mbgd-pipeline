#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM [-c CLUSTER] DIR
";

use DomRefine::Read;

### Arguments ###
my %OPT;
getopts('c:', \%OPT);

if (@ARGV != 1) {
    print STDERR $USAGE;
    exit 1;
}
my ($DIR) = @ARGV;
if (! -d $DIR) {
    print STDERR "Cannot find directory: $DIR\n";
    exit 1;
}

my $CLUSTER = $DIR;
if ($CLUSTER =~ /^(\S+)\.\S+$/) {
    $CLUSTER = $1;
}
if ($OPT{c}) {
    $CLUSTER = $OPT{c};
}
if (! -f $CLUSTER) {
    print STDERR "Cannot find cluster: $CLUSTER\n";
    exit 1;
}

### Main ###
my @CLUSTER = get_clusters($CLUSTER);

for my $cluster (@CLUSTER) {
    if (-f "$DIR/$cluster.out") {
	if (! -s  "$DIR/$cluster.out") {
	    print "WARNING: empty file for cluster $cluster\n";
	}
    } else {
	print "ERROR: no file for cluster $cluster\n";
    }
}
