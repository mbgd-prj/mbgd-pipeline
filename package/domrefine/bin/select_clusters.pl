#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: cat LINKS | $PROGRAM CLUSTER_LIST_FILE
";

my %OPT;
getopts('', \%OPT);

my %CHECK = ();
if (@ARGV != 0 and -f $ARGV[0]) {
    open(FILE, $ARGV[0]) || die;
    while (<FILE>) {
	chomp;
	my ($x) = split;
	my @cluster = split(/[-\s,]/, $x);
	for my $cluster (@cluster) {
	    $CHECK{$cluster} = 1;
	}
    }
    close(FILE);
}

while (<STDIN>) {
    chomp;
    my @cluster = split(/[-\s,]/, $_);
    if (check_cluster(@cluster)) {
    	print $_, "\n";
    }
}

sub check_cluster {
    my @cluster = @_;

    for my $cluster (@cluster) {
	if (! $CHECK{$cluster}) {
	    return 0;
	}
    }

    return 1;
}
