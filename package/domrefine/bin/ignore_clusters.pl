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

### Read cluster to ignore
my %CHECK = ();
if (@ARGV != 0 and -f $ARGV[0]) {
    open(FILE, $ARGV[0]) || die;
    while (<FILE>) {
	chomp;
	my ($cluster) = split;
	$CHECK{$cluster} = 1;
    }
    close(FILE);
}

### Filter input clusterset
while (<STDIN>) {
    chomp;
    my @cluster = split(/[-\s,]/, $_);
    if (check_cluster(@cluster)) {
    } else {
    	print $_, "\n";
    }
}

################################################################################
### Functions ##################################################################
################################################################################
sub check_cluster {
    my @cluster = @_;

    for my $cluster (@cluster) {
	if ($CHECK{$cluster}) {
	    return 1;
	}
    }

    return 0;
}
