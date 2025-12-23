#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: cat link | $PROGRAM
";

use DomRefine::General;
use DomRefine::Refine;

my %OPT;
getopts('mn', \%OPT);

!@ARGV && -t and die $USAGE;

my %LINK = ();
while (<>) {
    chomp;
    my @cluster = split(/[-,\s]/, $_);
    if (@cluster != 2) {
	die;
    }
    my ($cluster1, $cluster2) = @cluster;
    $LINK{$cluster1}{$cluster2} = 1;
    $LINK{$cluster2}{$cluster1} = 1;
}

for my $cluster (keys %LINK) {
    my @cluster = keys %{$LINK{$cluster}};
    print $cluster, "\t", scalar(@cluster), "\n";
}
