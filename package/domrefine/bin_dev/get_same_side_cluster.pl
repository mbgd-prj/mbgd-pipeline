#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM
";

my %OPT;
getopts('', \%OPT);

!@ARGV && -t and die $USAGE;
my %RIGHT_CLUSTER = ();
my %LEFT_CLUSTER = ();
my %CLUSTER = ();
while (<>) {
    my ($cluster1, $cluster2) = split;
    $RIGHT_CLUSTER{$cluster1}{$cluster2} = 1;
    $LEFT_CLUSTER{$cluster2}{$cluster1} = 1;
    $CLUSTER{$cluster1} = 1;
    $CLUSTER{$cluster2} = 1;
}

my %SAME_SIDE = ();
for my $cluster (keys %CLUSTER) {
    for my $right_cluster (keys %{$RIGHT_CLUSTER{$cluster}}) {
	for my $left_cluster (keys %{$LEFT_CLUSTER{$right_cluster}}) {
	    if ($left_cluster > $cluster) {
		$SAME_SIDE{$cluster}{$left_cluster} = 1;
		print STDERR "$right_cluster\t$cluster\t$left_cluster\n";
	    }
	}
    }
    for my $left_cluster (keys %{$LEFT_CLUSTER{$cluster}}) {
	for my $right_cluster (keys %{$RIGHT_CLUSTER{$left_cluster}}) {
	    if ($right_cluster > $cluster) {
		$SAME_SIDE{$cluster}{$right_cluster} = 1;
		print STDERR "$left_cluster\t$cluster\t$right_cluster\n";
	    }
	}
    }
}

for my $cluster1 (sort {$a<=>$b} keys %SAME_SIDE) {
    for my $cluster2 (sort {$a<=>$b} keys %{$SAME_SIDE{$cluster1}}) {
	print "$cluster1\t$cluster2\n";
    }
}
