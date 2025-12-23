#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM
";

my %OPT;
getopts('ahH', \%OPT);

### Main ###
!@ARGV && -t and die $USAGE;

my $HOMCLUSTER = "";
my $CLUSTER = "";
my $DESCRIPTION = "";
while (<>) {
    chomp;
    if (/^HomCluster (\d+)/) {
	$HOMCLUSTER = $1;
    }
    if (/^Cluster (\d+)/) {
	$CLUSTER = $1;
	if ($OPT{H}) {
	    print "$HOMCLUSTER\t$CLUSTER\n";
	}
    }
    if (/^\#CDescr\t(.*)$/){
	$DESCRIPTION = $1;
	if ($OPT{a}) {
	    print "$CLUSTER\t$DESCRIPTION\n";
	}
    }
    if (/^(\S+:\S+)\t(\d+)\t(\d+)\t(\d+)$/) {
	if ($OPT{a}) {
	    next;
	}
	if ($OPT{H}) {
	    next;
	}
	if ($OPT{h}) {
	    print "$HOMCLUSTER $1 $2 $3 $4\n";
	} else {
	    print "$CLUSTER $1 $2 $3 $4\n";
	}
    }
}
