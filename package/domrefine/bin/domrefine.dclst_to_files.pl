#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM [DCLST_FILE]
";

my %OPT;
getopts('', \%OPT);

### Main ###
!@ARGV && -t and die $USAGE;

open(CLUSTER, ">cluster") || die;
open(HOMCLUSTER, ">homcluster") || die;
open(HOMCLUSTER_CLUSTER, ">homcluster_cluster") || die;
open(DESCRIPTION, ">cluster.descr") || die;
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
	print HOMCLUSTER_CLUSTER "$HOMCLUSTER\t$CLUSTER\n";
    }
    if (/^\#CDescr\t(.*)$/){
	$DESCRIPTION = $1;
	print DESCRIPTION "$CLUSTER\t$DESCRIPTION\n";
    }
    if (/^(\S+:\S+)\t(\d+)\t(\d+)\t(\d+)/) {
	print CLUSTER "$CLUSTER $1 $2 $3 $4\n";
	print HOMCLUSTER "$HOMCLUSTER $1 $2 $3 $4\n";
    }
}
close(CLUSTER);
close(HOMCLUSTER);
close(HOMCLUSTER_CLUSTER);
close(DESCRIPTION);
