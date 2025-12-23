#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM [DCLST_FILE]
";

use DomRefine::Read;

my %OPT;
getopts('', \%OPT);

my $TMP_CLUSTER = define_tmp_file("$PROGRAM.cluster");
# my $TMP_HOMCLUSTER = define_tmp_file("$PROGRAM.homcluster");
END {
    remove_tmp_file($TMP_CLUSTER);
    # remove_tmp_file($TMP_HOMCLUSTER);
}

### Main ###
!@ARGV && -t and die $USAGE;

# open(CLUSTER, ">$TMP_CLUSTER") || die;
# open(HOMCLUSTER, ">$TMP_HOMCLUSTER") || die;
# open(HOMCLUSTER_CLUSTER, ">homcluster_cluster") || die;
# open(HOMCLUSTER_SCORE, ">homcluster.score") || die;
my $HOMCLUSTER = "";
my $CLUSTER = "";
while (<>) {
    chomp;
    if (/^HomCluster (\d+) : (\S+ \S+)/) {
	$HOMCLUSTER = $1;
	my $homcluster_score = $2;
	# print HOMCLUSTER_SCORE "$HOMCLUSTER\t$homcluster_score\n";
    } elsif (/^HomCluster/) {
	die;
    }
    if (/^Cluster (\d+)/) {
	$CLUSTER = $1;
	# print HOMCLUSTER_CLUSTER "$HOMCLUSTER\t$CLUSTER\n";
    }
    # if (/^L \d+ \S+ (\S+:\S+) (\d+) (\d+) (\d+)$/) {
    if (/^L \d+ -1 (\S+:\S+) (\d+) (\d+) (\d+)$/) {
	my ($gene, $begin, $end, $domain) = ($1, $2, $3, $4);
	print "$CLUSTER $gene $domain $begin $end\n";
	# print CLUSTER "$CLUSTER $gene $domain $begin $end\n";
	# print HOMCLUSTER "$HOMCLUSTER $gene $domain $begin $end\n";
    }
}
# close(CLUSTER);
# close(HOMCLUSTER);
# close(HOMCLUSTER_CLUSTER);
# close(HOMCLUSTER_SCORE);

# system "cat $TMP_CLUSTER | /bin/sort -k1,1n -k2,2 -k3,3n -k4,4n | uniq > cluster";
# system "cat $TMP_HOMCLUSTER | /bin/sort -k1,1n -k2,2 -k3,3n -k4,4n | uniq > homcluster";
