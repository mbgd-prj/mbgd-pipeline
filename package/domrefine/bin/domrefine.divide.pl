#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM CLUSTER_FILE [PREV_DIVIDE_LOG]
-1: only once (not recursive)

Divide each cluster into tighter ortholog groups using phyogenetic tree.

Intermediate output files:
CLUSTER_FILE*.divide_test(|.jobs|.jobs.watch|.summary|.to_divide)
CLUSTER_FILE*.divide(|.log)

Final output files:
CLUSTER_FILE.divided
";

print STDERR "\n";
print STDERR '@' . `date`;
print STDERR ">$PROGRAM @ARGV\n";

my %OPT;
getopts('1', \%OPT);

### Settings ###
my $QUEUE = "";

my $DIVIDE_TEST = "dom_divide_test";
# $DIVIDE_TEST .= " -r"; # align by region

### Main ###
if (@ARGV == 0) {
    print STDERR $USAGE;
    exit 1;
}
my ($CLUSTER, $PREV_DIVIDE_LOG) = @ARGV;
my $CLUSTER_ORG = $CLUSTER;

if ($OPT{1}) {
    divide($CLUSTER, $PREV_DIVIDE_LOG);
} else {
    my @created;
    do {
	($CLUSTER, $PREV_DIVIDE_LOG) = divide($CLUSTER, $PREV_DIVIDE_LOG);
	@created = `cat $PREV_DIVIDE_LOG | grep created`;
    } until (@created == 0);
    system "ln -s $CLUSTER ${CLUSTER_ORG}.divided";
}

################################################################################
### Function ###################################################################
################################################################################
sub divide {
    my ($cluster, $prev_divide_log) = @_;

    my $out_dir = $cluster . ".divide_test";
    
    my $start_time = time;
    system "cat $cluster | dom_cluster_size > $cluster.size";
    system "cat $cluster | cut.sh 1 | sort | uniq | ignore_clusters.pl cluster.to_ignore > $out_dir.to_check";
    if (defined $prev_divide_log) {
	system "cat $prev_divide_log | grep created | cut.sh 2 | ignore_clusters.pl cluster.to_ignore > $out_dir.to_check";
    }
    system "cat $out_dir.to_check | clusterset_sort_by_size -n 100 $cluster.size > $out_dir.to_check.multi_jobs";
    my $end_time = time;
    printf STDERR "pre:\t%.2f\tmin\n", ($end_time - $start_time)/60;

    # divide test
    # system "cat $out_dir.to_check | sge_script.pl $QUEUE -i $cluster -o $out_dir '$DIVIDE_TEST'";
    system "cat $out_dir.to_check.multi_jobs | sge_script_multi_jobs.pl $QUEUE -r -i $cluster -o $out_dir '$DIVIDE_TEST'";

    $start_time = time;
    # summarize
    system "cat $out_dir.to_check | summarize_divide.pl $out_dir > $out_dir.summary";
    # execute
    system "cat $out_dir.summary | perl -F\"\\t\" -lane '\$F[6]>=0.5 and print' > $out_dir.to_divide";
    system "cat $out_dir.to_divide | cluster_summary_divide -i $cluster -o $out_dir > $cluster.divide 2> $cluster.divide.log";

    my $prefix = $cluster;
    my $n = 2;
    if ($cluster =~ /^(\S*\D)(\d+)$/) {
	($prefix, $n) = ($1, $2);
	$n ++ ;
    }
    # my @created = `cat $cluster.divide.log | grep created`;
    # if (@created != 0) {
	system "ln -s $cluster.divide $prefix$n";
    # }

    $end_time = time;
    printf STDERR "post:\t%.2f\tmin\n", ($end_time - $start_time)/60;

    return ("$prefix$n", "$cluster.divide.log");
}
