#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM CLUSTER_FILE

Modify domain organizations (merge_divide_tree)

Intermediate output files:
CLUSTER_FILE.link
CLUSTER_FILE.merge_divide_test(|.jobs|.jobs.check|.summary|.to_divide)

Final output files:
CLUSTER_FILE.merge_divide(|.log)
";

my %OPT;
getopts('', \%OPT);

if (@ARGV == 0) {
    print STDERR $USAGE;
    exit 1;
}
my ($CLUSTER) = @ARGV;

### Settings ###
my $MERGE_DIVIDE = 'dom_merge_divide_tree';
$MERGE_DIVIDE .= ' -R'; # align only sub-sequences in specific regions

my $SCORE = "dom_score_c";
$SCORE .= " -r"; # align full sequences using specific regions

### Pre ###
print STDERR "\n";
print STDERR '@' . `date`;
print STDERR ">$PROGRAM @ARGV\n";

if (! -e "$CLUSTER.link") {
    my $start_time = time;
    system "cat $CLUSTER | dom_network -l > $CLUSTER.link";
    my $end_time = time;
    printf STDERR "pre:\t%.2f\tmin\n", ($end_time - $start_time)/60;
}

### Main ###
my $GET_CLUSTER = "cat $CLUSTER.link | ignore_clusters.pl cluster.to_ignore";
system "$GET_CLUSTER | sge_script_score.pl -r -i $CLUSTER -o $CLUSTER.merge_divide_test '$MERGE_DIVIDE' '$SCORE'";

### Post ###
{
    my $start_time = time;
    system "$GET_CLUSTER | summarize_score.pl -i $CLUSTER.merge_divide_test/score.before. -o $CLUSTER.merge_divide_test/score.after. > $CLUSTER.merge_divide_test.summary";
    system "cat $CLUSTER.merge_divide_test.summary | sort -k4,4gr | select_positive.pl -f 4 > $CLUSTER.merge_divide_test.to_divide";
    system "cat $CLUSTER.merge_divide_test.to_divide | clusterset_modify -i $CLUSTER -o $CLUSTER.merge_divide_test > $CLUSTER.merge_divide.out 2> $CLUSTER.merge_divide.log";
    if ($ENV{DOMREFINE_ALLOW_SPACING}) {
        system "cat $CLUSTER.merge_divide.out | dom_renumber > $CLUSTER.merge_divide";
    } else {
        system "cat $CLUSTER.merge_divide.out | dom_concat_domains > $CLUSTER.merge_divide";
    }
    my $end_time = time;
    printf STDERR "post:\t%.2f\tmin\n", ($end_time - $start_time)/60;
}
