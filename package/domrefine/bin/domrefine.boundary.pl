#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM [OPTION] CLUSTER_FILE
-e GAP_EXT
-o GAP_OPEN

Modify domain organizations: move_boundary, create_boundary

Intermediate output files:
CLUSTER_FILE.link(|_eval|_eval.score.jobs|_eval.score.jobs.check)
CLUSTER_FILE.move_test(|.jobs|.jobs.check|.score.jobs|.score.jobs.check|.summary)
CLUSTER_FILE.move.create_test(|.jobs|.jobs.check|.score.jobs|.score.jobs.check|.summary)

Final output files:
CLUSTER_FILE.move(|.log)
CLUSTER_FILE.move.create(|.log)
";

print STDERR "\n";
print STDERR '@' . `date`;
print STDERR ">$PROGRAM @ARGV\n";

my %OPT;
getopts('t:e:o:', \%OPT);

### Settings ###
my $OPTION = "";
if (defined $OPT{e}) {
    $OPTION .= " -e $OPT{e}";
}
if (defined $OPT{o}) {
    $OPTION .= " -o $OPT{o}";
}

my $SCORE = "dom_score_c";
# $SCORE .= " -r"; # align by region
$SCORE .= $OPTION;

my $MOVE_BOUNDARY = "dom_move_boundary";
# $MOVE_BOUNDARY .= " -r"; # align by region
$MOVE_BOUNDARY .= $OPTION;

my $CREATE_BOUNDARY = "dom_create_boundary";
# $CREATE_BOUNDARY .= " -r"; # align by region
$CREATE_BOUNDARY .= $OPTION;

### Main ###
if (@ARGV == 0) {
    print STDERR $USAGE;
    exit 1;
}
my ($CLUSTER) = @ARGV;

my $GET_CLUSTER = "cat $CLUSTER.link | ignore_clusters.pl cluster.to_ignore";

### prepare ###
my $START_TIME = time;
if (! -f "$CLUSTER.link") {
    system "cat $CLUSTER | dom_network -l > $CLUSTER.link";
}
my $END_TIME = time;
printf STDERR "pre:\t%.2f\tmin\n", ($END_TIME - $START_TIME)/60;

### move boundary ###
if (! -f "$CLUSTER.move") {
    my $start_time = time;
    system "$GET_CLUSTER | sge_script_score.pl -r -i $CLUSTER -o $CLUSTER.move_test '$MOVE_BOUNDARY' '$SCORE'";
    my $start_time2 = time;
    system "$GET_CLUSTER | summarize_score.pl -i $CLUSTER.move_test/score.before. -o $CLUSTER.move_test/score.after. > $CLUSTER.move_test.summary";
    system "cat $CLUSTER.move_test.summary | select_positive.pl -f 4 | sort -k4,4gr > $CLUSTER.move_test.to_move";
    system "cat $CLUSTER.move_test.to_move | clusterset_modify -i $CLUSTER -o $CLUSTER.move_test > $CLUSTER.move 2> $CLUSTER.move.log";
    system "cat $CLUSTER.move.log | grep modified_to >  $CLUSTER.moved";
    my $end_time = time;
    printf STDERR "post:\t%.2f\tmin\n", ($end_time - $start_time2)/60;
    printf STDERR "move:\t%.2f\tmin\n", ($end_time - $start_time)/60;
    printf STDERR "\n";
}

### create boundary ###
if (! -f "$CLUSTER.move.create") {
    if ($ENV{DOMREFINE_CREATE_BOUNDARY}) {
	my $start_time = time;
	system "cat $CLUSTER.moved | cut -f1 | sge_boundary_score.pl -r -i $CLUSTER.move_test -o $CLUSTER.move.create_test '$CREATE_BOUNDARY' '$SCORE'";
	my $start_time2 = time;
	system "cat $CLUSTER.moved | cut -f1 | summarize_score.pl -i $CLUSTER.move_test/score.after. -o $CLUSTER.move.create_test/score.after. > $CLUSTER.move.create_test.summary";
	system "cat $CLUSTER.move.create_test.summary | select_positive.pl -f 4 | sort -k4,4gr > $CLUSTER.move.create_test.to_create";
	system "cat $CLUSTER.move.create_test.to_create | clusterset_modify -i $CLUSTER.move -o $CLUSTER.move.create_test > $CLUSTER.move.create 2> $CLUSTER.move.create.log";
	my $end_time = time;
	printf STDERR "post:\t%.2f\tmin\n", ($end_time - $start_time2)/60;
	printf STDERR "create:\t%.2f\tmin\n", ($end_time - $start_time)/60;
	printf STDERR "\n";
    } else {
	system "ln -s $CLUSTER.move $CLUSTER.move.create";
    }
}
