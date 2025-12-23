#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: cat CLUSTER_LINK | $PROGRAM -i INPUT_CLUSTER -o DIR [OPTION] 'SCRIPT' 'SCORE'
-q QUEUE (by default, small)
-r : resume
";

use DomRefine::General;

### Settings ###
my %OPT;
getopts('i:o:q:r', \%OPT);

if (@ARGV != 2) {
    print STDERR $USAGE;
    exit 1;
}
my ($SCRIPT, $SCORE) = @ARGV;

my $DIR = $ENV{PWD};

my $INPUT_CLUSTER;
if ($OPT{i}) {
    $INPUT_CLUSTER = "$DIR/$OPT{i}";
    if (! -s $INPUT_CLUSTER) {
	die;
    }
} else {
    die;
}

my $OUT_DIR;
if ($OPT{o}) {
    $OUT_DIR = "$DIR/$OPT{o}";
} else {
    die $USAGE;
}

my $QUEUE = "";
if ($OPT{q}) {
    $QUEUE = "-q $OPT{q}";
}

if (! -e $OUT_DIR) {
    system "mkdir $OUT_DIR";
}

### Main ###
-t and die $USAGE;

my $SUBMITTED_JOBS = "$OUT_DIR.jobs";
open(SUBMITTED_JOBS, ">$SUBMITTED_JOBS") || die;
my $count = 0;
while (my $line = <STDIN>) {
    chomp($line);
    my ($cluster_set) = split(/\s+/, $line);
    my $command = "";
    if (! $OPT{r} || ! -s "$OUT_DIR/score.before.$cluster_set.out") {
	$command .= "cat $INPUT_CLUSTER | dom_extract $cluster_set | $SCORE > $OUT_DIR/score.before.$cluster_set.out 2> $OUT_DIR/score.before.$cluster_set.log ;";
    }
    if (! $OPT{r} || ! -s "$OUT_DIR/$cluster_set.out") {
	$command .= "cat $INPUT_CLUSTER | dom_extract $cluster_set | $SCRIPT > $OUT_DIR/$cluster_set.out 2> $OUT_DIR/$cluster_set.log ;";
    }
    if (! $OPT{r} || ! -s "$OUT_DIR/score.after.$cluster_set.out") {
	$command .= "cat $OUT_DIR/$cluster_set.out | $SCORE > $OUT_DIR/score.after.$cluster_set.out 2> $OUT_DIR/score.after.$cluster_set.log";
    }
    if ($command ne "") {
	if ($count % $ENV{DOMREFINE_QSUB_UNIT} == 0) { # periodically, check queue
	    check_queue();
	}
	my $submitted_job = `sge.pl $QUEUE -N md$cluster_set '$command'`;
	print SUBMITTED_JOBS $submitted_job;
	$count ++;
    }
}
close(SUBMITTED_JOBS);

system "sge_check_jobs.pl -N -n0 $SUBMITTED_JOBS > $SUBMITTED_JOBS.check";
