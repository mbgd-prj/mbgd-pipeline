#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: cat CLUSTER_LINK | $PROGRAM -i FILE_OR_DIR -o DIR [OPTION] SCRIPT
-q QUEUE (by default, small)
-r : forced resume
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

my $INPUT;
if ($OPT{i}) {
    $INPUT = "$ENV{PWD}/$OPT{i}";
    if (! -d $INPUT) {
	die $INPUT;
    }
} else {
    die;
}
my $OUT_DIR;
if ($OPT{o}) {
    $OUT_DIR = "$ENV{PWD}/$OPT{o}";
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
my $COUNT = 0;
while (my $cluster_set = <STDIN>) {
    $cluster_set =~ s/^(\S+?)\s.*/$1/;
    if (! -s "$INPUT/$cluster_set.out") {
	next;
    }
    my $j_boundary = `cat $INPUT/$cluster_set.log | grep j_max_score | cut.sh 1 | cut.sh 2 =`;
    chomp($j_boundary);
    if ($j_boundary !~ /^\d+$/) {
	next;
    }

    my $command = "";
    if (! $OPT{r} || ! -s "$OUT_DIR/$cluster_set.out") {
	$command .= "cat $INPUT/$cluster_set.out | $SCRIPT $j_boundary > $OUT_DIR/$cluster_set.out 2> $OUT_DIR/$cluster_set.log ;";
    }
    if (! $OPT{r} || ! -s "$OUT_DIR/score.after.$cluster_set.out") {
	$command .= "cat $OUT_DIR/$cluster_set.out | $SCORE > $OUT_DIR/score.after.$cluster_set.out 2> $OUT_DIR/score.after.$cluster_set.log";
    }
    if ($command ne "") {
	if ($COUNT % $ENV{DOMREFINE_QSUB_UNIT} == 0) {
	    check_queue();
	}
	my $submitted_job = `sge.pl $QUEUE -N s$cluster_set '$command'`;
	print SUBMITTED_JOBS $submitted_job;
	$COUNT ++;
    }
}
close(SUBMITTED_JOBS);

system "sge_check_jobs.pl -N -n0 $SUBMITTED_JOBS > $SUBMITTED_JOBS.check";
