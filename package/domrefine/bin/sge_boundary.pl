#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: cat CLUSTER_LINK | $PROGRAM -i FILE_OR_DIR -o DIR [OPTION] SCRIPT
-p PREFIX: prefix of output files in DIR
-q QUEUE (by default, small)
-r : resume
-R : forced resume
";

use DomRefine::General;

### Settings ###
my %OPT;
getopts('i:o:p:q:rR', \%OPT);

if (@ARGV != 1) {
    print STDERR $USAGE;
    exit 1;
}
my ($SCRIPT) = @ARGV;

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
my $PREFIX = "";
if ($OPT{p}) {
    $PREFIX = "$OPT{p}.";
}

if (! -e $OUT_DIR) {
    system "mkdir $OUT_DIR";
}

### Main ###
-t and die $USAGE;

my $SUBMITTED_JOBS = "$OUT_DIR.${PREFIX}jobs";
open(SUBMITTED_JOBS, ">$SUBMITTED_JOBS") || die;
my $count = 0;
while (my $cluster_set = <STDIN>) {
    $cluster_set =~ s/^(\S+?)\s.*/$1/;
    if (! -f "$INPUT/$cluster_set.out" or -z "$INPUT/$cluster_set.out") {
	next;
    }
    my $j_boundary = `cat $INPUT/$cluster_set.log | grep j_max_score | cut.sh 1 | cut.sh 2 =`;
    chomp($j_boundary);
    if ($j_boundary !~ /^\d+$/) {
	next;
    }
    my $out = "$OUT_DIR/$PREFIX$cluster_set";
    if ($OPT{r}) {
	if (-f "$out.out") {
	    next;
	}
    }
    if ($OPT{R}) {
	if (-f "$out.out" and ! -z "$out.out") {
	    next;
	}
    }
    if ($count % $ENV{DOMREFINE_QSUB_UNIT} == 0) {
	check_queue();
    }
    my $submitted_job = `sge.pl $QUEUE -N s$cluster_set 'cat $INPUT/$cluster_set.out | $SCRIPT $j_boundary > $out.out 2> $out.log'`;
    print SUBMITTED_JOBS $submitted_job;
    $count ++;
}
close(SUBMITTED_JOBS);

system "sge_check_jobs.pl -N -n0 $SUBMITTED_JOBS > $SUBMITTED_JOBS.check";
