#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: cat CLUSTER_SET | $PROGRAM -i INPUT -o DIR [OPTION] SCRIPT
-p PREFIX: prefix of output files in DIR
-q QUEUE (by default, small)
";

use DomRefine::General;

### Settings ###
my %OPT;
getopts('i:o:p:q:', \%OPT);

if (@ARGV != 1) {
    print STDERR $USAGE;
    exit 1;
}
my ($SCRIPT) = @ARGV;

my $DIR = $ENV{PWD};

my $IN_DIR;
if ($OPT{i}) {
    $IN_DIR = "$DIR/$OPT{i}";
} else {
    die;
}
if (! -d $IN_DIR) {
    die;
}

my $OUT_DIR;
if ($OPT{o}) {
    $OUT_DIR = "$DIR/$OPT{o}";
} else {
    die $USAGE;
}
if (! -e $OUT_DIR) {
    system "mkdir $OUT_DIR";
}

my $QUEUE = "";
if ($OPT{q}) {
    $QUEUE = "-q $OPT{q}";
}

my $PREFIX = "";
if ($OPT{p}) {
    $PREFIX = "$OPT{p}.";
}

### Main ###
-t and die $USAGE;

my $SUBMITTED_JOBS = "$OUT_DIR.${PREFIX}jobs";
open(SUBMITTED_JOBS, ">$SUBMITTED_JOBS") || die;
my $COUNT = 0;
while (my $line = <STDIN>) {
    chomp($line);
    my ($cluster) = $line;
    if (! -s "$OUT_DIR/$PREFIX$cluster.fa") {
	if ($COUNT % $ENV{DOMREFINE_QSUB_UNIT} == 0) {
	    check_queue();
	}
	my $submitted_job = `sge.pl $QUEUE -N s$cluster 'echo "$cluster" 1>&2; $SCRIPT $IN_DIR/$cluster.fa > $OUT_DIR/$PREFIX${cluster}.fa 2> $OUT_DIR/$PREFIX${cluster}.log'`;
	print SUBMITTED_JOBS $submitted_job;
	$COUNT ++;
    }
}
close(SUBMITTED_JOBS);

system "sge_check_jobs.pl -N -n0 $SUBMITTED_JOBS > $SUBMITTED_JOBS.check";
