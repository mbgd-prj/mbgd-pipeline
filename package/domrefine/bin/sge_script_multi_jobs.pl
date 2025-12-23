#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: cat CLUSTER_SET | $PROGRAM -i INPUT -o DIR [OPTION] SCRIPT
-q QUEUE (by default, small)
-r : resume
";

use DomRefine::General;

### Settings ###
my %OPT;
getopts('i:o:q:r', \%OPT);

if (@ARGV != 1) {
    print STDERR $USAGE;
    exit 1;
}
my ($SCRIPT) = @ARGV;

my $DIR = $ENV{PWD};

my $INPUT;
if ($OPT{i}) {
    $INPUT = "$DIR/$OPT{i}";
    if (! -s $INPUT) {
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

my @COMMAND = ();
while (my $line = <STDIN>) {
    chomp($line);
    my @cluster = split(/\s+/, $line);
    my $command = "";
    for my $cluster (@cluster) {
	if (! $OPT{r} || ! -s "$OUT_DIR/$cluster.out") {
	    $command .= "cat $INPUT | dom_extract $cluster | $SCRIPT > $OUT_DIR/$cluster.out 2> $OUT_DIR/$cluster.log ;";
	}
    }
    if ($command ne "") {
	push @COMMAND, $command;
    }
}

### Submit ###
my $SUBMITTED_JOBS = "$OUT_DIR.jobs";
open(SUBMITTED_JOBS, ">$SUBMITTED_JOBS") || die;
my $COUNT = 0;
for my $command (@COMMAND) {
    if ($COUNT % $ENV{DOMREFINE_QSUB_UNIT} == 0) {
	check_queue();
    }
    $COUNT ++;
    print SUBMITTED_JOBS `sge.pl $QUEUE -N mj$COUNT '$command'`;
}
close(SUBMITTED_JOBS);

system "sge_check_jobs.pl -N -n0 $SUBMITTED_JOBS > $SUBMITTED_JOBS.check";
