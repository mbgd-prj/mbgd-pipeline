#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM QSUB_LOG_FILE
";

### settings
my %OPT;
getopts('', \%OPT);

if (@ARGV != 1) {
    print STDERR $USAGE;
    exit 1;
}
my $QSUB_LOG = shift @ARGV;

### list of jobs
my @SUBMITTED_JOB = ();
for my $line (`cat $QSUB_LOG`) {
    if ($line =~ /^Your job (\d+) / or $line =~ /^(\d+)\.\S+/) { # SGE or PBS
	my $job_id = $1;
	push @SUBMITTED_JOB, $job_id;
    }
}

my @JOB_TO_DELETE = ();
for my $line (`qstat`) {
    if ($line =~ /^\s*(\d+)/) {
	my $job_id = $1;
	if (grep {/^$job_id$/} @SUBMITTED_JOB) {
	    push @JOB_TO_DELETE, $job_id;
	}
    }
}

### execute
if (@JOB_TO_DELETE) {
    system "qdel @JOB_TO_DELETE";
}
