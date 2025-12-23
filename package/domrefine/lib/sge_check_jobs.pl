#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM QSUB_LOG_FILE
-n NUMBER: 1 by default, 0 for no limit
-t INTERVAL: (sec)
-N: output job names
-J: output job IDs
-j: output the number of jobs
";
# -H: supress header

### Settings ###
my %OPT;
getopts('n:t:NJjH', \%OPT);

my $NUMBER_OF_TIMES = 1;
if (defined $OPT{n}) {
    if ($OPT{n} < 0) {
	die $USAGE;
    }
    $NUMBER_OF_TIMES = $OPT{n};
}

my $INTERVAL = $ENV{DOMREFINE_QUEUE_CHECK_INTERVAL};
if (defined $OPT{t}) {
    $INTERVAL = $OPT{t};
}

$| = 1; # set autoflush for STDOUT

if (@ARGV != 1) {
    print STDERR $USAGE;
    exit 1;
}
my @submitted_job = ();
for my $line (`cat $ARGV[0]`) {
    if ($line =~ /^Your job (\d+) / or $line =~ /^(\d+)\.\S+/) { # SGE or PBS
	my $job_id = $1;
	push @submitted_job, $job_id;
    }
}

### Main ###
my $COUNT = 0;
while (1) {
    ### qstat
    my @jobs_in_queue = ();
    my @ids_in_queue = ();
    my @names_in_queue = ();
    my @header = ();
    my @line = `qstat`;
    for my $line (@line) {
	if ($line =~ /^\s*(\d+)(\s+|\.)\S+\s+(\S+)\s/) { # SGE or PBS
	    my $job_id = $1;
	    my $job_name = $3;
	    if (grep {/^$job_id$/} @submitted_job) {
		push @jobs_in_queue, $line;
		push @ids_in_queue, $job_id;
		push @names_in_queue, $job_name;
	    }
	} else {
	    push @header, $line;
	}
    }
    if (@ids_in_queue == 0) {
	last;
    }

    ### print jobs
    if ($OPT{N}) {
	print "@names_in_queue\n";
    } elsif ($OPT{J}) {
	print "@ids_in_queue\n";
    } elsif ($OPT{j}) {
	print scalar(@ids_in_queue), " jobs in queue\n";
    } else {
	if (! $OPT{H}) {
	    print @header;
	}
	print @jobs_in_queue;
    }

    ### finish ?
    if ($NUMBER_OF_TIMES) {
	$COUNT++;
	if ($COUNT >= $NUMBER_OF_TIMES) {
	    last;
	}
    }

    sleep $INTERVAL;
}
