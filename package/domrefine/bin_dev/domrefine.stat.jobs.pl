#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM -b BEGIN_TIME -e END_TIME -u USER
-a: print all jobs of the user
-U: extract for all users
-v: for degug (print internal command and exit)
";

my %OPT;
getopts('b:e:u:Uav', \%OPT);

### Settings ###
my @KEY = ('jobname', 'jobnumber', 'hostname'
	   , 'owner'
	   # , 'qsub_time', 'start_time', 'end_time'
	   , 'ru_wallclock'
	   # , 'ru_utime', 'ru_stime'
	   , 'cpu', 'mem', 'io'
	   # , 'iow'
	   , 'maxvmem', 'failed');

my $COMMAND = "qacct -f /mnt/nfs/sge_account/accounting -j";
if ($OPT{b}) {
    my $begin_time = format_time($OPT{b});
    $COMMAND .= " -b $begin_time";
} else {
    die $USAGE;
}
if ($OPT{e}) {
    my $end_time = format_time($OPT{e});
    $COMMAND .= " -e $end_time";
}

if ($OPT{U}) {
    # all user
} else {
    my $user;
    if ($OPT{u}) {
	$user = $OPT{u};
    } elsif ($ENV{USER}) {
	$user = $ENV{USER};
    } else {
	print STDERR "WARNING: USER is not set\n";
    }
    if ($user) {
	$COMMAND .= " -u $user";
    }
}

### Main ###
if ($OPT{v}) {
    print STDERR "$COMMAND\n";
    system "$COMMAND";
    exit;
}

my @LINE = `$COMMAND`;
chomp(@LINE);

print join("\t", @KEY), "\n";
my %VALUE = ();
my $COUNT = 0;
for my $line (@LINE) {
    if ($line =~ /^(\w+)\s+(.*\S) *$/) {
	$VALUE{$1} = $2;
    } elsif ($line =~ /^=+$/) {
	$COUNT++;
	if ($COUNT >= 2) {
	    print_values(\%VALUE, %OPT);
	}
    } else {
	die;
    }
}

################################################################################
### Functions ##################################################################
################################################################################
sub print_values {
    my ($r_value, %opt) = @_;

    if (! %{$r_value}) {
	print STDERR "no value to print\n";
	return;
    }

    if (! $opt{a}) {
	if (defined ${$r_value}{failed}) {
	    if (${$r_value}{failed} eq '0') {
		return;
	    }
	} else {
	    print STDERR "no value for key failed\n";
	}
    }

    my @out = ();
    for my $key (@KEY) {
	if (defined ${$r_value}{$key}) {
	    push @out, ${$r_value}{$key};
	} else {
	    print STDERR "no value for key $key\n";
	}
    }
    %{$r_value} = ();

    print join("\t", @out), "\n";
}

sub format_time {
    my ($time) = @_;

    if ($time =~ /^\d{4}$/) {
	my $date = `date +'%Y%m%d'`;
	chomp($date);
	$time = "$date$time";
    } elsif ($time =~ /^\d{8}$/) {
	# MMDDhhmm
    } elsif ($time =~ /^\d{10}$/) {
	# YYMMDDhhmm
    } elsif ($time =~ /^\d{12}$/) {
	# CCYYMMDDhhmm
    }

    return $time;
}
