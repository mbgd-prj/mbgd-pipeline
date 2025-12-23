#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM
";

my %OPT;
getopts('', \%OPT);

### Main ###
!@ARGV && -t and die $USAGE;

my $CLUSTER = "";
while (<>) {
    chomp;
    if (/^#/) {
        next;
    }
    if (/^Cluster (\S+)/) {
	$CLUSTER = $1;
    } elsif (/^(\S+:\S+) (\-?\d+) (\-?\d+)$/) {
	my ($gene, $start, $end) = ($1, $2, $3);
	my $domain = 0;
	if ($gene =~ /(\S+)\((\d+)\)$/) {
	    $gene = $1;
	    $domain = $2;
	}
	print "$CLUSTER $gene $domain $start $end\n";
    } elsif (/^$/) {
    } elsif (/^(Hom|Sub)Cluster (\S+)/) {
    } else {
	die $_;
    }
}
