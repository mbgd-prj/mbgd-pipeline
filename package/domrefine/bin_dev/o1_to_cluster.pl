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

!@ARGV && -t and die $USAGE;

my $CLUSTER = "";
while (<>) {
    if (/^Cluster (\d+)$/) {
	$CLUSTER = $1;
    } elsif (/ (\S+) (\-?\d+ \-?\d+)$/) {
	my ($domain, $pos) = ($1, $2);
	my $domain_no = 0;
	if ($domain =~ /^(\S+)\((\d+)\)$/) {
	    ($domain, $domain_no) =  ($1, $2);
	}
	print "$CLUSTER $domain $domain_no $pos\n";
    }
}
