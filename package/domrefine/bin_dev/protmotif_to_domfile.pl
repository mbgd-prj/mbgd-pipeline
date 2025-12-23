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
while (<STDIN>) {
    chomp;
    my @f = split("\t", $_, -1);
    if ($f[0] =~ /^\d+$/) {
	my $sp = $f[2];
	my $name = $f[3];
	my $from = $f[4];
	my $to = $f[5];
	print "$sp:$name $from $to\n";
    } else {
	next;
    }
}
