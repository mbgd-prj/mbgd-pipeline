#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM -f N
";

my %OPT;
getopts('f:', \%OPT);

my $F = $OPT{f} || die $USAGE;
my $I = $F - 1;

!@ARGV && -t and die $USAGE;
while (<>) {
    my @x = split;
    if ($x[$I] > 0) {
	print;
    }
}
