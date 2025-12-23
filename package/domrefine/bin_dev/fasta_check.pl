#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: cat FASTA_FILE | $PROGRAM
";

my %OPT;
getopts('', \%OPT);

!@ARGV && -t and die $USAGE;
my $header = '';
while (my $line = <STDIN>) {
    chomp($line);
    if ($line =~ /^>(\S+)/) {
        $header = $1;
        next;
    }
    if ($line =~ /[^-A-Z]/) {
        print $header, "\n";
    }
}
