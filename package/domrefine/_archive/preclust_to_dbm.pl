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

my %PRECLUST = ();
dbmopen(%PRECLUST, "preclust.dbm", 0644) || die; # use DBM

my @arr = ();
while (<>) {
    chomp;
    if (/^$/) {
        if (@arr) {
            $PRECLUST{$arr[0]} = "@arr";
        } else {
            die;
        }
        @arr = ();
        next;
    }
    my @f = split(/\t/, $_);
    if (@f == 6) {
        if ($f[0] =~ /^ +(\w+:\w+)$/) {
            my $gene = $1;
            push(@arr, $gene);
        }
    } elsif (@f == 1) {
        if (/^\* (\w+:\w+)$/) {
            my $rep = $1;
            if (@arr) {
                die;
            }
            push @arr, $rep;
        } else {
            die;
        }
    } else {
        die $_;
    }
}
if (@arr) {
    $PRECLUST{$arr[0]} = "@arr";
}

dbmclose(%PRECLUST);
