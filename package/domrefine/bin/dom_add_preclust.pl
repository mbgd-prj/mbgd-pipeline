#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: cat cluster | $PROGRAM PRECLUST
-d: debug
";

my %OPT;
getopts('d', \%OPT);

if (!@ARGV) {
    print STDERR $USAGE;
    exit 1;
}
my ($PRECLUST_FILE) = @ARGV;

my $REP;
my %PRECLUST = ();
open(PRECLUST_FILE, "$PRECLUST_FILE") || die "$!";
while (<PRECLUST_FILE>) {
    chomp;
    if (/^$/) {
        if ($OPT{d}) {
            my @arr = keys %{$PRECLUST{$REP}};
            print "$REP @arr\n";
        }
        next;
    }
    my @f = split(/\t/, $_);
    if (@f == 6) {
        if ($f[0] =~ /^ +(\w+):(\S+)$/) {
            my ($org, $gene_id) = ($1, $2);
            $PRECLUST{$REP}{$org} = 1;
        } else {
            die;
        }
    } elsif (@f == 1) {
        if (/^\* (\w+:\S+)$/) {
            $REP = $1;
            if ($PRECLUST{$REP}) {
                die;
            }
        } else {
            die;
        }
    } else {
        die $_;
    }
}
if ($OPT{d}) {
    my @arr = keys %{$PRECLUST{$REP}};
    print "$REP @arr\n";
}
close(PRECLUST_FILE);

while (<STDIN>) {
    chomp;
    my @f = split;
    if (@f != 5) {
        die;
    }
    my $gene = $f[1];
    my @arr = ();
    if ($PRECLUST{$gene}) {
        @arr = sort keys(%{$PRECLUST{$gene}});
    }
    print "@f\t@arr\n";
}
