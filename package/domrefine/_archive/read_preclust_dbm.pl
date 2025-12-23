#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM [-k KEY] preclust.dbm
";

my %OPT;
getopts('k:', \%OPT);

if (!@ARGV) {
    print STDERR $USAGE;
    exit 1;
}
my ($FILE) = @ARGV;

my %PRECLUST = ();
dbmopen(%PRECLUST, $FILE, 0644) || die;
if ($OPT{k}) {
    my $key = $OPT{k};
    my $members = $PRECLUST{$key};
    print "$members\n";
} else {
    for my $key (keys(%PRECLUST)) {
        my $members = $PRECLUST{$key};
        print "$members\n";
    }
}
dbmclose(%PRECLUST);
