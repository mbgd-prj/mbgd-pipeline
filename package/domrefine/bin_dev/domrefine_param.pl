#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM ORIGINAL_CLUSTER_FILE SEQ_DB
creates directories for vaious parameters
";

my %OPT;
getopts('r', \%OPT);

if (@ARGV != 2) {
    print STDERR $USAGE;
    exit 1;
}
my ($ORIGINAL, $SEQ_DB) = @ARGV;
if ($ORIGINAL !~ /^\//) {
    $ORIGINAL = "$ENV{PWD}/$ORIGINAL";
}

# my @GAP_OPEN = (10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
# my @THRESHOLD = (-0.1, -0.05, 0);
# my @GAP_EXT = (0.5, 0.4, 0.3, 0.2, 0.1, 0);

my @GAP_OPEN = (7, 6, 5, 4, 3, 2, 1, 0);
my @THRESHOLD = (-0.05);
my @GAP_EXT = (0.5);

for my $o (@GAP_OPEN) {
    for my $t (@THRESHOLD) {
	for my $e (@GAP_EXT) {
	    my $dir = "t${t}e${e}o${o}";
	    if ($OPT{r}) {
		if (-e $dir) {
		    next;
		}
	    }
	    mkdir $dir;
	    chdir $dir;
	    system "domrefine.pl -t $t -e $e -o $o $ORIGINAL $SEQ_DB";
	    chdir "..";
	}
    }
}
