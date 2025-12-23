#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM CACHE_DIR
";

my %OPT;
getopts('', \%OPT);

if (@ARGV == 0) {
    print STDERR $USAGE;
    exit 1;
}
my ($CACHE_DIR) = @ARGV;

my @HEX_CHAR = (0..9,'a','b','c','d','e','f');

my $DATE = `date`;

### Main ###
if (! -e $CACHE_DIR) {
    mkdir $CACHE_DIR;
}
chdir $CACHE_DIR || die;

{
    my $count = 0;
    my $start_time = time;
    for my $char1 (@HEX_CHAR) {
	for my $char2 (@HEX_CHAR) {
	    if (! -e "$char1$char2") {
		mkdir "$char1$char2";
		$count ++;
	    }
	}
    }
    my $end_time = time;
    my $diff_time = $end_time - $start_time;

    if ($count != 0) {
	print STDERR "\n";
	print STDERR '@' . $DATE;
	printf STDERR "prepare_cache:\t%.2f\tmin\n", $diff_time/60;
    }
}
