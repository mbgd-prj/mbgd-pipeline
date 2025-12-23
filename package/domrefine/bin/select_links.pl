#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: cat LINKS | $PROGRAM CLUSTER_LIST_FILE
";

my %OPT;
getopts('', \%OPT);

my %CHECK = ();
if (@ARGV != 0 and -f $ARGV[0]) {
    open(FILE, $ARGV[0]) || die;
    while (<FILE>) {
	chomp;
	my ($link) = split;
	$CHECK{$link} = 1;
    }
    close(FILE);
}

while (<STDIN>) {
    chomp;
    my ($link) = split;
    if ($CHECK{$link}) {
    	print $_, "\n";
    }
}
