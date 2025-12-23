#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: cat LINKS | $PROGRAM [cluster_set_focused]
";

my %OPT;
getopts('', \%OPT);

if (@ARGV == 0) {
    system "cat";
    exit;
}
if (@ARGV > 1) {
    print STDERR $USAGE;
    exit 1;
}

my ($FILE) = @ARGV;
my %CHECK = ();
open(FILE, $FILE) || die;
while (<FILE>) {
    chomp;
    my ($n, $clusterset) = split(/\s+/, $_);
    my @cluster = split(",", $clusterset);
    if (@cluster == 0) {
	die;
    }
    $CHECK{$cluster[0]} = 1;
    # for my $cluster (@cluster) {
    # 	$CHECK{$cluster} = 1;
    # }
}
close(FILE);

while (<STDIN>) {
    chomp;
    my @cluster = split(/[-\s,]/, $_);
    if (@cluster != 2) {
	die;
    }
    if ($CHECK{$cluster[0]} || $CHECK{$cluster[1]}) {
	print $_, "\n";
    }
}
