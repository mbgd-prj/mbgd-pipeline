#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: cat CLUSTER_LINK | $PROGRAM -i REF_PREFIX -o OUT_PREFIX
";

use DomRefine::General;

### Settings ###
my %OPT;
getopts('i:o:', \%OPT);

my $DIR = $ENV{PWD};
my $OUT_DIR;
my $REF_DIR;
if ($OPT{i}) {
    $REF_DIR = "$DIR/$OPT{i}";
} else {
    die;
}
if ($OPT{o}) {
    $OUT_DIR = "$DIR/$OPT{o}";
} else {
    die;
}

my $SUFFIX = "out";

### Main ###
-t and die $USAGE;
my @cluster_pair = <STDIN>;
chomp(@cluster_pair);
for (my $i=0; $i<@cluster_pair; $i++) {
    my $ref_score;
    if (-f "$REF_DIR$cluster_pair[$i].out") {
	my $line = `cat $REF_DIR$cluster_pair[$i].out`;
	chomp($line);
	if ($line) {
	    ($ref_score) = split(/\s+/, $line);
	}
    }
    my ($score, $n_pos, $n_aa);
    if (-f "$OUT_DIR$cluster_pair[$i].$SUFFIX") {
	my $line = `cat $OUT_DIR$cluster_pair[$i].$SUFFIX`;
	chomp($line);
	if ($line) {
	    ($score, $n_pos, $n_aa) = split(/\s+/, $line);
	}
    }
    if ($ref_score and $score and $n_pos and $n_aa) {
	print "$cluster_pair[$i]\t$ref_score\t$score\t", ($score - $ref_score)/($n_pos * $n_aa), "\n";
    }
}
