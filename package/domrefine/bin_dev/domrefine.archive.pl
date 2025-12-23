#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: cat CLUSTER_LIST | $PROGRAM -[aN] -s CLUSTER_SIZE_FILE -n SIZE DIR_ORIGINAL DIR_NEW

Copy output for each cluster from directory to directory
 (only for clusters with cluster_size >= SIZE)
";

my %OPT;
getopts('aNs:n:', \%OPT);

### read cluster list ###
-t and die $USAGE;
my @CLUSTER_LIST = <STDIN>;
chomp(@CLUSTER_LIST);

### set file extensions ###
my $EXT = "out";

my $EXT_NEW;
if ($OPT{a}) {
    $EXT_NEW = "fa";
} elsif ($OPT{N}) {
    $EXT_NEW = "newick";
} else {
    die $USAGE;
}

### read cluster size ###
if (! $OPT{s} || ! $OPT{n}) {
    die $USAGE;
}
my %CLUSTER_SIZE = ();
open(SIZE, $OPT{s}) || die;
while (<SIZE>) {
    my ($cluster, $size) = split;
    $CLUSTER_SIZE{$cluster} = $size;
}
close(SIZE);

### set directories ###
if (@ARGV != 2) {
    print STDERR $USAGE;
    exit 1;
}
my ($DIR, $DIR_NEW) = @ARGV;

if (! -d $DIR) {
    die $USAGE;
}

if (! -e $DIR_NEW) {
    system "mkdir $DIR_NEW";
}

### copy files ###
for my $cluster (@CLUSTER_LIST) {
    if ($CLUSTER_SIZE{$cluster}) {
        if ($CLUSTER_SIZE{$cluster} >= $OPT{n}) {
            if (-s "$DIR/$cluster.$EXT") {
                if (! -s "$DIR_NEW/$cluster.$EXT_NEW") {
                    system "cp -a $DIR/$cluster.$EXT $DIR_NEW/$cluster.$EXT_NEW";
                }
            } else {
                print STDERR "WARNING: no output for cluster $cluster\n";
            }
        }
    } else {
        print STDERR "WARNING: size of cluster $cluster is unknown\n";
    }
}
