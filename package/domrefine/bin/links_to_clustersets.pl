#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: cat link | $PROGRAM
";

use DomRefine::General;
use DomRefine::Refine;

my %OPT;
getopts('mn', \%OPT);

!@ARGV && -t and die $USAGE;

my %GENE_TO_CLUSTER = ();
while (<>) {
    chomp;
    my @cluster = split(/[-,\s]/, $_);
    add_group_member(\%GENE_TO_CLUSTER, $., @cluster);
}

my %MODULE_TO_CLUSTER = ();
link_modules(\%GENE_TO_CLUSTER, \%MODULE_TO_CLUSTER);
print_module_members(\%MODULE_TO_CLUSTER, %OPT);
