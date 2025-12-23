#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM COG_SEQ.fa
-f GENE_LIST
-d : print domains
";

use DomRefine::COG;

my %OPT;
getopts('f:d', \%OPT);

my %GENE_DOMAIN_SEQ;
if (@ARGV != 1) {
    print STDERR $USAGE;
    exit 1;
}
my ($COG_SEQ_FILE) = @ARGV;
read_cog_seq($COG_SEQ_FILE, \%GENE_DOMAIN_SEQ);

my @GENE;
if ($OPT{f}) {
    @GENE = `cat $OPT{f}`;
    chomp(@GENE);
} else {
    @GENE  = sort { $a cmp $b } keys %GENE_DOMAIN_SEQ;
}

if ($OPT{d}) {
    my %GENE_DOMAIN_POS;
    get_cog_domain(\%GENE_DOMAIN_SEQ, \%GENE_DOMAIN_POS);
    print_cog_domain(\%GENE_DOMAIN_POS);
} else {
    print_cog_seq(\%GENE_DOMAIN_SEQ, \@GENE);
}
