#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: cat COG.txt | $PROGRAM COG_SEQ.fa
-a ANNOT_OUT_FILE
";

use DomRefine::COG;

my %OPT;
getopts('a:', \%OPT);

my %GENE_DOMAIN_SEQ;
if (@ARGV != 1) {
    print STDERR $USAGE;
    exit 1;
}
my ($COG_SEQ_FILE) = @ARGV;
read_cog_seq($COG_SEQ_FILE, \%GENE_DOMAIN_SEQ);

my %GENE_DOMAIN_POS;
get_cog_domain(\%GENE_DOMAIN_SEQ, \%GENE_DOMAIN_POS);

if ($OPT{a}) {
    my $annot_out_file = "$OPT{a}";
    open(ANNOT, ">$annot_out_file") || die;
}

!@ARGV && -t and die $USAGE;
my $COG;
my $ORG;
while (<STDIN>) {
    if (/^(\[[A-Z]+\]) (COG\d+) (.*)/) {
	my $cat = $1;
	$COG = $2;
	my $descr = $3;
	if ($OPT{a}) {
	    print ANNOT "$COG\t$cat $descr\n";
	}
    } elsif (/^  (\S+)  (.*)/) {
	$ORG = $1;
	my $genes = $2;
	print_cog_genes($COG, $ORG, $genes, \%GENE_DOMAIN_POS);
    } elsif (/^ +(.*)/) {
	my $genes = $1;
	print_cog_genes($COG, $ORG, $genes, \%GENE_DOMAIN_POS);
    } elsif (/^_/) {
    } elsif (/^$/) {
    } else {
	die $_;
    }
}

if ($OPT{a}) {
    close(ANNOT);
}
