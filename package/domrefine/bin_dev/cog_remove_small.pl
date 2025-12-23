#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: cat com.dom | $PROGRAM -t cog.tax
";

my %OPT;
getopts('t:', \%OPT);

use DomRefine::Read;
use DomRefine::General;
use DomRefine::COG;

my $TMP_INPUT = define_tmp_file("$PROGRAM.input");
END {
    remove_tmp_file($TMP_INPUT);
}

my $COG_TAX = $OPT{t} || die $USAGE;
my %COG_TAX;
read_cog_tax($COG_TAX, \%COG_TAX);

### Main ###
-t and die $USAGE;
save_stdin($TMP_INPUT);

my %cluster = ();
my %domain = ();
get_dclst_structure($TMP_INPUT, \%cluster, \%domain);

my @cluster = sort {$a cmp $b} keys %cluster;
for my $cluster (@cluster) {
    my @gene = sort {$a cmp $b} keys %{$cluster{$cluster}};
    my @tax = ();
    for my $gene (@gene) {
	my ($sp, $name) = decompose_gene_id($gene);
	my $tax = $COG_TAX{$sp};
	if ($tax) {
	    push @tax, $tax;
	} else {
	    push @tax, $sp;
	}
    }
    @tax =  uniq(@tax);
    if (@tax >= 3) {
	for my $gene (sort {$a cmp $b} keys %{$cluster{$cluster}}) {
	    for my $domain (sort {$a <=> $b} keys %{$domain{$gene}}) {
		my $cluster_for_this_domain = $domain{$gene}{$domain}{cluster};
		my $begin = $domain{$gene}{$domain}{begin};
		my $end = $domain{$gene}{$domain}{end};
		if ($cluster_for_this_domain) {
		    if ($cluster_for_this_domain eq $cluster) {
			print "$cluster $gene $domain $begin $end\n";
		    }
		} else {
		    print STDERR "error: [$cluster] $gene ($domain)\n";
		}
	    }
	}
    }
}
