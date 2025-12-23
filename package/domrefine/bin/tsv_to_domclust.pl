#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM
";

use DomRefine::Read;

### Settings ###
my %OPT;
getopts('', \%OPT);

my $TMP_INPUT = define_tmp_file("$PROGRAM.input");
END {
    remove_tmp_file($TMP_INPUT);
}

### Main ###
-t and die $USAGE;
my $DCLST = save_stdin($TMP_INPUT);

my %domain = ();
my %cluster = ();
get_dclst_structure($TMP_INPUT, \%cluster, \%domain);

for my $cluster (sort {$a cmp $b} keys %cluster) {
# for my $cluster (sort {$a <=> $b} keys %cluster) {
    print "Cluster $cluster\n";
    for my $gene (sort {$a cmp $b} keys %{$cluster{$cluster}}) {
	for my $domain (sort {$a <=> $b} keys %{$domain{$gene}}) {
	    my $cluster_for_this_domain = $domain{$gene}{$domain}{cluster};
	    my $begin = $domain{$gene}{$domain}{begin};
	    my $end = $domain{$gene}{$domain}{end};
	    if ($cluster_for_this_domain) {
		if ($cluster_for_this_domain eq $cluster) {
		    print "$gene";
		    if ($domain) {
			print "($domain)";
		    }
		    print " $begin $end\n";
		}
	    } else {
		print STDERR "error: [$cluster] $gene ($domain)\n";
	    }
	}
    }
    print "\n";
}
