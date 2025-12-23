#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: cat CLUSTER_SET | $PROGRAM -i CLUSTER_FILE -r REFERENCE_FILE
";

use DomRefine::General;
use DomRefine::Read;
use DomRefine::Motif;
use DomRefine::Refine;

### Settings ###
my %OPT;
getopts('i:r:', \%OPT);

if (! $OPT{i} and ! $OPT{r}) {
    die $USAGE;
}

### Get reference ###
my %REFERENCE_DOMAIN = ();
my %REFERENCE = ();
get_dclst_structure($OPT{r}, \%REFERENCE, \%REFERENCE_DOMAIN);

### Main ###
-t and die $USAGE;
while (my $line = <STDIN>) {
    chomp($line);
    my ($cluster_set) = split("\t", $line);
    my @cluster = split(/[-,\s]/, $cluster_set);
    my $dclst_part = extract_dclst($OPT{i}, @cluster);

    print STDERR "before:\n";
    my @before = dom_eval($dclst_part, \%REFERENCE_DOMAIN);
    print STDERR "after:\n";
    my @after = dom_eval(renumber_concat_domains(merge_all($dclst_part)), \%REFERENCE_DOMAIN);

    if ($before[0]) {
	my $change = ($after[1] - $before[1]) / $before[0];
	print "$line\t$change = ($after[1] - $before[1]) / $before[0]\n";
    }
}

################################################################################
### Functions ##################################################################
################################################################################

sub dom_eval {
    my ($dclst_part, $r_reference_domain) = @_;

    my %domain = ();
    my %cluster = ();
    parse_dclst_to_structure($dclst_part, \%cluster, \%domain);

    # hit of reference motifs for the selected domains
    my @hit_gene;
    my @hit_start;
    my @hit_end;
    my @hit_motif;
    my @gene = keys %domain;
    for (my $i=0; $i<@gene; $i++) {
	get_hit_pos($gene[$i], \@hit_start, \@hit_end, \@hit_motif, $i, $r_reference_domain, r_hit_gene => \@hit_gene);
    }
    # sum of overlap between the selected domains and reference motifs
    my %sum_overlap_len = ();
    for (my $k=0; $k<@hit_motif; $k++) {
	sum_overlap_len(\%domain, $hit_gene[$k], $hit_start[$k], $hit_end[$k], $hit_motif[$k], \%sum_overlap_len);
    }
    # identify pairs between domains and reference motifs
    my %corresp = ();
    get_correspondency(\%sum_overlap_len, \%corresp);

    # summarize overlap statistics
    my %overlap = ();
    for (my $k=0; $k<@hit_start; $k++) {
	$overlap{$hit_motif[$k]}{overlap} += get_overlap(\%domain, $hit_gene[$k], $hit_start[$k], $hit_end[$k], $corresp{motif}{$hit_motif[$k]});
	$overlap{$hit_motif[$k]}{count} ++;
    }

    my $count_total = 0;
    my $overlap_total = 0;
    for my $ref (keys %overlap) {
	print STDERR join("\t", $ref, $overlap{$ref}{count}, $overlap{$ref}{overlap}, $overlap{$ref}{overlap} / $overlap{$ref}{count}), "\n"; 
	$count_total += $overlap{$ref}{count};
	$overlap_total += $overlap{$ref}{overlap};
    }
    print STDERR join("\t", $count_total, $overlap_total), "\n";

    return($count_total, $overlap_total);
}

sub get_overlap {
    my ($r_domain, $gene, $hit_start, $hit_end, $cluster) = @_;

    if (! defined $cluster) {
	return 0;
    }

    for my $domain (keys %{${$r_domain}{$gene}}) {
	if (${$r_domain}{$gene}{$domain}{cluster} eq $cluster) {
	    my $overlap_len = overlap_len($hit_start, $hit_end, ${$r_domain}{$gene}{$domain}{begin}, ${$r_domain}{$gene}{$domain}{end});
	    my $overlap_ratio_motif = $overlap_len / ($hit_end - $hit_start + 1);
	    return  $overlap_ratio_motif;
	}
    }
    return 0;
}

sub sum_overlap_len {
    my ($r_domain, $gene, $hit_start, $hit_end, $hit_motif, $r_overlap_len) = @_;
    
    for my $domain (keys %{${$r_domain}{$gene}}) {
	my $cluster = ${$r_domain}{$gene}{$domain}{cluster};
	${$r_overlap_len}{$hit_motif}{$cluster} += overlap_len($hit_start, $hit_end, ${$r_domain}{$gene}{$domain}{begin}, ${$r_domain}{$gene}{$domain}{end});
    }
}

sub get_correspondency {
    my ($r_overlap_len, $r_corresp) = @_;

    my @motif = ();
    my @cluster = ();
    my @overlap = ();

    for my $motif (keys %{$r_overlap_len}) {
	for my $cluster (keys %{${$r_overlap_len}{$motif}}) {
	    push @motif, $motif;
	    push @cluster, $cluster;
	    push @overlap, ${$r_overlap_len}{$motif}{$cluster};
	}
    }

    my @sorted_i = sort {$overlap[$b] <=> $overlap[$a]} (0..$#motif);
    my %check = ();
    for my $i (@sorted_i) {
	if (! $check{motif}{$motif[$i]} and ! $check{cluster}{$cluster[$i]} and $overlap[$i]>0) {
	    ${$r_corresp}{cluster}{$cluster[$i]} = $motif[$i];
	    ${$r_corresp}{motif}{$motif[$i]} = $cluster[$i];
	    print STDERR "$motif[$i] <-> [$cluster[$i]] ($overlap[$i] aa overlap)\n";
	    $check{motif}{$motif[$i]} = 1;
	    $check{cluster}{$cluster[$i]} = 1;
	}
    }
}
