package DomRefine::COG;
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(read_cog_seq get_cog_domain print_cog_genes print_cog_domain print_cog_seq read_cog_tax);

use strict;

sub decompose_cog_member_id {
    my ($gene_domain) = @_;

    # if ($gene_domain =~ /^(\S+)_(\d)$/) {
    if ($gene_domain =~ /^(\S+)_(0?\d)$/) { # aqq_01 ... aqq_07
	my ($gene, $domain) = ($1, $2);
	return ($gene, $domain);
    } else {
	return ($gene_domain, 0);
    }
}

sub get_cog_domain {
    my ($r_gene_domain_seq, $r_gene_domain_pos) = @_;

    for my $gene (sort { $a cmp $b } keys %{$r_gene_domain_seq}) {
	my @domain = sort {$a <=> $b} keys %{${$r_gene_domain_seq}{$gene}};
	if (@domain == 1 and $domain[0] == 0) {
	    my $seq = ${$r_gene_domain_seq}{$gene}{0};
	    my $len = length($seq);
	    ${$r_gene_domain_pos}{$gene}{0}{start} = 1;
	    ${$r_gene_domain_pos}{$gene}{0}{end} = $len;
	} else {
	    my $cum_len = 0;
	    for (my $i=0; $i<@domain; $i++) {
		if ($domain[$i] != $i + 1) {
		    die;
		}
		my $domain_seq = ${$r_gene_domain_seq}{$gene}{$domain[$i]};
		my $domain_len = length($domain_seq);
		my $start = $cum_len + 1;
		$cum_len += $domain_len;
		my $end = $cum_len;
		${$r_gene_domain_pos}{$gene}{$domain[$i]}{start} = $start;
		${$r_gene_domain_pos}{$gene}{$domain[$i]}{end} = $end;
	    }
	}
    }
}

sub print_cog_seq {
    my ($r_gene_domain_seq, $r_gene) = @_;

    for my $gene (@{$r_gene}) {
	my @domain = sort {$a <=> $b} keys %{${$r_gene_domain_seq}{$gene}};
	if (@domain == 1 and $domain[0] == 0) {
	    my $seq = ${$r_gene_domain_seq}{$gene}{0};
	    print_a_gene_seq($gene, $seq);
	} else {
	    my $seq = "";
	    for (my $i=0; $i<@domain; $i++) {
		if ($domain[$i] != $i + 1) {
		    die;
		}
		my $domain_seq = ${$r_gene_domain_seq}{$gene}{$domain[$i]};
		$seq .= $domain_seq;
	    }
	    print_a_gene_seq($gene, $seq);
	}
    }
}

sub print_a_gene_seq {
    my ($gene, $seq) = @_;
    
    if ($seq eq "") {
	print STDERR "no seq for $gene\n";
    }
    print ">$gene\n";
    $seq =~ s/(.{1,60})/$1\n/g;
    print "$seq";
}

sub print_cog_domain {
    my ($r_gene_domain_pos) = @_;

    my @gene = sort { $a cmp $b } keys %{$r_gene_domain_pos};
    for my $gene (@gene) {
	my @domain = sort {$a <=> $b} keys %{${$r_gene_domain_pos}{$gene}};
	for my $domain (@domain) {
	    my $start = ${$r_gene_domain_pos}{$gene}{$domain}{start};
	    my $end = ${$r_gene_domain_pos}{$gene}{$domain}{end};
	    print "$gene $domain $start $end\n";
	}
    }
}

sub print_cog_genes {
    my ($cog, $org, $genes, $r_gene_domain_pos) = @_;

    my @genes = split(/\s+/, $genes);
    for my $gene_domain (@genes) {
	my ($gene, $domain) = decompose_cog_member_id($gene_domain);
	my $start = ${$r_gene_domain_pos}{$gene}{$domain}{start};
	my $end = ${$r_gene_domain_pos}{$gene}{$domain}{end};
	print "$cog $org$gene $domain $start $end\n";
    }
}

sub read_cog_seq {
    my ($file, $r_gene_domain_seq) = @_;

    open(FILE, $file) || die;
    my $gene;
    my $domain;
    while (<FILE>) {
	chomp;
	if (/^>(\S+)/) {
	    my $gene_domain = $1;
	    ($gene, $domain) = decompose_cog_member_id($gene_domain);
	    ${$r_gene_domain_seq}{$gene}{$domain} = "";
	} elsif (/^[A-Z]+$/) {
	    ${$r_gene_domain_seq}{$gene}{$domain} .= $_;
	} else {
	    die $_;
	}
    }
    close(FILE);
}

sub read_cog_tax {
    my ($file, $r_cog_tax) = @_;
    
    my $index = 1;
    open(FILE, $file) || die;
    while (<FILE>) {
	if (/^{(.*)}$/) {
	    my $species = $1;
	    my @species = split(",", $species);
	    for my $sp (@species) {
		${$r_cog_tax}{$sp} = $index;
	    }
	    $index ++;
	}
    }
    close(FILE);
}

1;
