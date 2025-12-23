package DomRefine::Motif;
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(read_motif_cutoff sort_hit_motifs get_cog_hits_from_file
	     get_hit_positions get_hit_pos read_motif_list get_hit_positions_sparql
	     read_cog_hits get_overlaps reverse_overlaps print_cluster_overlaps print_cluster_overlap_max count_overlap count_overlap_fast get_overlaps_check print_cluster_overlaps_1to1
             calc_overlap_aa print_overlap_with_homcluster
             check_merge_by_hom check_merge_by_hom_mysql check_hom_mysql find_hom_mysql
	     );

use strict;
use DomRefine::General;
use DomRefine::Read;

my $TMP_PRINT_CLUSTER_OVERLAP_MAX = define_tmp_file("print_cluster_overlap_max");
END {
    remove_tmp_file($TMP_PRINT_CLUSTER_OVERLAP_MAX);
}

sub read_motif_list {
    my ($motif_list_file, $r_motif_list) = @_;
    
    open(MOTIF_LIST, $motif_list_file) || die;
    while (<MOTIF_LIST>) {
	my ($motif) = split;
	${$r_motif_list}{$motif} = 1;
    }
    close(MOTIF_LIST);

}

sub read_motif_len {
    my ($motif_len_file, $r_motif_len) = @_;
    
    open(MOTIF, $motif_len_file) || die;
    while (<MOTIF>) {
	my @f = split("\t", $_);
	my ($motif, $len) = ($f[1], $f[3]);
	${$r_motif_len}{$motif} = $len;
    }
    close(MOTIF)
}

sub read_motif_cutoff {
    my ($motif_cutoff_file, $r_motif_cutoff) = @_;
    
    open(MOTIF, $motif_cutoff_file) || die;
    while (<MOTIF>) {
	my ($motif, $cutoff) = split;
	${$r_motif_cutoff}{$motif} = $cutoff;
    }
    close(MOTIF)
}

sub get_motif_len_mysql {
    my ($dbh, $motif) = @_;

    my $r_r_len = $dbh->selectall_arrayref("select length from motif where motid='$motif'");
    if (@{$r_r_len} != 1 or @{${$r_r_len}[0]} != 1) {
	die;
    }
    my $len = ${$r_r_len}[0][0];
    
    return $len;
}

sub get_cog_hits_from_file {
    my ($r_cogs, $r_motif_cutoff, $r_hit_gene, $r_hit_start, $r_hit_end, $r_hit_score, $file) = @_;

    open(FILE, $file) || die;
    while (<FILE>) {
	chomp;
	my @f = split("\t", $_);
	if (match($f[7], @{$r_cogs})) {
	    if ($f[12] > ${$r_motif_cutoff}{$f[7]}) {
		push @{$r_hit_gene}, "$f[2]:$f[3]";
		push @{$r_hit_start}, $f[4];
		push @{$r_hit_end}, $f[5];
		push @{$r_hit_score}, $f[12];
	    }
	}
    }
    close(FILE);
}

### Does not accept fused genes. Should be modified (c.f. get_hit_positions)
sub get_hit_positions_sparql {
    my ($r_gene, $r_get_j, $reference, $r_hit_motif, $r_hit_gene_i, $r_hit_start_j, $r_hit_end_j, $r_hit_evalue, %opt) = @_;

    my @motifs_to_ignore = ();
    if ($opt{motifs_to_ignore}) {
	@motifs_to_ignore = split(",", $opt{motifs_to_ignore});
    }

    my %get_i = ();
    for (my $i=0; $i<@{$r_gene}; $i++) {
	my $gene = ${$r_gene}[$i];
	$get_i{$gene} = $i;
    }

    my $option = "";
    if ($reference =~ /tigr/) {
	$option = "-t";
    } elsif ($reference =~ /pfam/) {
	$option = "-p";
    }

    my @line = `gene2motif $option @{$r_gene}`;
    chomp(@line);
    for my $line (@line) {
	my @f = split("\t", $line);
	my ($gene, $hit_start, $hit_end, $hit_motif, $evalue) = @f;
	my $i = $get_i{$gene};
	my $hit_start_j = ${${$r_get_j}[$i]}{$hit_start};
	my $hit_end_j = ${${$r_get_j}[$i]}{$hit_end};
	# print join("\t", $i, $hit_motif, $hit_start_j), "\n";
	if (grep { /^$hit_motif$/ } @motifs_to_ignore) {
	    next;
	}
	push @{$r_hit_motif}, $hit_motif;
	push @{$r_hit_gene_i}, $i;
	push @{$r_hit_start_j}, $hit_start_j;
	push @{$r_hit_end_j}, $hit_end_j;
	push @{$r_hit_evalue}, $evalue;
    }
}

### Will be obsolete. Use SPARQL version instead.
sub get_hit_positions {
    my ($fused_gene, $r_hit_start_j, $r_hit_end_j, $r_hit_motif, $r_hit_gene_i, $r_hit_evalue, $r_get_j_row, $i, $r_reference_domain, %opt) = @_;

    my @fused_gene = split(/\|/, $fused_gene);
    my $previous_gene = "";
    my $offset = 0;
    for my $gene (@fused_gene) {
	if ($previous_gene ne "" and $previous_gene ne $gene) {
	    $offset += get_gene_length($opt{r_seq}, $previous_gene);
	}
	get_hit_pos($gene, $r_hit_start_j, $r_hit_end_j, $r_hit_motif, $i, $r_reference_domain, %opt, offset => $offset, r_get_j_row => $r_get_j_row, r_hit_gene_i => $r_hit_gene_i);
	$previous_gene = $gene;
    }
}

sub read_cog_hits {
    my ($file, $r_cog_hits) = @_;
    
    print STDERR "Reading COG hits..\n";
    open(COG_HITS, "$file") || die;
    while (<COG_HITS>) {
	chomp;
	my @f = split("\t", $_);
	my ($sp, $name) = @f[2,3];
	my $gene = "$sp:$name";

	${$r_cog_hits}{line}{$.} = $_;
	if (${$r_cog_hits}{gene}{$gene}) {
	    push @{${$r_cog_hits}{gene}{$gene}}, $.;
	} else {
	    ${$r_cog_hits}{gene}{$gene} = [$.];
	}
    }
    close(COG_HITS);
}

sub get_hit_pos {
    my ($gene, $r_hit_start_j, $r_hit_end_j, $r_hit_motif, $i, $r_cog_domain, %opt) = @_;
    my $offset = $opt{offset};

    for my $domain (keys %{${$r_cog_domain}{$gene}}) {
	my $hit_start = ${$r_cog_domain}{$gene}{$domain}{begin};
	my $hit_end = ${$r_cog_domain}{$gene}{$domain}{end};
	my $motif = ${$r_cog_domain}{$gene}{$domain}{cluster};
    
	if ($offset) {
	    $hit_start += $offset;
	    $hit_end += $offset;
	}
	push @{$r_hit_motif}, $motif;
	if ($opt{r_get_j_row}) {
	    my $r_get_j_row = $opt{r_get_j_row};
	    push @{$r_hit_start_j}, ${$r_get_j_row}{$hit_start};
	    push @{$r_hit_end_j}, ${$r_get_j_row}{$hit_end};
	} else {
	    push @{$r_hit_start_j}, $hit_start;
	    push @{$r_hit_end_j}, $hit_end;
	}
	if ($opt{r_hit_gene_i}) {
	    my $r_hit_gene_i = $opt{r_hit_gene_i};
	    push @{$r_hit_gene_i}, $i;
	}
	if ($opt{r_hit_gene}) {
	    my $r_hit_gene = $opt{r_hit_gene};
	    push @{$r_hit_gene}, $gene;
	}
	
    }
}

sub sort_hit_motifs {
    my @motif = @_;

    my %count = ();
    for my $motif (@motif) {
	$count{$motif}++;
    }

    return sort {$count{$b} <=> $count{$a}} keys %count;
}

sub get_overlaps {
    my ($r_domain, $r_motif, $r_overlap, %opt) = @_;
    
    for my $gene (keys %{$r_motif}) {
	for my $motif_no (keys %{${$r_motif}{$gene}}) {
	    my $motif = ${$r_motif}{$gene}{$motif_no}{cluster};
	    my $motif_begin = ${$r_motif}{$gene}{$motif_no}{begin};
	    my $motif_end = ${$r_motif}{$gene}{$motif_no}{end};
	    my $motif_len = ($motif_end - $motif_begin + 1);
	    for my $domain_no (keys %{${$r_domain}{$gene}}) {
		my $cluster = ${$r_domain}{$gene}{$domain_no}{cluster};
		my $begin = ${$r_domain}{$gene}{$domain_no}{begin};
		my $end = ${$r_domain}{$gene}{$domain_no}{end};
		my $domain_len = ($end - $begin + 1);
		my $overlap_len = overlap_len($motif_begin, $motif_end, $begin, $end);
		if ($overlap_len > 0) {
		    my $r_over = $overlap_len / max($motif_len, $domain_len);
		    my $r_over1 = $overlap_len / $domain_len;
		    my $r_over2 = $overlap_len / $motif_len;
		    if (! ${$r_overlap}{$motif}{$cluster}{common}) {
			${$r_overlap}{$motif}{$cluster}{common} = 0;
		    }
		    if (defined $opt{r_over}) {
			unless ($r_over > $opt{r_over}) {
			    next;
			}
		    } elsif (defined $opt{r_over2}) {
			unless ($r_over2 > $opt{r_over2}) {
			    next;
			}
		    } elsif (defined $opt{r_over_t}) {
			unless ($r_over >= $opt{r_over_t}) {
			    next;
			}
		    } elsif (defined $opt{r_over_T}) {
			unless ($r_over2 >= $opt{r_over_T}) {
			    next;
			}
		    }
		    ${$r_overlap}{$motif}{$cluster}{common} ++;
		    ${$r_overlap}{$motif}{$cluster}{aa} += $overlap_len;
		}
	    }
	}
    }
}

sub count_overlap {
    my ($r_motif, $r_domain, $r_cluster, $r_count) = @_;

    for my $gene (keys %{$r_motif}) {
	for my $motif_no (keys %{${$r_motif}{$gene}}) {
	    my $motif = ${$r_motif}{$gene}{$motif_no}{cluster};
	    my $motif_begin = ${$r_motif}{$gene}{$motif_no}{begin};
	    my $motif_end = ${$r_motif}{$gene}{$motif_no}{end};
	    my $motif_len = ($motif_end - $motif_begin + 1);
	    for my $domain_no (keys %{${$r_domain}{$gene}}) {
		my $cluster = ${$r_domain}{$gene}{$domain_no}{cluster};
		my $begin = ${$r_domain}{$gene}{$domain_no}{begin};
		my $end = ${$r_domain}{$gene}{$domain_no}{end};
		my $domain_len = ($end - $begin + 1);
		my $overlap_len = overlap_len($motif_begin, $motif_end, $begin, $end);
		my $overlap_len_ratio = $overlap_len / max($motif_len, $domain_len);
		if ($overlap_len > 0) {
		    if (${$r_cluster}{$cluster}) {
			${$r_count}{$motif}{$gene}{$motif_no} ++;
		    }
		}
	    }
	}
    }
}

sub count_overlap_fast {
    my ($r_motif_domain, $r_dclst_domain, $r_cluster, $r_count, %opt) = @_;

    for my $gene (keys %{$r_dclst_domain}) {
	for my $dclst_no (keys %{${$r_dclst_domain}{$gene}}) {
	    my $dclst = ${$r_dclst_domain}{$gene}{$dclst_no}{cluster};
	    if ($opt{skip}) {
	    } else {
		if (! ${$r_cluster}{$dclst}) {
		    next;
		}
	    }
	    my $dclst_begin = ${$r_dclst_domain}{$gene}{$dclst_no}{begin};
	    my $dclst_end = ${$r_dclst_domain}{$gene}{$dclst_no}{end};
	    my $dclst_len = ($dclst_end - $dclst_begin + 1);
	    for my $motif_no (keys %{${$r_motif_domain}{$gene}}) {
		my $motif = ${$r_motif_domain}{$gene}{$motif_no}{cluster};
		my $motif_begin = ${$r_motif_domain}{$gene}{$motif_no}{begin};
		my $motif_end = ${$r_motif_domain}{$gene}{$motif_no}{end};
		my $motif_len = ($motif_end - $motif_begin + 1);
		my $overlap_len = overlap_len($motif_begin, $motif_end, $dclst_begin, $dclst_end);
		my $overlap_len_ratio = $overlap_len / max($motif_len, $dclst_len);
		if ($overlap_len > 0) {
		    ${$r_count}{$motif}{$gene}{$motif_no} ++;
		}
	    }
	}
    }
}

sub get_overlaps_check {
    my ($r_dclst_domain, $r_motif_domain, $cluster_pair) = @_;

    my %cluster = ();
    for my $cluster (split(/[-,\s]/, $cluster_pair)) {
	$cluster{$cluster} = 1;
    }

    my %count = ();
    count_overlap_fast($r_motif_domain, $r_dclst_domain, \%cluster, \%count);

    my $n_splitted_motif = 0;
    my $n_motif = 0;
    for my $motif (keys %count) {
	for my $gene (keys %{$count{$motif}}) {
	    for my $no (keys %{$count{$motif}{$gene}}) {
		my $count = $count{$motif}{$gene}{$no};
		if ($count >= 2) {
		    $n_splitted_motif ++;
		}
		$n_motif ++;
	    }
	}
    }

    return $n_motif, $n_splitted_motif, $n_motif ? $n_splitted_motif / $n_motif : "NA", 
}

sub reverse_overlaps {
    my ($r_hash, $r_hash2) = @_;
    
    for my $key (keys %{$r_hash}) {
	for my $key2 (keys %{${$r_hash}{$key}}) {
	    ${$r_hash2}{$key2}{$key}{aa} = ${$r_hash}{$key}{$key2}{aa};
	    ${$r_hash2}{$key2}{$key}{common} = ${$r_hash}{$key}{$key2}{common};
	}
    }
}

sub print_cluster_overlaps_1to1 {
    my ($r_ref_clus, $r_input_count, $r_reference_count, %opt) = @_;

    my %clus_ref = ();
    reverse_overlaps($r_ref_clus, \%clus_ref);
    for my $reference (sort {$a cmp $b} keys %{$r_ref_clus}) {
	my @cluster = keys %{${$r_ref_clus}{$reference}};
	my $n_corresp_clusters = @cluster;
	for my $cluster (@cluster) {
	    my @reference = keys %{$clus_ref{$cluster}};
	    my $n_corresp_references = @reference;
	    if ($opt{one2one}) {
		unless ($n_corresp_references == 1 and $n_corresp_clusters == 1) {
		    next;
		}
	    }
	    my $n_common_member = ${$r_ref_clus}{$reference}{$cluster}{common};
	    my $n_member = ${$r_input_count}{$cluster};
	    my $n_reference_member = ${$r_reference_count}{$reference};
	    my $r_com = $n_common_member / max($n_member, $n_reference_member);
	    my $r_com1 = $n_common_member / $n_member;
	    my $r_com2 = $n_common_member / $n_reference_member;
	    my $f;
	    if ($r_com1 + $r_com2) {
		$f = 2 * $r_com1 * $r_com2 / ($r_com1 + $r_com2);
	    } else {
		$f = -1;
	    }
	    my $threshold = 0.7;
	    my $type = "others";
	    if ($f >= $threshold) {
	    	$type = "equiv";
	    } elsif ($r_com1 >= $threshold) {
	    	$type = "super";
	    } elsif ($r_com2 >= $threshold) {
	    	$type = "sub";
	    }
	    # my $type = "";
	    # if ($r_com1 >= $threshold and $r_com2 < $threshold) {
	    # 	$type = "super";
	    # } elsif ($r_com1 < $threshold and $r_com2 >= $threshold) {
	    # 	$type .= "sub";
	    # }
	    print join("\t", $cluster, "n=$n_member", "pair=${n_corresp_references}",
	    	       $reference, "n=$n_reference_member", "pair=${n_corresp_clusters}",
	    	       "n_com=$n_common_member", "r_com=$r_com", "r_com1=$r_com1", "r_com2=$r_com2", 
		       "f=$f", "$type"), "\n";
	}
    }
}

sub print_cluster_overlaps {
    my ($r_ref_clus, $r_input_count, $r_reference_count, %opt) = @_;

    my %print_ref_clus = ();
    for my $reference (sort {$a cmp $b} keys %{$r_ref_clus}) {
	for my $cluster (keys %{${$r_ref_clus}{$reference}}) {
	    my $n_common_member = ${$r_ref_clus}{$reference}{$cluster}{common};
	    my $n_member = ${$r_input_count}{$cluster};
	    my $n_reference_member = ${$r_reference_count}{$reference};
	    my $r_com = $n_common_member / max($n_member, $n_reference_member);
	    my $r_com1 = $n_common_member / $n_member;
	    my $r_com2 = $n_common_member / $n_reference_member;
	    if (defined $opt{r_com}) {
		if ($r_com > $opt{r_com}) {
		} else {
		    next;
		}
	    }
	    if (defined $opt{r_com2}) {
		if ($r_com2 > $opt{r_com2}) {
		} else {
		    next;
		}
	    }
	    $print_ref_clus{$reference}{$cluster}{common} = $n_common_member;
	}
    }

    print_cluster_overlaps_1to1(\%print_ref_clus, $r_input_count, $r_reference_count, %opt);
}

sub calc_overlap_aa {
    my ($r_ref_clus, $r_input_count, $r_reference_count, %opt) = @_;

    my $overlap_cluster = 0;
    my $overlap_member = 0;
    my $overlap_aa = 0;
    for my $reference (sort {$a cmp $b} keys %{$r_ref_clus}) {
	for my $cluster (keys %{${$r_ref_clus}{$reference}}) {
	    my $n_common_member = ${$r_ref_clus}{$reference}{$cluster}{common};
	    my $n_member = ${$r_input_count}{$cluster};
	    my $n_reference_member = ${$r_reference_count}{$reference};
	    my $r_com = $n_common_member / max($n_member, $n_reference_member);
	    if (defined $opt{r_com}) {
		if ($r_com > $opt{r_com}) {
		} else {
		    next;
		}
	    }
	    $overlap_cluster ++;
	    $overlap_member += ${$r_ref_clus}{$reference}{$cluster}{common};
	    $overlap_aa += ${$r_ref_clus}{$reference}{$cluster}{aa};
	}
    }
    return ($overlap_cluster, $overlap_member, $overlap_aa);
}

sub print_overlap_with_homcluster {
    my ($r_ref_clus, $r_input_count, $r_reference_count) = @_;

    open(TMP_PRINT_CLUSTER_OVERLAP_MAX, ">$TMP_PRINT_CLUSTER_OVERLAP_MAX") || die;
    for my $reference (sort {$a cmp $b} keys %{$r_ref_clus}) {
	for my $cluster (keys %{${$r_ref_clus}{$reference}}) {
	    my $aa = ${$r_ref_clus}{$reference}{$cluster}{aa};
	    print TMP_PRINT_CLUSTER_OVERLAP_MAX "$reference\t$cluster\t$aa\n";
	}
    }
    close(TMP_PRINT_CLUSTER_OVERLAP_MAX);
    my @result = `cat $TMP_PRINT_CLUSTER_OVERLAP_MAX | sort -t '\t' -k3,3gr`;

    open(TMP_PRINT_CLUSTER_OVERLAP_MAX, ">$TMP_PRINT_CLUSTER_OVERLAP_MAX") || die;
    my %selected_cluster = ();
    for my $result (@result) {
	chomp($result);
	my ($reference, $cluster, $aa) = split("\t", $result);
	if (! $selected_cluster{$cluster}) {
	    print TMP_PRINT_CLUSTER_OVERLAP_MAX join("\t", $reference, $cluster, 
		       $aa,
		       # ${$r_ref_clus}{$reference}{$cluster}{common} / max(${$r_input_count}{$cluster}, ${$r_reference_count}{$reference}), 
		       # ${$r_ref_clus}{$reference}{$cluster}{common} / ${$r_reference_count}{$reference}
		       ${$r_ref_clus}{$reference}{$cluster}{common} / ${$r_input_count}{$cluster}
		), "\n";
	    $selected_cluster{$cluster} = 1;
	}
    }
    close(TMP_PRINT_CLUSTER_OVERLAP_MAX);

    system "cat $TMP_PRINT_CLUSTER_OVERLAP_MAX | sort -t '\t' -k1,1n -k2,2n";
}

sub print_cluster_overlap_max {
    my ($r_overlap, $r_cluster_count, $r_motif_count) = @_;

    # sort
    open(TMP_PRINT_CLUSTER_OVERLAP_MAX, ">$TMP_PRINT_CLUSTER_OVERLAP_MAX") || die;
    for my $motif (sort {$a cmp $b} keys %{$r_overlap}) {
	for my $cluster (keys %{${$r_overlap}{$motif}}) {
	    my $aa = ${$r_overlap}{$motif}{$cluster}{aa};
	    print TMP_PRINT_CLUSTER_OVERLAP_MAX "$motif\t$cluster\t$aa\n";
	}
    }
    close(TMP_PRINT_CLUSTER_OVERLAP_MAX);
    my @result = `cat $TMP_PRINT_CLUSTER_OVERLAP_MAX | sort -t '\t' -k3,3gr`;

    # select max
    my %selected_motif = ();
    my %selected_cluster = ();
    for my $result (@result) {
	chomp($result);
	my ($motif, $cluster, $aa) = split("\t", $result);
	if (! $selected_motif{$motif} and ! $selected_cluster{$cluster}) {
	    print join("\t", 
		       $motif, $cluster, $aa,
		       ${$r_overlap}{$motif}{$cluster}{common} / max(${$r_cluster_count}{$cluster}, ${$r_motif_count}{$motif}), 
		       ${$r_overlap}{$motif}{$cluster}{common} / ${$r_motif_count}{$motif}
		), "\n";
	    $selected_motif{$motif} = 1;
	    $selected_cluster{$cluster} = 1;
	}
    }
}

sub check_merge_by_hom {
    my ($r_member, $cluster1, $cluster2, $r_hom, %opt) = @_;

    my @gene1 = keys %{${$r_member}{$cluster1}};
    my @gene2 = keys %{${$r_member}{$cluster2}};

    my $count_all = 0;
    my $count_homology = 0;
    for my $gene1 (@gene1) {
        for my $gene2 (@gene2) {
            if (${$r_member}{$cluster1}{$gene1} && ${$r_member}{$cluster1}{$gene2} and 
                ${$r_member}{$cluster2}{$gene1} && ${$r_member}{$cluster2}{$gene2}) {
            } elsif ($gene1 eq $gene2) {
            } else {
                $count_all ++;
                if (overlap_domain_hom($gene1, $gene2, $r_hom, $r_member, $cluster1, $cluster2, %opt) ||
                    overlap_domain_hom($gene2, $gene1, $r_hom, $r_member, $cluster1, $cluster2, %opt)) {
                    $count_homology ++;
                }
            }
        }
    }

    print STDERR "$cluster1-$cluster2\t$count_homology/$count_all\t";
    print STDERR $count_homology/$count_all if $count_all;
    print STDERR "\t", scalar(@gene1), "\t", scalar(@gene2), "\n";

    if ($count_all) {
        if ($count_homology/$count_all >= $opt{R}) {
            return 1;
        } else {
            return 0;
        }
    } else {
        print STDERR "WARNING: no links were evaluated\n";
        return 1;
    }
}

### Incomplete
sub find_hom_mysql {
    my ($r_gene, %opt) = @_;

    my %seq = ();
    read_seq(\%seq, $r_gene);

    my %chksum2gene = ();
    get_multi_seq_checksum_mysql($opt{dbh}, \%chksum2gene, @{$r_gene});
    my @chksum = keys %chksum2gene;

    my @condition1 = ();
    for my $chksum (@chksum) {
	push @condition1, "name1='$chksum'";
    }
    my $condition1 = join(" or ", @condition1);

    my @condition2 = ();
    for my $chksum (@chksum) {
	push @condition2, "name2='$chksum'";
    }
    my $condition2 = join(" or ", @condition2);

    my $r_r_homol = $opt{dbh_accum}->selectall_arrayref("SELECT from1, to1, from2, to2, eval, score, name1, name2 FROM homology WHERE ($condition1) or ($condition2)");
    for (my $i=0; $i<@{$r_r_homol}; $i++) {
	if (@{${$r_r_homol}[$i]} == 8) {
	    my ($start1, $end1, $start2, $end2, $e_value, $score, $name1, $name2) = @{${$r_r_homol}[$i]};
	    # print join("\t", ($start1, $end1, $start2, $end2, $e_value, $score, $name1, $name2)), "\n";
	    my $gene1 = get_genes_from_checksum($opt{dbh}, $name1);
	    my $gene2 = get_genes_from_checksum($opt{dbh}, $name2);
	    # my $gene1_len = get_gene_length(\%seq, $gene1);
	    # my $gene2_len = get_gene_length(\%seq, $gene2);
	    # print join("\t", $gene1, "$start1-$end1/$gene1_len", $gene2, "$start2-$end2/$gene2_len", $e_value, $score), "\n";

	    # my @gene1 = @{$chksum2gene{$name1}};
	    # my @gene2 = @{$chksum2gene{$name2}};
	    # if (@gene1 != 1) {
	    # 	print STDERR "WARNING: $name1 corresponds to @gene1\n";
	    # }
	    # if (@gene2 != 1) {
	    # 	print STDERR "WARNING: $name2 corresponds to @gene2\n";
	    # }
	    # for my $gene1 (@gene1) {
	    # 	for my $gene2 (@gene2) {
	    # 	    if (! $opt{all} and $gene1 eq $gene2) {
	    # 		next;
	    # 	    }
	    # 	    my $gene1_len = get_gene_length(\%seq, $gene1);
	    # 	    my $gene2_len = get_gene_length(\%seq, $gene2);
	    # 	    print join("\t", $gene1, "$start1-$end1/$gene1_len", $gene2, "$start2-$end2/$gene2_len", $e_value, $score), "\n";
	    # 	}
	    # }
	}
    }
}

sub check_hom_mysql {
    my ($r_gene, %opt) = @_;

    my %seq = ();
    read_seq(\%seq, $r_gene);

    my %chksum2gene = ();
    get_multi_seq_checksum_mysql($opt{dbh}, \%chksum2gene, @{$r_gene});
    my @chksum = keys %chksum2gene;

    my @condition1 = ();
    for my $chksum (@chksum) {
	push @condition1, "name1='$chksum'";
    }
    my $condition1 = join(" or ", @condition1);

    my @condition2 = ();
    for my $chksum (@chksum) {
	push @condition2, "name2='$chksum'";
    }
    my $condition2 = join(" or ", @condition2);

    my $r_r_homol = $opt{dbh_accum}->selectall_arrayref("SELECT from1, to1, from2, to2, eval, score, name1, name2 FROM homology WHERE ($condition1) and ($condition2)");
    for (my $i=0; $i<@{$r_r_homol}; $i++) {
	if (@{${$r_r_homol}[$i]} == 8) {
	    my ($start1, $end1, $start2, $end2, $e_value, $score, $name1, $name2) = @{${$r_r_homol}[$i]};
	    my @gene1 = @{$chksum2gene{$name1}};
	    my @gene2 = @{$chksum2gene{$name2}};
	    if (@gene1 != 1) {
		print STDERR "WARNING: $name1 corresponds to @gene1\n";
	    }
	    if (@gene2 != 1) {
		print STDERR "WARNING: $name2 corresponds to @gene2\n";
	    }
	    for my $gene1 (@gene1) {
		for my $gene2 (@gene2) {
		    if (! $opt{all} and $gene1 eq $gene2) {
			next;
		    }
		    my $gene1_len = get_gene_length(\%seq, $gene1);
		    my $gene2_len = get_gene_length(\%seq, $gene2);
		    print join("\t", $gene1, "$start1-$end1/$gene1_len", $gene2, "$start2-$end2/$gene2_len", $e_value, $score), "\n";
		}
	    }
	}
    }
}

sub check_merge_by_hom_mysql {
    my ($r_member, $cluster1, $cluster2, %opt) = @_;

    my @gene1 = keys %{${$r_member}{$cluster1}};
    my @gene2 = keys %{${$r_member}{$cluster2}};

    my %chksum2gene = ();
    get_multi_seq_checksum_mysql($opt{dbh}, \%chksum2gene, @gene1);
    get_multi_seq_checksum_mysql($opt{dbh}, \%chksum2gene, @gene2);
    # my @chksum = map { "'$_'" } keys %chksum2gene;
    my @chksum = keys %chksum2gene;
    my @condition1 = ();
    for my $chksum (@chksum) {
	push @condition1, "name1='$chksum'";
    }
    my @condition2 = ();
    for my $chksum (@chksum) {
	push @condition2, "name2='$chksum'";
    }

    my %homology = ();
    # my $r_r_homol = $opt{dbh_accum}->selectall_arrayref("SELECT from1, to1, from2, to2, eval, score, name1, name2 FROM homology WHERE name1 IN(" . join(",", @chksum) . ") and name2 IN(" . join(",", @chksum) . ")");
    my $r_r_homol = $opt{dbh_accum}->selectall_arrayref("SELECT from1, to1, from2, to2, eval, score, name1, name2 FROM homology WHERE (" . join(" or ", @condition1) . ") and (" . join(" or ", @condition2) . ")");
    for (my $i=0; $i<@{$r_r_homol}; $i++) {
	if (@{${$r_r_homol}[$i]} == 8) {
	    my ($start1, $end1, $start2, $end2, $e_value, $score, $name1, $name2) = @{${$r_r_homol}[$i]};
	    for my $gene1 (@{$chksum2gene{$name1}}) {
		for my $gene2 (@{$chksum2gene{$name2}}) {
		    # print STDERR join(" ", $name1, $name2, $gene1, $gene2, $start1, $end1, $start2, $end2, $e_value, $score), "\n";
		    if ($e_value < 0.001 and $score >= 60) {
			if (defined $opt{r}) {
			    my $overlap_gene1 = check_overlap_with_clusters($start1, $end1, $r_member, $cluster1, $cluster2, $gene1, %opt);
			    my $overlap_gene2 = check_overlap_with_clusters($start2, $end2 ,$r_member, $cluster1, $cluster2, $gene2, %opt);
			    if ($overlap_gene1 && $overlap_gene2) {
				$homology{$gene1}{$gene2} = 1;
				$homology{$gene2}{$gene1} = 1;
			    }
			} else {
			    $homology{$gene1}{$gene2} = 1;
			    $homology{$gene2}{$gene1} = 1;
			}
		    }
		}
	    }
	}
    }

    my $count_all = 0;
    my $count_homology = 0;
    for my $gene1 (@gene1) {
	for my $gene2 (@gene2) {
	    if (${$r_member}{$cluster1}{$gene1} && ${$r_member}{$cluster1}{$gene2} and 
		${$r_member}{$cluster2}{$gene1} && ${$r_member}{$cluster2}{$gene2}) {
	    } elsif ($gene1 eq $gene2) {
	    } else {
		$count_all ++;
		# if (overlap_domain_hom_mysql($gene1, $gene2, $r_member, $cluster1, $cluster2, %opt)) {
		if ($homology{$gene1}{$gene2}) {
		    $count_homology ++;
		}
	    }
	}
    }

    print STDERR "$cluster1-$cluster2\t$count_homology/$count_all";
    print STDERR "\t", $count_homology/$count_all if $count_all;
    print STDERR "\t", scalar(@gene1), "\t", scalar(@gene2), "\n";

    if ($count_all) {
	if ($count_homology/$count_all >= $opt{R}) {
	    return 1;
	} else {
	    return 0;
	}
    } else {
	print STDERR "WARNING: no links were evaluated\n";
	return 1;
    }
}

################################################################################
### Private functions ##########################################################
################################################################################

sub overlap_domain_hom {
    my ($gene1, $gene2, $r_hom, $r_member, $cluster1, $cluster2, %opt) = @_;

    # if (defined ${$r_hom}{$gene1}{$gene2}) {
    if (defined ${$r_hom}{"$gene1 $gene2"}) {
	if (defined $opt{r}) {
	    # my ($start1, $end1) = (${$r_hom}{$gene1}{$gene2}{start1}, ${$r_hom}{$gene1}{$gene2}{end1});
	    # my ($start2, $end2) = (${$r_hom}{$gene1}{$gene2}{start2}, ${$r_hom}{$gene1}{$gene2}{end2});
        # my $start_end = ${$r_hom}{$gene1}{$gene2};
        my $start_end = ${$r_hom}{"$gene1 $gene2"};
        my ($start1, $end1, $start2, $end2);
        if ($start_end =~ /^s1:(\S+) e1:(\S+) s2:(\S+) e2:(\S+)$/) {
            ($start1, $end1, $start2, $end2) = ($1, $2, $3, $4);
        }
	    my $overlap_gene1 = check_overlap_with_clusters($start1, $end1, $r_member, $cluster1, $cluster2, $gene1, %opt);
	    my $overlap_gene2 = check_overlap_with_clusters($start2, $end2, $r_member, $cluster1, $cluster2, $gene2, %opt);
	    # print STDERR "$gene1\t$gene2\t$start1\t$end1\t$start2\t$end2\n";
	    if ($overlap_gene1 && $overlap_gene2) {
		return 1;
	    }
	} else {
	    return 1;
	}
    }
    return 0;
}

sub overlap_domain_hom_mysql {
    my ($gene1, $gene2, $r_member, $cluster1, $cluster2, %opt) = @_;

    my $chksum1 = get_seq_checksum_mysql($opt{dbh}, $gene1);
    my $chksum2 = get_seq_checksum_mysql($opt{dbh}, $gene2);

    my @homology = get_homology_mysql($opt{dbh_accum}, $chksum1, $chksum2);
    my @homology_rev = get_homology_mysql($opt{dbh_accum}, $chksum2, $chksum1);
    if (@homology == 4 && @homology_rev == 4) {
    	die;
    }

    my ($start1, $end1, $start2, $end2);
    if (@homology == 4) {
	($start1, $end1, $start2, $end2) = @homology;
    } elsif (@homology_rev == 4) {
	($start2, $end2, $start1, $end1) = @homology_rev;
    }

    if (@homology == 4 or @homology_rev == 4) {
	if (defined $opt{r}) {
	    my $overlap_gene1 = check_overlap_with_clusters($start1, $end1, $r_member, $cluster1, $cluster2, $gene1, %opt);
	    my $overlap_gene2 = check_overlap_with_clusters($start2, $end2 ,$r_member, $cluster1, $cluster2, $gene2, %opt);
	    if ($overlap_gene1 && $overlap_gene2) {
		return 1;
	    }
	} else {
	    return 1;
	}
    }
    return 0;
}

sub check_overlap_with_clusters {
    my ($hom_start, $hom_end, $r_member, $cluster1, $cluster2, $gene, %opt) = @_;

    my ($overlap_len1, $domain_len1) = check_overlap($hom_start, $hom_end, $r_member, $cluster1, $gene);
    my ($overlap_len2, $domain_len2) = check_overlap($hom_start, $hom_end, $r_member, $cluster2, $gene);
	
    my $overlap_ratio = ($overlap_len1 + $overlap_len2) / ($domain_len1 + $domain_len2);
    
    if ($overlap_ratio > $opt{r}) {
	return 1;
    } else {
	return 0;
    }
}

sub check_overlap {
    my ($hom_start, $hom_end, $r_member, $cluster, $gene) = @_;

    my $overlap_len = 0;
    my $domain_len = 0;

    if (defined ${$r_member}{$cluster}{$gene}) {
	for my $domain (keys %{${$r_member}{$cluster}{$gene}}) {
	    my ($domain_start, $domain_end) = (${$r_member}{$cluster}{$gene}{$domain}{start} , ${$r_member}{$cluster}{$gene}{$domain}{end});
	    if ($domain_start > $domain_end) {
		die "$cluster $gene $domain $domain_start $domain_end";
	    }
	    $overlap_len += overlap_len($hom_start, $hom_end, $domain_start, $domain_end);
	    $domain_len += $domain_end - $domain_start + 1;
	}
    }

    return ($overlap_len, $domain_len);
}

### Incomplete
sub get_genes_from_checksum {
    my ($dbh, $chksum) = @_;

    my $r_r_id = $dbh->selectall_arrayref("select id from proteinseq where chksum='$chksum'");
    # my @id = ();
    # for (my $i=0; $i<@{$r_r_id}; $i++) {
    # 	if (@{${$r_r_id}[$i]} != 1) {
    # 	    die $i;
    # 	}
    # 	my ($id) = @{${$r_r_id}[$i]};
    # 	push @id, $id;
    # }
    # print "@id\n";

    if (@{$r_r_id} == 0) {
	return;
    }
    if (@{$r_r_id} != 1 or @{${$r_r_id}[0]} != 1) {
    	die;
    }
    my ($id) = @{${$r_r_id}[0]};

    # my $r_r_gene = $dbh->selectall_arrayref("select sp, name from gene where aaseq=$id");
    # if (@{$r_r_gene} != 1 or @{${$r_r_gene}[0]} != 2) {
    # 	die $id;
    # }
    # my ($sp, $name) = @{${$r_r_gene}[0]};

    # return "$sp:$name";
}

sub get_multi_seq_checksum_mysql {
    my ($dbh, $r_chksum2gene, @gene) = @_;

    my @condition = ();
    for my $gene (@gene) {
	my ($sp, $name) = decompose_gene_id($gene);
	push @condition, "sp='$sp' and name='$name'";
    }

    my $r_r_id = $dbh->selectall_arrayref("SELECT sp, name, aaseq FROM gene WHERE " . join(" or ", @condition));
    if (@gene != @{$r_r_id}) {
	print STDERR "WARNING: number of genes ", scalar(@gene), " != ", scalar(@{$r_r_id}), "\n";
    }

    my @id = ();
    my %id2gene = ();
    for (my $i=0; $i<@{$r_r_id}; $i++) {
	if (@{${$r_r_id}[$i]} != 3) {
	    die $i;
	}
	my ($sp, $name, $aaseq) = @{${$r_r_id}[$i]};
	if (defined $id2gene{$aaseq}) {
	    push @{$id2gene{$aaseq}}, "$sp:$name"
	} else {
	    $id2gene{$aaseq} = ["$sp:$name"];
	}
	push @id, $aaseq;
    }

    my $r_r_seq = $dbh->selectall_arrayref("SELECT id, chksum FROM proteinseq WHERE id IN(" . join(",", @id) . ")");
    # if (@gene != @{$r_r_seq}) {
    # 	print STDERR "WARNING: number of seq ", scalar(@gene), " != ", scalar(@{$r_r_seq}), "\n";
    # }

    my %id2chksum = ();
    for (my $i=0; $i<@{$r_r_seq}; $i++) {
	if (@{${$r_r_seq}[$i]} != 2) {
	    die $i;
	}
	my ($id, $chksum) = @{${$r_r_seq}[$i]};
	$id2chksum{$id} = $chksum;
	@{${$r_chksum2gene}{$chksum}} = @{$id2gene{$id}};
    }
}

sub get_seq_checksum_mysql {
    my ($dbh, $gene) = @_;

    my ($sp, $name) = decompose_gene_id($gene);

    my $r_r_id = $dbh->selectall_arrayref("select aaseq from gene where sp='$sp' and name='$name'");
    if (@{$r_r_id} != 1 or @{${$r_r_id}[0]} != 1) {
	die $gene;
    }
    my ($id) = @{${$r_r_id}[0]};

    my $r_r_seq = $dbh->selectall_arrayref("select chksum from proteinseq where id=$id");
    if (@{$r_r_seq} != 1 or @{${$r_r_seq}[0]} != 1) {
	die;
    }
    my ($chksum) = @{${$r_r_seq}[0]};

    return $chksum;
}

sub get_homology_mysql {
    my ($dbh, $chksum1, $chksum2) = @_;

    my $r_r_homol = $dbh->selectall_arrayref("select from1, to1, from2, to2, eval, score from homology where name1='$chksum1' and name2='$chksum2'");

    if (@{$r_r_homol} == 1) {
	if (@{${$r_r_homol}[0]} == 6) {
	    my ($start1, $end1, $start2, $end2, $e_value, $score) = @{${$r_r_homol}[0]};
	    if ($e_value < 0.001 and $score >= 60) {
		return ($start1, $end1, $start2, $end2);
	    }
	} else {
	    die;
	}
    }
}

1;
