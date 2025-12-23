package DomRefine::Refine;
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(merge_all partial_move
             get_adjacency_information_detail get_adjacency_information get_adjacent_organism
             renumber_domain renumber_cluster renumber_concat_domains assign_new_domain_number
             get_patterns
             get_patterns_detail
             print_in_simple_text
             print_in_text
             extract_most_adjacent_cluster_pair
             print_module_members link_modules
             update_domain calculate_sp_score_of_clusters
             update_domain_mapping update_domain_mapping_local boundary_move boundary_move_gene
             get_terminal_j
             get_clusters_subset get_new_cluster_ids replace_obsolete_clusters update_clusters choose_new_cluster_id
             );

use strict;
use DomRefine::General;
use DomRefine::Read;
use DomRefine::Score;

my $TMP_MOVE_DCLST = define_tmp_file("move_dclst");
my $TMP_MOVE2_DCLST = define_tmp_file("move2_dclst");
my $TMP_REFINE_PL_DCLST = define_tmp_file("refine");
END {
    remove_tmp_file($TMP_MOVE_DCLST);
    remove_tmp_file($TMP_MOVE2_DCLST);
    remove_tmp_file($TMP_REFINE_PL_DCLST);
}

################################################################################
### Function ###################################################################
################################################################################

### Merge ###
sub merge_all {
    my ($dclst, $rep_cluster) = @_;

    if (! defined $rep_cluster) {
        $rep_cluster = get_rep_cluster($dclst);
    }

    $dclst =~ s/^(\S+) /$rep_cluster /gm;

    return $dclst;
}

### Re-number ###

sub renumber_domain {
    my ($dclst_table_file) = @_;
    
    # read
    my %domain = ();
    open(DCLST, $dclst_table_file) || die;
    while (my $line = <DCLST>) {
        chomp($line);
        my ($cluster, $gene, @domain_info) = decompose_dclst_line($line);
        unless (@domain_info and @domain_info % 3 == 0) {
            die;
        }
        for (my $i=0; $i<@domain_info; $i+=3) {
            my ($domain_no, $begin_pos, $end_pos) = ($domain_info[$i], $domain_info[$i+1], $domain_info[$i+2]);
            my $domain_id = "${begin_pos}-${end_pos}:$cluster";
            if (defined $domain{$gene}{$domain_id}) {
                print STDERR "WARNING: duplicated domains exist. $gene $domain_id\n";
            }
            $domain{$gene}{$domain_id}{cluster} = $cluster;
            $domain{$gene}{$domain_id}{begin} = $begin_pos;
            $domain{$gene}{$domain_id}{end} = $end_pos;
        }
    }
    close(DCLST);

    # renumber
    my $out = "";
    for my $gene (keys %domain) {
        my @domain_id = sort { $domain{$gene}{$a}{begin} <=> $domain{$gene}{$b}{begin} or 
                                   $domain{$gene}{$a}{end} <=> $domain{$gene}{$b}{end} or
                                   $domain{$gene}{$a}{cluster} cmp $domain{$gene}{$b}{cluster} } keys %{$domain{$gene}};
        if (@domain_id == 1) {
            my $cluster = $domain{$gene}{$domain_id[0]}{cluster};
            my $begin_pos = $domain{$gene}{$domain_id[0]}{begin};
            my $end_pos = $domain{$gene}{$domain_id[0]}{end};
            $out .= "$cluster $gene 0 $begin_pos $end_pos\n";
        } else {
            my $i = 1;
            for my $domain_id (@domain_id) {
                my $cluster = $domain{$gene}{$domain_id}{cluster};
                my $begin_pos = $domain{$gene}{$domain_id}{begin};
                my $end_pos = $domain{$gene}{$domain_id}{end};
                $out .= "$cluster $gene $i $begin_pos $end_pos\n";
                $i ++;
            }
        }
    }

    return $out;
}

sub renumber_cluster {
    my ($dclst) = @_;

    # parse
    my %member = ();
    read_cluster_members($dclst, \%member);
    
    my %count = ();
    for my $cluster (keys %member) {
        for my $gene (sort {$a cmp $b} keys %{$member{$cluster}}) {
            for my $domain (sort {$a <=> $b} keys %{$member{$cluster}{$gene}}) {
                $count{$cluster} ++;
            }
        }
    }

    # output
    my $out = "";
    my $cluster_id = 1;
    for my $cluster (sort {$count{$b} <=> $count{$a} or $a cmp $b} keys %member) {
        for my $gene (sort {$a cmp $b} keys %{$member{$cluster}}) {
            for my $domain (sort {$a <=> $b} keys %{$member{$cluster}{$gene}}) {
                my $start = $member{$cluster}{$gene}{$domain}{start};
                my $end = $member{$cluster}{$gene}{$domain}{end};
                $out .= "$cluster_id $gene $domain $start $end\n";
            }
        }
        $cluster_id ++;
    }

    return $out;
}

sub renumber_concat_domains {
    my ($dclst) = @_;

    # parse
    my %member = ();
    read_cluster_members($dclst, \%member);

    # re-number
    concat_members_in_cluster(\%member);
    my %new_domain_no = ();
    assign_new_domain_number(\%member, \%new_domain_no);

    # output
    my $out = "";
    for my $cluster (sort {$a cmp $b} keys %member) {
        for my $gene (sort {$a cmp $b} keys %{$member{$cluster}}) {
            for my $domain (sort {$a <=> $b} keys %{$member{$cluster}{$gene}}) {
                my $start = $member{$cluster}{$gene}{$domain}{start};
                my $end = $member{$cluster}{$gene}{$domain}{end};
                my $new_domain_no = $new_domain_no{$gene}{$domain};
                $out .= "$cluster $gene $new_domain_no $start $end\n";
            }
        }
    }

    return $out;
}

sub assign_new_domain_number {
    my ($r_member, $r_new_domain_no) = @_;

    for my $cluster (keys %{$r_member}) {
        for my $gene (keys %{${$r_member}{$cluster}}) {
            for my $old_domain_no (keys %{${$r_member}{$cluster}{$gene}}) {
                ${$r_new_domain_no}{$gene}{$old_domain_no} = $old_domain_no; # initialize by old domain number
            }
        }
    }

    for my $gene (keys %{$r_new_domain_no}) {
        my @old_domain_no = sort {$a<=>$b} keys %{${$r_new_domain_no}{$gene}};
        if (@old_domain_no == 1) {
            ${$r_new_domain_no}{$gene}{$old_domain_no[0]} = 0; # only one domain
        } else {
            my $count = 0;
            for my $old_domain_no (@old_domain_no) {
                $count ++;
                my $new_domain_no = $count;
                ${$r_new_domain_no}{$gene}{$old_domain_no} = $new_domain_no; # new domain number
            }
        }
    }
}

# sub renumber_concat_domains {
#     my ($dclst) = @_;

#     # parse
#     my %member = ();
#     read_cluster_members($dclst, \%member);

#     # re-number
#     concat_members_in_cluster(\%member);
#     my %new_domain_no = ();
#     assign_new_domain_number(\%member, \%new_domain_no);

#     # output
#     my $out = "";
#     for my $cluster (sort {$a <=> $b} keys %member) {
# 	for my $gene (sort {$a cmp $b} keys %{$member{$cluster}}) {
# 	    for my $domain (sort {$a <=> $b} keys %{$member{$cluster}{$gene}}) {
# 		my $start = $member{$cluster}{$gene}{$domain}{start};
# 		my $end = $member{$cluster}{$gene}{$domain}{end};
# 		$out .= "$cluster $gene $domain $start $end\n";
# 	    }
# 	}
#     }

#     return $out;
# }

# sub assign_new_domain_number {
#     my ($r_member, $r_new_domain_no) = @_;

#     for my $cluster (keys %{$r_member}) {
# 	for my $gene (keys %{${$r_member}{$cluster}}) {
# 	    for my $old_domain_no (keys %{${$r_member}{$cluster}{$gene}}) {
# 		${$r_new_domain_no}{$gene}{$old_domain_no} = $old_domain_no; # initialize by old domain number
# 	    }
# 	}
#     }

#     for my $gene (keys %{$r_new_domain_no}) {
# 	my @old_domain_no = sort {$a<=>$b} keys %{${$r_new_domain_no}{$gene}};
# 	if (@old_domain_no == 1) {
# 	    ${$r_new_domain_no}{$gene}{$old_domain_no[0]} = 0; # only one domain
# 	} else {
# 	    my $count = 0;
# 	    for my $old_domain_no (@old_domain_no) {
# 		$count ++;
# 		my $new_domain_no = $count;
# 		${$r_new_domain_no}{$gene}{$old_domain_no} = $new_domain_no; # new domain number
# 	    }
# 	}
#     }
# }

### Adjacency ###
sub get_adjacency_information_detail {
    my ($h_domain, $r_cluster_adjacency) = @_;

    for my $gene (keys %{$h_domain}) {
        my ($sp, $name) = decompose_gene_id($gene);
        my @domain = sort {$a<=>$b} keys %{${$h_domain}{$gene}};
        for (my $i=1; $i<@domain; $i++) {
            my $domain_before = $domain[$i-1];
            my $domain_after = $domain[$i];
            if ($domain_before + 1 == $domain_after) {
                my $cluster_before = ${$h_domain}{$gene}{$domain_before}{cluster};
                my $cluster_after = ${$h_domain}{$gene}{$domain_after}{cluster};
                ${$r_cluster_adjacency}{$cluster_before}{$cluster_after}{n}++;
                ${$r_cluster_adjacency}{$cluster_before}{$cluster_after}{org}{$sp} ++;
                ${$r_cluster_adjacency}{$cluster_before}{$cluster_after}{domain}{"$gene($domain_before)"} = 1;
                ${$r_cluster_adjacency}{$cluster_before}{$cluster_after}{domain}{"$gene($domain_after)"} = 1;
            }
        }
    }
}

sub get_adjacency_information {
    my ($h_domain, $r_cluster_adjacency) = @_;

    for my $gene (keys %{$h_domain}) {
        my @domain = sort {$a<=>$b} keys %{${$h_domain}{$gene}};
        for (my $i=1; $i<@domain; $i++) {
            my $domain_before = $domain[$i-1];
            my $domain_after = $domain[$i];
            if ($domain_before + 1 == $domain_after) {
                my $cluster_before = ${$h_domain}{$gene}{$domain_before}{cluster};
                my $cluster_after = ${$h_domain}{$gene}{$domain_after}{cluster};
                ${$r_cluster_adjacency}{$cluster_before}{$cluster_after}++;
            }
        }
    }
}

sub get_adjacent_organism {
    my ($h_domain, $r_cluster_adjacency) = @_;

    for my $gene (keys %{$h_domain}) {
        my ($sp, $name) = decompose_gene_id($gene);
        my @domain = sort {$a<=>$b} keys %{${$h_domain}{$gene}};
        for (my $i=1; $i<@domain; $i++) {
            my $domain_before = $domain[$i-1];
            my $domain_after = $domain[$i];
            if ($domain_before + 1 == $domain_after) {
                my $cluster_before = ${$h_domain}{$gene}{$domain_before}{cluster};
                my $cluster_after = ${$h_domain}{$gene}{$domain_after}{cluster};
                ${$r_cluster_adjacency}{$cluster_before}{$cluster_after}{$sp} ++;
            }
        }
    }
}

sub get_patterns {
    my ($r_cluster_adjacency) = @_;

    my $n_loop = 0;
    my $n_reverse = 0;
    my %neighbor = ();
    for my $cluster1 (keys %{$r_cluster_adjacency}) {
        for my $cluster2 (keys %{${$r_cluster_adjacency}{$cluster1}}) {
            if ($cluster1 eq $cluster2) {
                $n_loop++;
            } else {
                if ($neighbor{right}{$cluster2}{$cluster1}) {
                    $n_reverse++;
                }
                $neighbor{right}{$cluster1}{$cluster2} = 1;
                $neighbor{left}{$cluster2}{$cluster1} = 1;
            }
        }
    }

    my @out = ();

    my @node_open = sort {$a cmp $b} keys(%{$neighbor{right}});
    for my $node (@node_open) {
        my $n_branch = keys %{$neighbor{right}{$node}};
        if ($n_branch > 2) {
            push @out, "<$n_branch";
        } elsif ($n_branch == 2) {
            push @out, "<";
        }
    }

    my @node_close = sort {$a cmp $b} keys(%{$neighbor{left}});
    for my $node (@node_close) {
        my $n_branch = keys %{$neighbor{left}{$node}};
        if ($n_branch > 2) {
            push @out, "$n_branch>";
        } elsif ($n_branch == 2) {
            push @out, ">";
        }
    }

    for (my $i=0; $i<$n_loop; $i++) {
        push @out, "O";
    }

    for (my $i=0; $i<$n_reverse; $i++) {
        push @out, "X";
    }

    if (!@out) {
        push @out, "-";
    }
    
    return @out;
}

sub get_patterns_detail {
    my ($r_cluster_adjacency, @cluster) = @_;

    my %neighbor = ();
    for my $cluster1 (keys %{$r_cluster_adjacency}) {
        for my $cluster2 (keys %{${$r_cluster_adjacency}{$cluster1}}) {
            if ($cluster1 eq $cluster2) {
            } else {
                $neighbor{$cluster1}{right}{$cluster2} = 1;
                $neighbor{$cluster2}{left}{$cluster1} = 1;
            }
        }
    }

    my @out = ();
    for my $cluster (@cluster) {
        my @left = keys %{$neighbor{$cluster}{left}};
        my $left = "";
        if (@left == 1) {
            $left = "-";
        } elsif (@left == 2) {
            $left = ">";
        } elsif (@left >= 3) {
            $left = scalar(@left) . ">";
        }
        my @right = keys %{$neighbor{$cluster}{right}};
        my $right = "";
        if (@right == 1) {
            $right = "-";
        } elsif (@right == 2) {
            $right = "<";
        } elsif (@right >= 3) {
            $right = "<" . scalar(@right);
        }
        push @out, $left . "o" . $right;
    }
    
    return @out;
}

sub extract_most_adjacent_cluster_pair {
    my ($r_cluster_adjacency) = @_;

    my ($cluster1, $cluster2);
    my $adjacency_max;
    for my $c1 (keys %{$r_cluster_adjacency}) {
        for my $c2 (keys %{${$r_cluster_adjacency}{$c1}}) {
            my $adjacency = ${$r_cluster_adjacency}{$c1}{$c2};
            if (! defined $adjacency_max or $adjacency > $adjacency_max) {
                $adjacency_max = $adjacency;
                $cluster1 = $c1;
                $cluster2 = $c2;
            }
        }
    }
    
    return ($cluster1, $cluster2);
}

sub print_in_simple_text {
    my ($r_cluster, $r_cluster_count, $r_cluster_adjacency) = @_;

    for my $cluster (@{$r_cluster}) {
        print "$cluster(${$r_cluster_count}{$cluster})\n";
    }
    print "\n";
    for my $cluster1 (keys %{$r_cluster_adjacency}) {
        for my $cluster2 (keys %{${$r_cluster_adjacency}{$cluster1}}) {
            print "$cluster1-$cluster2(${$r_cluster_adjacency}{$cluster1}{$cluster2})\n";
        }
    }
    print "\n";
}

sub print_in_text {
    my ($r_cluster_adjacency, $r_annotation, @cluster) = @_;

    my %neighbor = ();
    for my $cluster1 (keys %{$r_cluster_adjacency}) {
        for my $cluster2 (keys %{${$r_cluster_adjacency}{$cluster1}}) {
            if ($cluster1 eq $cluster2) {
            } else {
                $neighbor{$cluster1}{right}{$cluster2} = 1;
                $neighbor{$cluster2}{left}{$cluster1} = 1;
            }
        }
    }

    my @out = ();
    for my $cluster (@cluster) {
        my @left = keys %{$neighbor{$cluster}{left}};
        my $left = "";
        if (@left == 1) {
            $left = "-";
        } elsif (@left == 2) {
            $left = ">";
        } elsif (@left >= 3) {
            $left = scalar(@left) . ">";
        }
        my @right = keys %{$neighbor{$cluster}{right}};
        my $right = "";
        if (@right == 1) {
            $right = "-";
        } elsif (@right == 2) {
            $right = "<";
        } elsif (@right >= 3) {
            $right = "<" . scalar(@right);
        }
        print $cluster, "\t", $left . "o" . $right, "\t", ${$r_annotation}{$cluster}, "\t", join(",", @left), "\t", join(",", @right), "\n";
    }
}

sub print_module_members {
    my ($h_module_to_clusters, %opt) = @_;
    
    for my $module (sort {$a<=>$b} keys(%{$h_module_to_clusters})) {
        my @cluster = sort {$a cmp $b} @{${$h_module_to_clusters}{$module}};
        if ($opt{m}) {
            print $module, "\t";
        }
        if ($opt{n}) {
            print scalar(@cluster), "\t";
        }
        print join(",", @cluster), "\n";
    }
}

sub link_modules {
    my ($h_gene_to_clusters, $h_module_to_cluster) = @_;

    my %cluster_to_module;
    for my $key (keys %{$h_gene_to_clusters}) {
        my @cluster = uniq @{${$h_gene_to_clusters}{$key}};
        if (@cluster <= 1) {
            delete ${$h_gene_to_clusters}{$key};
            next;
        }
        for my $cluster (@cluster) {
            $cluster_to_module{$cluster} = $cluster;
            ${$h_module_to_cluster}{$cluster} = [$cluster];
        }
    }

    for my $gene (keys %{$h_gene_to_clusters}) {
        my @cluster = uniq @{${$h_gene_to_clusters}{$gene}};

        my @related_module = ();
        for my $cluster (@cluster) {
            push @related_module, $cluster_to_module{$cluster};
        }
        @related_module = uniq(@related_module);

        my @related_cluster = ();
        for my $related_module (@related_module) {
            push @related_cluster, @{${$h_module_to_cluster}{$related_module}};
            delete ${$h_module_to_cluster}{$related_module}
        }
        @related_cluster = uniq(@related_cluster);

        my ($rep_module) = sort {$a <=> $b} @related_module;
        ${$h_module_to_cluster}{$rep_module} = \@related_cluster;
        for my $related_cluster (@related_cluster) {
            $cluster_to_module{$related_cluster} = $rep_module;
        }
    }
}

### Move ###
sub get_candidates {
    my ($tmp_cluster, $cluster) = @_;

    my @candidates = `cat $tmp_cluster | grep '^$cluster '`;
    chomp(@candidates);

    my @out = ();
    for my $candidate (@candidates) {
        if ($candidate =~ /^\S+ \S+ 0 \d+ \d+$/) {
        } elsif ($candidate =~ /^\S+ \S+ \S+ \d+ \d+$/) {
            push @out, $candidate;
        } else {
            die "$candidate";
        }
    }

    return @out;
}

sub get_candidate_neighbor {
    my ($dclst, $candidate) = @_;

    my ($cluster, $gene, $domain) = split(" ", $candidate);

    my @neighbor = `cat $dclst | grep -P '^\\S+ $gene \\S+ \\d+ \\d+\$' | grep -v '^$cluster '`;
    chomp(@neighbor);

    return @neighbor;
}

sub only_adjacent_neighbor {
    my ($candidate, @neighbor) = @_;

    my ($cluster, $gene, $domain) = split(" ", $candidate);
    
    my @out = ();
    for my $neighbor (@neighbor) {
        my ($cluster_n, $gene_n, $domain_n) = split(" ", $neighbor);
        if ($domain_n == $domain - 1 or $domain_n == $domain + 1) {
            push @out, $neighbor;
        }
    }

    return @out;
}

sub move_one {
    my ($dclst, $candidate, $neighbor) = @_;
    
    my @dclst = split("\n", $dclst);
    chomp(@dclst);

    my ($cluster, $gene, $domain) = split(" ", $candidate);
    my ($cluster_n) = split(" ", $neighbor);
    for (my $i=0; $i<@dclst; $i++) {
        $dclst[$i] =~ s/^$cluster ($gene $domain \S+ \S+)$/$cluster_n $1/;
    }

    return join("\n",@dclst)."\n";
}

sub calc_gain_by_move {
    my ($dclst_file, $candidate, $neighbor, $flg_uniq) = @_;

    # restrict the calculation
    my ($cluster1) = split(/\s+/, $candidate);
    my ($cluster2) = split(/\s+/, $neighbor);
    my $dclst_of_cluster_pairs = extract_dclst($dclst_file, $cluster1, $cluster2);
    
    my $score_before;
    my $score_after;
    if ($flg_uniq) {
        $score_before = score_dclst_one_alignment($dclst_of_cluster_pairs, uniq => 1); # old; should use *_c
        $score_after = score_dclst_one_alignment(move_one($dclst_of_cluster_pairs, $candidate, $neighbor), uniq => 1); # old; should use *_c
    } else {
        $score_before = score_dclst($dclst_of_cluster_pairs);
        $score_after = score_dclst(move_one($dclst_of_cluster_pairs, $candidate, $neighbor));
    }
    return $score_after - $score_before;
}

sub partial_move {
    my ($dclst, $cluster, $flg_uniq, $adjacent_neighbor_optional) = @_;

    save_contents($dclst, $TMP_MOVE_DCLST);

    print STDERR "[$cluster]\n";
    my %neighbor_to_join = ();
    for my $candidate (get_candidates($TMP_MOVE_DCLST, $cluster)) {
        my @neighbor = get_candidate_neighbor($TMP_MOVE_DCLST, $candidate);
        if ($adjacent_neighbor_optional) {
            @neighbor = only_adjacent_neighbor($candidate, @neighbor); # domain numbers should be integer
        }
        my $max_gain = 0;
        for my $neighbor (@neighbor) {
            my $gain = calc_gain_by_move($TMP_MOVE_DCLST, $candidate, $neighbor, $flg_uniq);
            print STDERR "TEST $candidate\t->\t$neighbor\tGAIN = $gain\n";
            if ($gain > $max_gain) {
                $max_gain = $gain;
                $neighbor_to_join{$candidate} = $neighbor;
            }
        }
    }
    
    for my $candidate (sort {$a cmp $b} keys %neighbor_to_join) {
        print STDERR "MOVE $candidate\t->\t$neighbor_to_join{$candidate}\n";
        $dclst = move_one($dclst, $candidate, $neighbor_to_join{$candidate}); # update dclst table
    }

    return $dclst;
}

sub concat_members_in_cluster {
    my ($r_member) = @_;

    for my $cluster (keys %{$r_member}) {
        for my $gene (keys %{${$r_member}{$cluster}}) {
            my @domains = sort {$a<=>$b} keys(%{${$r_member}{$cluster}{$gene}});
            for (my $i=$#domains; $i>=1; $i--) {
                if (${$r_member}{$cluster}{$gene}{$domains[$i-1]}{end} + 1 >= ${$r_member}{$cluster}{$gene}{$domains[$i]}{start}) { # BUG?
                    if (${$r_member}{$cluster}{$gene}{$domains[$i-1]}{start} > ${$r_member}{$cluster}{$gene}{$domains[$i]}{start}) {
                        die "$cluster $gene";
                    }
                    ${$r_member}{$cluster}{$gene}{$domains[$i-1]}{end} = ${$r_member}{$cluster}{$gene}{$domains[$i]}{end};
                    delete ${$r_member}{$cluster}{$gene}{$domains[$i]};
                }
            }
        }
    }
}

### boundary move ###
sub boundary_move {
    my ($r_a, $r_gene_idx, $r_domain, $r_get_pos, $j_boundary, $r_gene_to_change, $r_domain_to_change) = @_;

    for (my $g=0; $g<@{$r_gene_to_change}; $g++) {
        my $gene = ${$r_gene_to_change}[$g];
        my $i = ${$r_gene_idx}{$gene};
        my $j = $j_boundary;
        while (${$r_a}[$i][$j] eq '-') {
            $j ++;
        }
        my $domain = ${$r_domain_to_change}[$g];
        ${$r_domain}{$gene}{$domain-1}{end} = ${$r_get_pos}[$i][$j] - 1;
        ${$r_domain}{$gene}{$domain}{begin} = ${$r_get_pos}[$i][$j];
    }
}

sub boundary_move_gene {
    my ($r_a, $r_domain, $r_get_pos, $i, $gene, $domain, $j_boundary) = @_;

    my $j = $j_boundary;
    while (${$r_a}[$i][$j] and ${$r_a}[$i][$j] eq '-') {
        $j ++;
    }

    my $pos = ${$r_get_pos}[$i][$j];
    if ($pos) {
        print STDERR "pos= ", $pos, "\n";
        ${$r_domain}{$gene}{$domain-1}{end} = $pos - 1;
        ${$r_domain}{$gene}{$domain}{begin} = $pos;
        return $pos;
    }
}

### Refine ###
sub update_domain {
    my ($r_d_row, $from_pos, $to_pos, $begin_pos, $end_pos, $r_get_j_row) = @_;

    my $begin_j = ${$r_get_j_row}{$begin_pos};
    my $end_j = ${$r_get_j_row}{$end_pos};
    my $from_j = ${$r_get_j_row}{$from_pos};
    my $to_j = ${$r_get_j_row}{$to_pos};

    for (my $j=$from_j; $j<=$to_j; $j++) {
        if ($j >= $begin_j && $j <= $end_j) {
            ${$r_d_row}[$j] = 1;
        } else {
            ${$r_d_row}[$j] = 0;
        }
    }
}

sub update_domain_mapping {
    my ($r_d_row, $begin_pos, $end_pos, $r_get_j_row, $cluster) = @_;

    my $begin_j = ${$r_get_j_row}{$begin_pos};
    my $end_j = ${$r_get_j_row}{$end_pos};

    for (my $j=$begin_j; $j<=$end_j; $j++) {
        ${$r_d_row}[$j] = $cluster;
    }
}

sub update_domain_mapping_local {
    my ($r_d_row, $from_pos, $to_pos, $begin_pos, $end_pos, $r_get_j_row, $cluster) = @_;

    my $from_j = ${$r_get_j_row}{$from_pos};
    my $to_j = ${$r_get_j_row}{$to_pos};
    my $begin_j = ${$r_get_j_row}{$begin_pos};
    my $end_j = ${$r_get_j_row}{$end_pos};

    for (my $j=$from_j; $j<=$to_j; $j++) {
        if ($j >= $begin_j && $j <= $end_j) {
            ${$r_d_row}[$j] = $cluster;
        }
    }
}

sub get_clusters_subset {
    my ($r_cluster_members, @cluster) = @_;
    
    my $dclst = "";
    my $merged_members = "";

    for my $cluster (@cluster) {
        if (${$r_cluster_members}{$cluster}) {
            $merged_members .= ${$r_cluster_members}{$cluster};
            my $dclst_tmp = ${$r_cluster_members}{$cluster};
            $dclst_tmp =~ s/^(\S)/$cluster $1/gm;
            $dclst .= $dclst_tmp;
        } else {
            die "$cluster not contained in the input.";
        }
    }

    return ($dclst, $merged_members);
}

sub get_new_cluster_ids {
    my ($r_old2new_cluster, @cluster) = @_;

    my @updated_cluster = ();
    my $flg_changed = 0;
    for my $cluster (@cluster) {
        if (${$r_old2new_cluster}{$cluster}) {
            my $new_cluster = ${$r_old2new_cluster}{$cluster};
            print STDERR "[$cluster] is changed to [$new_cluster].\n";
            push @updated_cluster, $new_cluster;
            $flg_changed = 1;
        } else {
            push @updated_cluster, $cluster;
        }
    }
    @updated_cluster = uniq(@updated_cluster);

    return($flg_changed, @updated_cluster);
}

sub update_clusters {
    my ($r_cluster_members_txt, $r_member, $r_old2new_cluster, $r_new2old_cluster, @target_cluster) = @_;

    my $new_cluster = choose_new_cluster_id($r_cluster_members_txt, @target_cluster);

    replace_clusters($r_cluster_members_txt, @target_cluster);

    for my $target_cluster (@target_cluster) {
        for my $gene (keys %{${$r_member}{$target_cluster}}) {
            for my $domain (keys %{${$r_member}{$target_cluster}{$gene}}) {
                ${$r_member}{$new_cluster}{$gene}{$domain}{start} = ${$r_member}{$target_cluster}{$gene}{$domain}{start};
                ${$r_member}{$new_cluster}{$gene}{$domain}{end} = ${$r_member}{$target_cluster}{$gene}{$domain}{end};
            }
        }
    }

    # Update mappings between old/new IDs
    my @original_cluster = ();
    for my $target_cluster (@target_cluster) {
        if (${$r_new2old_cluster}{$target_cluster}) {
            push @original_cluster, @{${$r_new2old_cluster}{$target_cluster}};
        }
    }
    @original_cluster = uniq(@original_cluster, @target_cluster);

    @{${$r_new2old_cluster}{$new_cluster}} = @original_cluster;
    for my $original_cluster (@original_cluster) {
        ${$r_old2new_cluster}{$original_cluster} = $new_cluster;
    }

}

sub replace_clusters {
    my ($r_cluster_members_txt, @target_cluster) = @_;

    my ($dclst, $merged_members) = get_clusters_subset($r_cluster_members_txt, @target_cluster);

    # Delete
    for my $target_cluster (@target_cluster) {
        delete ${$r_cluster_members_txt}{$target_cluster};
    }

    # Insert
    my $new_cluster = choose_new_cluster_id($r_cluster_members_txt, @target_cluster);
    ${$r_cluster_members_txt}{$new_cluster} = $merged_members;
}

sub choose_new_cluster_id {
    my ($r_cluster_members, @cluster) = @_;

    my $flg_int = 1;
    my %prefix = ();
    for my $cluster (@cluster) {
        if ($cluster =~ /^(\d+)/) {
            my $int = $1;
            $prefix{$cluster} = $int;
        } else {
            $flg_int = 0;
        }
    }

    my @cluster_sorted = ();
    if ($flg_int) {
        @cluster_sorted = sort {$prefix{$a} <=> $prefix{$b}} @cluster;
    } else {
        @cluster_sorted = sort {$a cmp $b} @cluster;
    }

    my $new_cluster_id = $cluster_sorted[0] . "a";

    if (${$r_cluster_members}{$new_cluster_id}) {
        die; # possibly bug?: duplicated IDs might be produced in the case of string ?
    }

    return $new_cluster_id;
}

1;
