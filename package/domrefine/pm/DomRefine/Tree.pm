package DomRefine::Tree;
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(new_tree phylogenetic_tree find_min_tree find_min_tree_fast tree_children
	     save_support_values my_reroot_at_midpoint
	     check_species_overlap get_sub_tree_leaves sum_of_duplication get_sub_tree_branch_length
	     put_support_values put_node_labels get_species_from_leaves
	     print_tree tree_node tree_children_dist tree_children_o11
             reset_bar tree_padding
	     );

use strict;
use Bio::TreeIO;
use DomRefine::General;
use DomRefine::Align;
use Digest::MD5 qw(md5_hex);
use DomRefine::Read;

################################################################################
### Public functions ###########################################################
################################################################################

sub phylogenetic_tree {
    my ($tmp_alignment, $tmp_tree, $tmp_tree_log, %opt) = @_;
    
    my $program = $ENV{DOMREFINE_TREE} || die;

    escape_fasta_headers($tmp_alignment, %opt);

    my $cache_file = cache_file_path($ENV{DOMREFINE_CACHE}, md5_hex(fasta_to_sorted_tsv($tmp_alignment)), ".$program");

    if ($ENV{DOMREFINE_READ_TREE_CACHE} and found_tree_cache($tmp_alignment, $cache_file)) {
	print STDERR " found tree cache $cache_file\n";
	system "cp $cache_file $tmp_tree";
	return;
    }

    my $start_time = time;

    my $command = "cat $tmp_alignment | $program > $tmp_tree 2> $tmp_tree_log";
    print STDERR "command: $command\n";
    my ($n_seq, $n_aa, $mean_len) = seq_file_stat($tmp_alignment);
    print STDERR "$program: n_seq=$n_seq len=$mean_len ";
    system "$command";

    my $end_time = time;
    printf STDERR "%d sec\n", $end_time - $start_time;

    if ($ENV{DOMREFINE_WRITE_TREE_CACHE}) {
	print STDERR " saving tree cache $cache_file\n";
	print STDERR " saving tree cache $cache_file.err\n";
	system "cp -a $tmp_tree $cache_file";
	system "cp -a $tmp_tree_log $cache_file.err";
    }
}

sub found_tree_cache {
    my ($alignment_file, $tree_cache) = @_;

    if (! -f $tree_cache) {
	return 0;
    }

    my @ids_in_alignment = ();
    open(ALIGNMENT_FILE, $alignment_file) || die;
    while (my $line = <ALIGNMENT_FILE>) {
	chomp($line);
	if ($line =~ /^>(\S+)$/) {
	    push @ids_in_alignment, $1;
	}
    }
    close(ALIGNMENT_FILE);
    @ids_in_alignment = sort @ids_in_alignment;
    my $ids_in_alignment = join("\n", @ids_in_alignment);
    
    open(TREE_CACHE, $tree_cache) || die;
    my $tree = <TREE_CACHE>;
    close(TREE_CACHE);
    $tree =~ s/:/\n:/g;
    $tree =~ s/,/\n/g;
    $tree =~ s/\(//g;
    $tree =~ s/:.*\n//g;
    my @ids_in_tree = split("\n", $tree);
    @ids_in_tree = sort @ids_in_tree;
    my $ids_in_tree = join("\n", @ids_in_tree);
    
    if ($ids_in_alignment ne $ids_in_tree) {
    	print STDERR "WARNING: inconsistent tree_cache $tree_cache\n";
    	return 0;
    }

    return 1;
}

sub node_info {
    my ($node, $padding, $r_height, $r_inverse_height, $r_current_min_height, $r_min_tree_node) = @_;

    my $id = $node->internal_id;
    my $length = $node->branch_length || 0;
    my $height = $node->height || 0;
    my $depth = $node->depth;
    ${$r_height}{$id} = $height;

    ### Calc minimum height ###
    my $min_height;
    my $diff_child_height = abs(${$r_height}{$id} - ${$r_inverse_height}{$id});
    # if ($length >= $diff_child_height) {
    # 	$min_height = (${$r_height}{$id} + ${$r_inverse_height}{$id} + $length) / 2; # balanced tee
    # } else {
    # 	$min_height = max(${$r_height}{$id}, ${$r_inverse_height}{$id}); # unblanced tree
    # }
    $min_height = max(${$r_height}{$id}, ${$r_inverse_height}{$id}) + $length / 2;

    if ($padding eq "") {
	# root
    } elsif (! defined ${$r_current_min_height}) {
	${$r_current_min_height} = $min_height;
	${$r_min_tree_node} = $node;
    } elsif ($min_height < ${$r_current_min_height}) {
	${$r_current_min_height} = $min_height;
	${$r_min_tree_node} = $node;
    }

    # print STDERR "$padding($length) $id (d=$depth h=$height) $min_height ${$r_current_min_height}\n";

    ### Recursion ###
    if ($node->is_Leaf) {
	return;
    }
    for my $child ($node->each_Descendent) {
	node_info($child, "$padding    ", $r_height, $r_inverse_height, $r_current_min_height, $r_min_tree_node);
    }
}

sub calc_inverse_heights {
    my ($node, $parent_inverse_height, $r_inverse_height) = @_;

    if ($node->is_Leaf) {
	return;
    }

    my $node_id = $node->internal_id;
    my $length = $node->branch_length || 0;
    my @child = $node->each_Descendent;
    my @id;
    my @height;
    my @length;
    for (my $i=0; $i<@child; $i++) {
	$id[$i] = $child[$i]->internal_id;
	$height[$i] = $child[$i]->height;
	$length[$i] = $child[$i]->branch_length;
    }

    # if (@child == 2) {
    # 	${$r_inverse_height}{$id[0]} = max($parent_inverse_height + $length, $length[1] + $height[1]);
    # 	${$r_inverse_height}{$id[1]} = max($parent_inverse_height + $length, $length[0] + $height[0]);
    # 	calc_inverse_heights($child[0], ${$r_inverse_height}{$id[0]}, $r_inverse_height);
    # 	calc_inverse_heights($child[1], ${$r_inverse_height}{$id[1]}, $r_inverse_height);
    # } elsif (@child == 3) {
    # 	${$r_inverse_height}{$id[0]} = max($height[1] + $length[1], $height[2] + $length[2]);
    # 	${$r_inverse_height}{$id[1]} = max($height[0] + $length[0], $height[2] + $length[2]);
    # 	${$r_inverse_height}{$id[2]} = max($height[0] + $length[0], $height[1] + $length[1]);
    # 	for (my $i=0; $i<@child; $i++) {
    # 	    calc_inverse_heights($child[$i], ${$r_inverse_height}{$id[$i]}, $r_inverse_height);
    # 	}
    # } else {
    # 	die;
    # }

    if (@child >= 2) {
	my @length_plus_height = ();
	for (my $i=0; $i<@child; $i++) {
	    $length_plus_height[$i] = $length[$i] + $height[$i];
	}
	for (my $i=0; $i<@child; $i++) {
	    my @optional_edges;
	    for (my $j=0; $j<@child; $j++) {
		if ($i != $j) {
		    push @optional_edges, $length_plus_height[$j];
		}
	    }
	    ${$r_inverse_height}{$id[$i]} = max($parent_inverse_height + $length, @optional_edges);
	    calc_inverse_heights($child[$i], ${$r_inverse_height}{$id[$i]}, $r_inverse_height);
	}
    } else {
	die;
    }

}

sub reroot_with_bootstrap {
    my ($tree,$new_root) = @_;
    unless (defined $new_root && $new_root->isa("Bio::Tree::NodeI")) {
        $tree->warn("Must provide a valid Bio::Tree::NodeI when rerooting");
        return 0;
    }

    my $old_root = $tree->get_root_node;
    if( $new_root == $old_root ) {
        $tree->warn("Node requested for reroot is already the root node!");
        return 0;
    }
    my $anc = $new_root->ancestor;
    unless( $anc ) {
        # this is already the root
        $tree->warn("Node requested for reroot is already the root node!");
        return 0;
    }
    my $tmp_node = $new_root->create_node_on_branch(-position=>0,-force=>1);
    # reverse the ancestor & children pointers
    my $former_anc = $tmp_node->ancestor;
    my $new_first_child;
    my @path_from_oldroot = ($tree->get_lineage_nodes($tmp_node), $tmp_node);
    for (my $i = 0; $i < $#path_from_oldroot; $i++) {
        my $current = $path_from_oldroot[$i];
        my $next = $path_from_oldroot[$i + 1];
        $current->remove_Descendent($next);
        $current->branch_length($next->branch_length);

	### move IDs (bootstrap values) between old and new roots ###
	if ($i == $#path_from_oldroot - 1) {
	    $new_first_child = $current;
	    $current->id("");
	} else {
	    $current->id($next->id) if defined $next->id;
	}

        $next->remove_tag('B');
        $next->add_Descendent($current);
    }

    $new_root->add_Descendent($former_anc);
    $tmp_node->remove_Descendent($former_anc);
    
    $tmp_node = undef;
    $new_root->branch_length(undef);

    $old_root = undef;
    $tree->set_root_node($new_root);

    ### copy IDs (bootstrap values) between the new children ###
    # my @new_child = $new_root->each_Descendent;
    # if (@new_child != 2) {
    # 	die scalar(@new_child);
    # }

    # my $new_second_child;
    # if ($new_child[0] eq $new_first_child and $new_child[1] ne $new_first_child) {
	# $new_second_child = $new_child[1];
    # } elsif ($new_child[0] ne $new_first_child and $new_child[1] eq $new_first_child) {
	# $new_second_child = $new_child[0];
    # } else {
	# die;
    # }

    # if (! $new_second_child->is_Leaf) {
	# $new_first_child->id($new_second_child->id);
    # }

}

sub find_min_tree_fast {
    my ($tree) = @_;

    my $root = $tree->get_root_node;
    my $root_id = $root->internal_id;

    my %inverse_height;
    $inverse_height{$root_id} = 0;
    calc_inverse_heights($root, $inverse_height{$root_id}, \%inverse_height);

    my %height = ();
    my $current_min_height;
    my $min_tree_node;
    node_info($root, "", \%height, \%inverse_height, \$current_min_height, \$min_tree_node);

    # print STDERR "current_min_height = $current_min_height\n";
    # print STDERR "min_tree_node: ", $min_tree_node->internal_id, "\n";

    ### Slow here ? ###
    my @node = $root->get_all_Descendents;
    my $min_tree_idx;
    for (my $i=0; $i<@node; $i++) {
	if ($node[$i]->internal_id eq $min_tree_node->internal_id) {
	    $min_tree_idx = $i;
	}
    }

    my $midpoint = $min_tree_node->create_node_on_branch(-FRACTION=>0.5, -FORCE=>1);
    # $tree->reroot($midpoint);
    reroot_with_bootstrap($tree, $midpoint);

    return $min_tree_idx;
}

sub find_min_tree {
    my ($tree_file, $r_nodes) = @_;

    my @height = ();
    for (my $i=0; $i<@{$r_nodes}; $i++) {
	my $internal_id = ${$r_nodes}[$i]->internal_id;
	my $branch_length = ${$r_nodes}[$i]->branch_length;
	my $tree = new_tree($tree_file);
	my @nodes = $tree->get_root_node->get_all_Descendents;
	my_reroot_at_midpoint(\$tree, $nodes[$i]);

	my @child = $tree->get_root_node->each_Descendent;
	my $sum_child_height = -1;
	my $child1_height = -1;
	my $child2_height = -1;
	my $min_height = -1;
	if (@child == 2) {
	    $child1_height = $child[0]->height;
	    $child2_height = $child[1]->height;
	    $sum_child_height = $child1_height + $child2_height;
	    my $diff_child_height = abs($child1_height - $child2_height);
	    if ($branch_length >= $diff_child_height) {
		$min_height = ($child1_height + $child2_height + $branch_length) / 2;
	    } else {
		$min_height = max($child1_height, $child2_height);
	    }
	} else {
	    print STDERR "WARNING: child of new root == ", scalar(@child) ,"\n";
	}

	my $height = sprintf("%.5f", $tree->get_root_node->height);
	# $height[$i] = $height;
	$height[$i] = $min_height;
	# print STDERR "[$i] internal=$internal_id branch=$branch_length\th=$height, $min_height", "\t", ${$r_nodes}[$i]->to_string, "\n"; # print for debug
    }

    if (@height) {
	my $min_i = min_i(@height);
	print STDERR "minimum height: [$min_i]\n";
	return $min_i
    } else {
	print STDERR "WARNING: cannot find heights nor minimum heights; thus, return node[0].\n";
	return 0;
    }
}

sub tree_children {
    my ($node, $offset, @bar) = @_;

    my @child = $node->each_Descendent;
    if (@child == 2) {
	my @bar1 = reset_bar(-$offset, @bar, -($offset+2));
	my @bar2 = reset_bar($offset, @bar, $offset+2);
	return
	    tree_children($child[0], $offset+2, @bar1) .
	    tree_padding($offset, @bar) . "+-|\n" .
	    tree_children($child[1], $offset+2, @bar2);
    } elsif (@child == 0) {
	return
	    tree_padding($offset, @bar) . "+- " . tree_node($node) . "\n";
    } else {
	print STDERR scalar(@child), " children.\n";
	die;
    }
}

sub sum_of_duplication {
    my ($r_gene_species, $r_tree) = @_;

    my @node = grep {! $_->is_Leaf} ${$r_tree}->get_nodes();

    my $sum_of_duplication = 0;
    my $sum_of_sp_disappeared = 0;
    for my $node (@node) {
	my ($is_diuplication, $sp_disappeared) = check_duplication($r_gene_species, $node);
	if ($is_diuplication) {
	    $sum_of_duplication += $is_diuplication;
	    $sum_of_sp_disappeared += $sp_disappeared;
	}
    }

    return ($sum_of_duplication, $sum_of_sp_disappeared);
}

sub save_support_values {
    my ($r_tree, $r_support_value) = @_;
    
    my @nodes = grep {! $_->is_Leaf} ${$r_tree}->get_root_node->get_all_Descendents;
    
    for (my $i=0; $i<@nodes; $i++) {
	my $id1 = $nodes[$i]->internal_id;
	my $id2 = $nodes[$i]->ancestor->internal_id;
	my $support_value = $nodes[$i]->id;
	
	${$r_support_value}[$id1][$id2] = $support_value;
	${$r_support_value}[$id2][$id1] = $support_value;
    }
    
}

sub put_support_values {
    my ($r_tree, $r_support_value) = @_;

    my $root_node = ${$r_tree}->get_root_node;
    my @nodes = grep {! $_->is_Leaf} $root_node->get_all_Descendents;
    for (my $i=0; $i<@nodes; $i++) {
	my $id1 = $nodes[$i]->internal_id;
	my $id2 = $nodes[$i]->ancestor->internal_id;
	$nodes[$i]->id("");
# 	if (${$r_support_value}[$id1][$id2]) {
# 	    $nodes[$i]->id(${$r_support_value}[$id1][$id2]);
# 	}
    }
    my ($node1, $node2) = $root_node->each_Descendent;
    my $id1 = $node1->internal_id;
    my $id2 = $node2->internal_id;
    $node1->id(${$r_support_value}[$id1][$id2]);
    $node2->id(${$r_support_value}[$id1][$id2]);
}

sub check_species_overlap {
    my ($r_gene_species, $tree_file, $i, $node1, $node2) = @_;

    my $tree = new_tree($tree_file);
    my @nodes = $tree->get_root_node->get_all_Descendents;
    my $n_seq = grep {$_->is_Leaf} @nodes;
    my $branch_length = $nodes[$i]->branch_length;
    $tree->move_id_to_bootstrap;
    my $boot = $nodes[$i]->bootstrap || "";

    # calculation
    my_reroot_at_midpoint(\$tree, $nodes[$i]);
    my $root_node = $tree->get_root_node;
    my ($sub_tree1_node, $sub_tree2_node) = $root_node->each_Descendent;
    my @leaves1 = get_sub_tree_leaves($sub_tree1_node);
    my @leaves2 = get_sub_tree_leaves($sub_tree2_node);
#     print STDERR "leaves1 @leaves1\n";
#     print STDERR "leaves2 @leaves2\n";

    my @species1 = get_species_from_leaves($r_gene_species, @leaves1);
    my @species2 = get_species_from_leaves($r_gene_species, @leaves2);
#     print STDERR "species1 @species1\n";
#     print STDERR "species2 @species2\n";
    my @all_species = uniq(@species1, @species2);
    my @common_species = check_redundancy(@species1, @species2);

    my $sum_of_duplication = 0;
    my $sum_of_sp_disappeared = 0;
    ($sum_of_duplication, $sum_of_sp_disappeared) = sum_of_duplication($r_gene_species, \$tree);

    # print
    $node1 ||= "";
    $node2 ||= "";
    my $height = sprintf("%.5f", $root_node->height);
    my @len1 = get_sub_tree_branch_length($sub_tree1_node);
    my @len2 = get_sub_tree_branch_length($sub_tree2_node);
    my $len1 = mean(@len1) || 0;
    my $len2 = mean(@len2) || 0;
    $len1 = sprintf("%.5f", $len1);
    $len2 = sprintf("%.5f", $len2);
    my $len12 = mean(@len1, @len2);
    my $len_relative = 0;
    if ($len12) {
	$len_relative = $branch_length / $len12;
    }
    $len_relative = sprintf("%.5f", $len_relative);
    my $log_ratio_len = "";
    if ($len2 != 0 && $len1/$len2) {
	$log_ratio_len = log($len1/$len2)/log(2);
	$log_ratio_len = abs($log_ratio_len);
	$log_ratio_len = sprintf("%.5f", $log_ratio_len);
    }
    my $height1 = sprintf("%.5f", $sub_tree1_node->height);
    my $height2 = sprintf("%.5f", $sub_tree2_node->height);
    my $n1 = scalar(@leaves1);
    my $n2 = scalar(@leaves2);
    my $n_sp = scalar(@all_species);
    my $n_sp_common = scalar(@common_species);
    my $sp_overlap = @common_species/@all_species;
    $sp_overlap = sprintf("%.5f", $sp_overlap);
    my $sp_overlap1 = 0;
    if (@species1) {
	$sp_overlap1 =  @common_species/@species1;
    }
    my $sp_overlap2 = 0;
    if (@species2) {
	$sp_overlap2 =  @common_species/@species2;
    }
    my $max_sp_overlap_part = max($sp_overlap1, $sp_overlap2);

    my $min_height = ($height1 + $height2 + $branch_length) / 2;
    if ($min_height < $height1 and $min_height < $height2) {
	die;
    } elsif ($min_height < $height1) {
	$min_height = $height1;
    } elsif ($min_height < $height2) {
	$min_height = $height2;
    }

    my $detail = "i=$i($node1,$node2)\tb=$boot,\tl= $branch_length"
	. ", l_rel= $len_relative"
	. ", l1= $len1, l2= $len2"
	. ", |log2(l1/l2)|= $log_ratio_len,"
	. "\th= $min_height"
	# . "\th= $height"
	# . ", h1= $height1, h2= $height2"
	. ", n= $n_seq = $n1 + $n2, "
	. "o_sp=$n_sp_common/$n_sp=\t$sp_overlap\t, o_sp_part=\t$max_sp_overlap_part, n_dup=$sum_of_duplication, n_dis=$sum_of_sp_disappeared\n";

    return $sp_overlap, \@leaves1, \@leaves2, $detail
	# , $height
	, $min_height
	;
}

sub get_species_from_leaves {
    my ($r_gene_species, @leaves) = @_;
    
    my %species;
    for my $leaf (@leaves) {
	if ($leaf =~ /^\d+_([A-Za-z0-9]+)_(\S+)$/) {
	    my ($sp, $gene) = ($1, $2);
	    $species{$sp} = 1;
            expand_species_for_gene("$sp:$gene", $r_gene_species, \%species);
	} elsif ($leaf =~ /^([A-Za-z0-9]+)_(\S+)$/) {
	    my ($sp, $gene) = ($1, $2);
	    $species{$sp} = 1;
            expand_species_for_gene("$sp:$gene", $r_gene_species, \%species);
	} else {
	    die $leaf;
	}
    }

    return keys %species;
}

sub expand_species_for_gene {
    my ($gene, $r_gene_species, $r_species) = @_;
    if (${$r_gene_species}{$gene}) {
        my @expanded_species = keys %{${$r_gene_species}{$gene}};
        for my $e (@expanded_species) {
            ${$r_species}{$e} = 1;
        }
    } else {
        die "cannot get expanded species for gene $gene";
    }
}

sub get_sub_tree_leaves {
    my ($sub_tree_node) = @_;

    my @leaves = grep {$_->is_Leaf} ($sub_tree_node, $sub_tree_node->get_Descendents);
    my @leaves_id = ();
    for my $leaf (@leaves) {
	push @leaves_id, $leaf->id;
    }

    return @leaves_id;
}

sub get_sub_tree_branch_length {
    my ($sub_tree_node) = @_;

    my @node = $sub_tree_node->get_Descendents;

    my @len = ();
    for my $node (@node) {
    	push @len, $node->branch_length;
    }
    return @len;

    # my $sum = 0;
    # for my $node (@node) {
    # 	$sum += $node->branch_length;
    # }

    # my $mean = 0;
    # if (@node != 0) {
    # 	$mean = $sum/@node;
    # }

    # return sprintf("%.5f", $mean);
}

sub print_tree {
    my ($r_tree) = @_;

    my $tree_io = Bio::TreeIO->new(-format=>'newick');
    $tree_io->write_tree(${$r_tree});
}

sub new_tree {
    my ($tree_file) = @_;

    my $tree_io = Bio::TreeIO->new(-format => "newick", -file => "$tree_file");
    my $tree = $tree_io->next_tree;
    
    return $tree;
}

sub my_reroot_at_midpoint {
    my ($r_tree, $node) = @_;

    my $midpoint = $node->create_node_on_branch(-FRACTION=>0.5, -FORCE=>1);
    ${$r_tree}->reroot($midpoint);
}

sub tree_children_dist {
    my ($h_domain, $node, $offset, $r_bar, %opt) = @_;

    my $branch_len = 0;
    if ($node->branch_length) {
	$branch_len = $node->branch_length;
    }
    my @child = $node->each_Descendent;
    if (@child == 2) {
	my @bar1 = reset_bar(-$offset, @{$r_bar}, -($offset+2));
	my @bar2 = reset_bar($offset, @{$r_bar}, $offset+2);
	my ($tree_txt1, $height1) = tree_children_dist($h_domain, $child[0], $offset+2, \@bar1, %opt);
	my ($tree_txt2, $height2) = tree_children_dist($h_domain, $child[1], $offset+2, \@bar2, %opt);
	my $dist = $height1 + $height2;
	my $edge = tree_padding($offset, @{$r_bar}) . "+-| " . $dist * 100 .
	    # " ($branch_len)" .
	    "\n";
	return("${tree_txt1}${edge}${tree_txt2}", $branch_len + $dist/2, $dist);
    } elsif (@child == 1) {
	my ($tree_txt, $height, $dist) = tree_children_dist($h_domain, $child[0], $offset, $r_bar, %opt);
	print STDERR "WARNING: only one children.\n";
	return("$tree_txt", $branch_len + $height, $dist);
    } elsif (@child == 0) {
	my ($tree_txt) = tree_padding($offset, @{$r_bar}) . "+- " . 
	    tree_node_with_optional_domain_number($h_domain, $node, %opt) . 
	    # tree_node_with_domain_number($h_domain, $node) . 
	    # tree_node($node) . 
	    # " ($branch_len)" .
	    "\n";
	return("$tree_txt", $branch_len);
    } else {
	print STDERR "node ", $node->to_string, " has ", scalar(@child), " children.\n";
	for (my $i=0; $i<@child; $i++) {
	    if ($child[$i]->branch_length != 0) {
		print STDERR "child[$i] branch=", $child[$i]->branch_length, "\n";
		die;
	    }
	}
	my ($tree_txt, $height, $dist) = tree_children_dist_multi($h_domain, \@child, 0, $offset, $r_bar, %opt);
	return("$tree_txt", $branch_len + $height, $dist);
    }
}

################################################################################
### Private functions ##########################################################
################################################################################

sub tree_children_dist_multi {
    my ($h_domain, $r_node, $i, $offset, $r_bar, %opt) = @_;

    my @child = @{$r_node};

    my @bar1 = reset_bar(-$offset, @{$r_bar}, -($offset+2));
    my @bar2 = reset_bar($offset, @{$r_bar}, $offset+2);
    my ($tree_txt1, $height1) = tree_children_dist($h_domain, $child[$i], $offset+2, \@bar1, %opt);
    my ($tree_txt2, $height2);
    if ($i == @child - 2) {
	($tree_txt2, $height2) = tree_children_dist($h_domain, $child[$i+1], $offset+2, \@bar2, %opt);
    } else {
	($tree_txt2, $height2) = tree_children_dist_multi($h_domain, $r_node, $i+1, $offset+2, \@bar2, %opt);
    }
    my $dist = $height1 + $height2;
    my $edge = tree_padding($offset, @{$r_bar}) . "+-| " . $dist * 100 .
	# " ($branch_len)" .
	"\n";
    return("${tree_txt1}${edge}${tree_txt2}", 0, 0);
}

sub tree_children_o11_multi {
    my ($h_domain, $r_node, $i, $level, $parent_id, $max_id, %opt) = @_;

    my @child = @{$r_node};

    $max_id ++;
    my $id = $max_id;
    my ($tree_txt1, $height1, $dist1, $gene1) = tree_children_o11($h_domain, $child[$i], $level+1, $id, $max_id, %opt);
    my ($tree_txt2, $height2, $dist2, $gene2);
    if ($i == @child - 2) {
	($tree_txt2, $height2, $dist2, $gene2) = tree_children_o11($h_domain, $child[$i+1], $level+1, $id, $max_id, %opt);
    } else {
	($tree_txt2, $height2, $dist2, $gene2) = tree_children_o11_multi($h_domain, $r_node, $i+1, $level+1, $id, $max_id, %opt);
    }
    my $dist = $height1 + $height2;
    if ($level == 0) {
	$parent_id = -1;
    }
    return("I $id $parent_id $gene1 " . $dist * 100 . " " . $dist * 100 . "\n" .
	   $tree_txt1 . $tree_txt2, 0, 0, $gene1);
}

sub tree_children_o11 {
    my ($h_domain, $node, $level, $parent_id, $max_id, %opt) = @_;

    my $branch_len = 0;
    if ($node->branch_length) {
	$branch_len = $node->branch_length;
    }
    my $id = $node->internal_id;
    if (! $parent_id) {
	my $ancestor = $node->ancestor;
	if ($ancestor) {
	    $parent_id = $ancestor->internal_id;
	}
    }
    my @child = $node->each_Descendent;
    if (@child == 2) {
	my ($tree_txt1, $height1, $dist1, $gene1) = tree_children_o11($h_domain, $child[0], $level+1, $id, $max_id, %opt);
	my ($tree_txt2, $height2, $dist2, $gene2) = tree_children_o11($h_domain, $child[1], $level+1, $id, $max_id, %opt);
	my $dist = $height1 + $height2;
	if ($level == 0) {
	    $parent_id = -1;
	}
	return("I $id $parent_id $gene1 " . $dist * 100 . " " . $dist * 100 . "\n" .
	       $tree_txt1 . $tree_txt2, $branch_len + $dist/2, $dist, $gene1);
    } elsif (@child == 1) {
	my ($tree_txt, $height, $dist, $gene) = tree_children_o11($h_domain, $child[0], $level+1, $parent_id, $max_id, %opt);
	print STDERR "WARNING: only one children.\n";
	return($tree_txt, $branch_len + $height, $dist, $gene);
    } elsif (@child == 0) {
	my @leaf = tree_node_o11($h_domain, $node, %opt);
	my $gene = $leaf[0];
	return("L $id $parent_id @leaf\n", $branch_len, 0, $gene);
    } else {
	print STDERR "node ", $node->to_string, " has ", scalar(@child), " children.\n";
	for (my $i=0; $i<@child; $i++) {
	    if ($child[$i]->branch_length != 0) {
		print STDERR "child[$i] branch=", $child[$i]->branch_length, "\n";
		die;
	    }
	}
	my ($tree_txt, $height, $dist, $gene) = tree_children_o11_multi($h_domain, \@child, 0, $level+1, $parent_id, $max_id, %opt);
	return("$tree_txt", $branch_len + $height, $dist, $gene);
	die;
    }
}

sub tree_node_o11 {
    my ($h_domain, $node, %opt) = @_;
    
    my $string = $node->to_string;
    my ($gene) = split(":",  $string);

    if ($opt{DOMAIN}) {
	my ($organism, $name, $domain);
	if ($gene =~ /^(\d+)_(\S+?)_(\S+)$/) {
	    ($organism, $name, $domain) = ($2, $3, $1);
	} elsif ($gene =~ /^(\S+?)_(\S+)$/) {
	    ($organism, $name, $domain) = ($1, $2, 0);
	} else {
	    die;
	}
	if (${$h_domain}{"$organism:$name"}{$domain}) {
	    my $begin = ${$h_domain}{"$organism:$name"}{$domain}{begin};
	    my $end = ${$h_domain}{"$organism:$name"}{$domain}{end};
	    return("$organism:$name", $begin, $end, $domain);
	} else {
	    print STDERR "WARNING: no domain information";
	    return("$organism:$name", $domain);
	}
    } else {
	if ($gene =~ /^(\S+?)_(\S+)$/) {
	    my ($organism, $name) = ($1, $2);
	    return("$organism:$name");
	} else {
	    die;
	}
    }
}

sub tree_node_with_optional_domain_number {
    my ($h_domain, $node, %opt) = @_;
    
    my $string = $node->to_string;
    my ($gene) = split(":",  $string);

    if ($opt{DOMAIN}) {
	my ($organism, $name, $domain);
	if ($gene =~ /^(\d+)_(\S+?)_(\S+)$/) {
	    ($organism, $name, $domain) = ($2, $3, $1);
	} elsif ($gene =~ /^(\S+?)_(\S+)$/) {
	    ($organism, $name, $domain) = ($1, $2, 0);
	} else {
	    die;
	}
	if (${$h_domain}{"$organism:$name"}{$domain}) {
	    my $begin = ${$h_domain}{"$organism:$name"}{$domain}{begin};
	    my $end = ${$h_domain}{"$organism:$name"}{$domain}{end};
	    if ($domain eq "0") {
		return "$organism:$name $begin $end";
	    } else {
		return "$organism:$name($domain) $begin $end";
	    }
	} else {
	    print STDERR "WARNING: no domain information";
	    return "$organism:$name $domain";
	}
    } else {
	if ($gene =~ /^(\S+?)_(\S+)$/) {
	    my ($organism, $name) = ($1, $2);
	    return "$organism:$name";
	} else {
	    die;
	}
    }
}

sub tree_node_with_domain_number {
    my ($h_domain, $node) = @_;
    
    my $string = $node->to_string;
    my ($gene) = split(":",  $string);
    if ($gene =~ /^(\d+)_(\S+?)_(\S+)$/) {
	my ($domain, $organism, $name) = ($1, $2, $3);
	if (${$h_domain}{"$organism:$name"}{$domain}) {
	    my $begin = ${$h_domain}{"$organism:$name"}{$domain}{begin};
	    my $end = ${$h_domain}{"$organism:$name"}{$domain}{end};
	    if ($domain eq "0") {
		return "$organism:$name $begin $end";
	    } else {
		return "$organism:$name($domain) $begin $end";
	    }
	} else {
	    print STDERR "WARNING: no domain information";
	    return "$organism:$name $domain";
	}
    } elsif ($gene =~ /^(\S+?)_(\S+)$/) {
	my ($organism, $name) = ($1, $2);
	# print STDERR "WARNING: no domain information";
	return "$organism:$name";
    } else {
	die;
    }
}

sub tree_node {
    my ($node) = @_;
    
    my $string = $node->to_string;
    my ($gene) = split(":",  $string);
    if ($gene =~ /^(\d+)_(\S+?)_(\S+)$/) {
	my ($domain, $organism, $name) = ($1, $2, $3);
	return "$organism:$name $domain";
    } elsif ($gene =~ /^(\S+?)_(\S+)$/) {
	my ($organism, $name) = ($1, $2);
	return "$organism:$name";
    } else {
	die;
    }
}

sub reset_bar {
    my ($bar_reset, @bar) = @_;
    
    my @new_bar = ();
    for my $bar (@bar) {
	if ($bar != $bar_reset) {
	    push @new_bar, $bar;
	}
    }

    return @new_bar;
}

sub tree_padding {
    my ($offset, @bar) = @_;

    my @padding = ();
    for (my $i=0; $i<$offset; $i++) {
	$padding[$i] = " ";
    }

    for my $bar (@bar) {
	my $pos = abs($bar);
	if ($pos < $offset) {
	    $padding[$pos] = "|";
	}
    }

    return join("", @padding);
}

sub check_duplication {
    my ($r_gene_species, $node) = @_;

    my ($sub_tree1_node, $sub_tree2_node) = $node->each_Descendent;
    
    my @leaves1 = get_sub_tree_leaves($sub_tree1_node);
    my @leaves2 = get_sub_tree_leaves($sub_tree2_node);
    my @species1 = get_species_from_leaves($r_gene_species, @leaves1);
    my @species2 = get_species_from_leaves($r_gene_species, @leaves2);
    my @all_species = uniq(@species1, @species2);
    my @common_species = check_redundancy(@species1, @species2);

    my $sp_disappeared = @all_species - @common_species;

    if (@common_species) {
	return (1, $sp_disappeared);
    } else {
	return 0;
    }
}

sub put_node_labels {
    my ($r_tree) = @_;

    my @node = grep {! $_->is_Leaf} ${$r_tree}->get_root_node->get_all_Descendents;
    for my $node (@node) {
	my $id = $node->internal_id;
	$node->id($id);
    }
}

sub escape_fasta_headers {
    my ($tmp_fasta, %opt) = @_;

    my @line = `cat $tmp_fasta`;
    chomp(@line);
    for (my $i=0; $i<@line; $i++) {
	if ($line[$i] =~ /^>\S+/) {
	    if ($opt{format} and $opt{format} eq "domain_no_with_0") {
		if ($line[$i] =~ /^>(\S+)\((\S+)\)/) {
		    $line[$i] = ">${2}_${1}";
		} elsif ($line[$i] =~ /^>(\S+)/) {
		    $line[$i] = ">0_$1";
		}
	    } elsif ($opt{format} and $opt{format} eq "output_for_newick") {
		if ($line[$i] =~ /^>(\S+)\((\d+)\)/) {
		    $line[$i] = ">${1}__${2}";
		} elsif ($line[$i] =~ /^>(\S+)/) {
		    $line[$i] = ">$1";
		}
	    } else {
		if ($line[$i] =~ /^>(\S+)\((\S+)\)/) {
		    $line[$i] = ">${2}_${1}";
		} elsif ($line[$i] =~ /^>(\S+)/) {
		    $line[$i] = ">$1";
		}
	    }
	    $line[$i] =~ s/:/_/g;
	    $line[$i] =~ s/\|/_/g;
	}
    }

    open(TMP_FASTA, ">$tmp_fasta") || die;
    print TMP_FASTA join("\n", @line), "\n";
    close(TMP_FASTA);
}

1;
