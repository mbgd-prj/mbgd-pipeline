package DomRefine::Score;
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(score_dclst check_score_cache calculate_scores calculate_total_score
	     score_for_lines score_dclst_one_alignment calculate_sp_score pair_score calculate_sp_score_of_clusters
             score_c score_one_c dsp_score_dclst
             count_open_gaps_local_light calculate_scores_local_sum_light calculate_scores_local calculate_scores_local_i
	     set_extention_gap_penalty set_terminal_extention_gap_penalty set_open_gap_penalty set_terminal_open_gap_penalty
             scores_dclst_of_a_cluster stats_dclst_of_a_cluster
	     score_between_two_domains only_stats_between_two_domains
	     score_within_two_domains only_stats_within_two_domains
	     );

use strict;
use Digest::MD5 qw(md5_hex);
use DomRefine::Read;
use DomRefine::Align;
use DomRefine::General;

my $USE_SCORE_CACHE = 0;
my $IGNORE_SAME_SPECIES = 0;
my $DOMAIN_EXISTENCE_THRESHOLD = 0;
my $IGNORE_TERMINAL_GAPS = 0;

# creating tmp files here => OK ??

my $TMP_ALIGNMENT_TO_GET_SCORE = define_tmp_file("alignment_to_get_score");
my $TMP_DCLST_TO_GET_SCORE_1 = define_tmp_file("dclst_to_get_score1");
my $TMP_DCLST_TO_GET_SCORE_2 = define_tmp_file("dclst_to_get_score2");
my $TMP_DCLST_TO_GET_SCORE_3 = define_tmp_file("dclst_to_get_score3");
my $TMP_DCLST_TO_GET_SCORE_4 = define_tmp_file("dclst_to_get_score4");
my $TMP_DCLST_TO_GET_SCORE_5 = define_tmp_file("dclst_to_get_score5");
my $TMP_DCLST_TO_GET_SCORE_6 = define_tmp_file("dclst_to_get_score6");
my $TMP_DCLST_TO_GET_SCORE_8 = define_tmp_file("dclst_to_get_score8");
END {
    remove_tmp_file($TMP_ALIGNMENT_TO_GET_SCORE);
    remove_tmp_file("$TMP_ALIGNMENT_TO_GET_SCORE.err");
    remove_tmp_file($TMP_DCLST_TO_GET_SCORE_1);
    remove_tmp_file($TMP_DCLST_TO_GET_SCORE_2);
    remove_tmp_file($TMP_DCLST_TO_GET_SCORE_3);
    remove_tmp_file($TMP_DCLST_TO_GET_SCORE_4);
    remove_tmp_file($TMP_DCLST_TO_GET_SCORE_5);
    remove_tmp_file($TMP_DCLST_TO_GET_SCORE_6);
    remove_tmp_file($TMP_DCLST_TO_GET_SCORE_8);
}

################################################################################
### Parameters #################################################################
################################################################################

my $OPEN_GAP_PENALTY = 10;
my $TERMINAL_OPEN_GAP_PENALTY = 10;
my $EXTENTION_GAP_PENALTY = 0.5;
my $TERMINAL_EXTENTION_GAP_PENALTY = 0.5;

my %M = ();

sub set_open_gap_penalty {
    my ($value) = @_;
    $OPEN_GAP_PENALTY = $value;
}

sub set_extention_gap_penalty {
    my ($value) = @_;
    $EXTENTION_GAP_PENALTY = $value;
}

sub set_terminal_open_gap_penalty {
    my ($value) = @_;
    $TERMINAL_OPEN_GAP_PENALTY = $value;
}

sub set_terminal_extention_gap_penalty {
    my ($value) = @_;
    $TERMINAL_EXTENTION_GAP_PENALTY = $value;
}

sub read_score_matrix {
    my ($file) = @_;

    my $lines = `cat $file | grep -v '#'`;
    my @line = split("\n", $lines);
    my @a = split(" ", $line[0]);
    my @b = ();
    for (my $i=1; $i<@line; $i++) {
	my ($b, @x) = split(" ", $line[$i]);
	for (my $j=0; $j<@x; $j++) {
	    $M{$b}{$a[$j]} = $x[$j];
	}
    }
}

################################################################################
### C version of functions #####################################################
################################################################################

sub score_c {
    my ($dclst_file, %opt) = @_;

    my $score_sum = 0;
    my $n_aa = 0;
    for my $cluster (get_clusters($dclst_file)) {
	my $dclst = extract_dclst($dclst_file, $cluster);
	my @score = score_one_c($dclst, %opt);
	$score_sum += $score[0];
	$n_aa += $score[2];
    }

    return ($score_sum, $n_aa);
}

sub score_one_c {
    my ($dclst, %opt) = @_;

    save_contents($dclst, $TMP_DCLST_TO_GET_SCORE_6);

    my @gene = ();
    create_alignment($TMP_DCLST_TO_GET_SCORE_6, \@gene, $TMP_ALIGNMENT_TO_GET_SCORE, %opt);

    return dsp_score($TMP_DCLST_TO_GET_SCORE_6, $TMP_ALIGNMENT_TO_GET_SCORE, %opt);
}

sub dsp_score_dclst {
    my ($dclst, $tmp_alignment, %opt) = @_;

    save_contents($dclst, $TMP_DCLST_TO_GET_SCORE_6);

    return dsp_score($TMP_DCLST_TO_GET_SCORE_6, $tmp_alignment, %opt);
}

sub dsp_score {
    my ($tmp_dclst, $tmp_alignment, %opt) = @_;

    # my $program = "dsp_score.c";
    my $program = "dsp_score-1.1.0";

    my $option = "";
    if (defined $opt{lines}) {
	$option .= " -l";
    }
    if (defined $opt{between}) {
	$option .= " -b";
    }
    if (defined $opt{from}){
	$option .= " -F $opt{from}";
    }
    if (defined $opt{to}){
	$option .= " -T $opt{to}";
    }
    if (defined $opt{i}) {
	$option .= " -i $opt{i}";
    }
    my @clusters = get_clusters($tmp_dclst);
    if ($opt{clusters}) {
	my $r_clusters = $opt{clusters};
    	@clusters = @{$r_clusters};
    }

    if (@clusters > 256) {
	die; # dsp_score.c can not treat.
    }

    my $command = "$program -m $ENV{DOMREFINE_DIR}/lib/BLOSUM45 -t $tmp_dclst -a $tmp_alignment -e $EXTENTION_GAP_PENALTY -o $OPEN_GAP_PENALTY $option @clusters";
    # print STDERR "command: $command\n";

    print STDERR "$program: ";
    my $start_time = time;
    my $score = `$command`;
    chomp($score);
    my $end_time = time;
    printf STDERR "%d sec\n", $end_time - $start_time;

    return split(/\s/, $score);
}

################################################################################
### Functions ##################################################################
################################################################################

sub score_for_lines {
    my ($dclst_file, $a_gene, %opt) = @_;

    my ($a_a, $a_b) = ([], []);
    my $h_p;
    %{$h_p} = ();
    my ($n, $m) = get_alignment_structure($dclst_file, $a_gene, $a_a, $a_b, $h_p, %opt);
    print STDERR "  matrix: n=$n m=$m\n";

    my %cluster = ();
    my %domain = ();
    get_dclst_structure($dclst_file, \%cluster, \%domain);
    my $a_d = [];
    map_domains($a_gene, \%domain, $a_b, $a_d);

    my @score_for_line = ();
    calculate_scores($a_a, $a_b, $a_d, $h_p, \@score_for_line);
    return @score_for_line;
}

sub calculate_scores { # various scores, # for just one domain,  # cannot treat multi-clusters
    my ($r_a, $r_b, $r_d, $r_p, $r_score_for_line) = @_;
    
    my $N = @{$r_a};
    my $M = @{${$r_a}[0]};

    print STDERR "  domain overlap score..\n";
    my @domain_existence_ratio = ();
    my @domain_exists = ();
    my $n_amino_acids_in_domain = 0;
    my @n_amino_acids_in_domain = ();
    for (my $j=0; $j<$M; $j++) {
	my $domain_existence_count = 0;
	for (my $i=0; $i<$N; $i++) {
	    if (${$r_d}[$i][$j] and ${$r_b}[$i][$j] == 1) {
		$n_amino_acids_in_domain[$i] ++;
		$n_amino_acids_in_domain ++;
		$domain_existence_count ++;
	    }
	}
	$domain_existence_ratio[$j] = $domain_existence_count / $N;
	if ($domain_existence_ratio[$j] > $DOMAIN_EXISTENCE_THRESHOLD) {
	    $domain_exists[$j] = 1;
	} else {
	    $domain_exists[$j] = 0;
	}
    }

    my $domain_overlap_score = sum_sq(@domain_existence_ratio) / sum(@domain_existence_ratio);
    
    print STDERR "  sum of conservation rates score..\n";
    my $sum_of_cons_raw = 0;
    my $sum_of_cons_potential = 0;
    my $sum_of_cons_norm = 0;
    for (my $j=0; $j<$M; $j++) {
	my $f = 0;
	for (my $i=0; $i<$N; $i++) {
	    my $a = ${$r_a}[$i][$j];
	    if (${$r_d}[$i][$j]) {
		$sum_of_cons_raw += ${$r_p}{$j}{$a} * $domain_exists[$j];
	    }
	    $sum_of_cons_potential += ${$r_p}{$j}{$a} * $domain_exists[$j];
	    if (${$r_b}[$i][$j] == 1) {
		$f ++;
	    }
	}
	$sum_of_cons_norm += $f ** 2 * $domain_exists[$j];
    }

    print STDERR "  sum of pairs score";
    my $score = 0;
    for (my $i=0; $i<$N; $i++) {
	print STDERR ".";
	my $sum_of_pairs = 0;
	for (my $i2=0; $i2<$N; $i2++) {
	    if ($i != $i2) {
		for (my $j=0; $j<$M; $j++) {
		    $sum_of_pairs += pair_score(${$r_a}[$i][$j], ${$r_a}[$i2][$j], ${$r_b}[$i][$j], ${$r_b}[$i2][$j], ${$r_d}[$i][$j], ${$r_d}[$i2][$j]); # cannot treat multi-clusters
		}
	    }
	}
	my $gap_penalty = count_open_gaps_for_line($r_b, $r_d, $i) * $OPEN_GAP_PENALTY; # cannot treat multi-clusters
	if ($r_score_for_line) {
	    my $score_for_line = ($sum_of_pairs/2);
	    # $score_for_line -= $gap_penalty;
	    $score_for_line /= $N * $n_amino_acids_in_domain[$i];
	    push @{$r_score_for_line}, $score_for_line;
	}
	$score += $sum_of_pairs/2 - $gap_penalty;
    }

    print STDERR "\n";

    return ($N,
	    $n_amino_acids_in_domain,
	    max(@domain_existence_ratio), 
	    $domain_overlap_score,
	    $sum_of_cons_raw/$sum_of_cons_potential,
	    $score/$N/$n_amino_acids_in_domain,
	    $score
	    );
}

sub calculate_sp_score_of_clusters {
    my ($r_a_a, $r_a_b, $r_a_d) = @_;
    
    my $score_sum = 0;
    for (my $c=0; $c<@{$r_a_a}; $c++) {
	my $score = calculate_sp_score(${$r_a_a}[$c], ${$r_a_b}[$c], ${$r_a_d}[$c]);
	$score_sum += $score;
    }

    return $score_sum;
}

sub calculate_sp_score { # cannot treat multi-clusters
    my ($r_a, $r_b, $r_d) = @_;
    
    my $N = @{$r_a};
    my $M = @{${$r_a}[0]};

    my $sum_of_pairs = 0;
    my $n_open_gap = 0;
    for (my $i=0; $i<$N; $i++) {
	for (my $i2=$i+1; $i2<$N; $i2++) {
	    for (my $j=0; $j<$M; $j++) {
		$sum_of_pairs += pair_score(${$r_a}[$i][$j], ${$r_a}[$i2][$j], ${$r_b}[$i][$j], ${$r_b}[$i2][$j], ${$r_d}[$i][$j], ${$r_d}[$i2][$j]); # cannot treat multi-clusters
	    }
	}
	$n_open_gap += count_open_gaps_for_line($r_b, $r_d, $i); # cannot treat multi-clusters
    }

    my $score = $sum_of_pairs - $n_open_gap * $OPEN_GAP_PENALTY;
    print STDERR "$score (= $sum_of_pairs - $n_open_gap x $OPEN_GAP_PENALTY)\n";
    return $score;
}

sub calculate_sp_score_with_sp { # cannot treat multi-clusters
    my ($r_gene, $r_a, $r_b, $r_d) = @_;
    
    my $N = @{$r_a};
    my $M = @{${$r_a}[0]};

    my $sum_of_pairs = 0;
    my $n_open_gap = 0;
    for (my $i=0; $i<$N; $i++) {
	for (my $i2=$i+1; $i2<$N; $i2++) {
	    for (my $j=0; $j<$M; $j++) {
		$sum_of_pairs += pair_score_with_sp(${$r_gene}[$i], ${$r_gene}[$i2], ${$r_a}[$i][$j], ${$r_a}[$i2][$j], ${$r_b}[$i][$j], ${$r_b}[$i2][$j], ${$r_d}[$i][$j], ${$r_d}[$i2][$j]);
	    }
	}
	$n_open_gap += count_open_gaps_for_line_with_sp($r_gene, $r_b, $r_d, $i); # cannot treat multi-clusters
    }

    my $score = $sum_of_pairs - $n_open_gap * $OPEN_GAP_PENALTY;
    print STDERR "$score (= $sum_of_pairs - $n_open_gap x $OPEN_GAP_PENALTY)\n";
    return $score;
}

sub pair_score_with_sp {
    my ($gene1, $gene2, $a1, $a2, $b1, $b2, $d1, $d2, $cluster1, $cluster2) = @_;

    my ($sp1) = decompose_gene_id($gene1);
    my ($sp2) = decompose_gene_id($gene2);
    if ($sp1 ne $sp2) {
	return pair_score($a1, $a2, $b1, $b2, $d1, $d2, $cluster1, $cluster2);
    } else {
	return 0;
    }
}

sub pair_score {
    my ($a1, $a2, $b1, $b2, $d1, $d2, $cluster1, $cluster2) = @_;
    
    if (! %M) {
	read_score_matrix("$ENV{DOMREFINE_DIR}/lib/BLOSUM45");
    }

    if ($b1==1 && domain_belongs_to_cluster($d1, $cluster1) and
	$b2==1 && domain_belongs_to_cluster($d2, $cluster2)) {
	if (defined $M{$a1}{$a2}){
	    return $M{$a1}{$a2};
	} else {
	    print STDERR "Warning: no score for $a1 and $a2, so return 0\n"; # BUG ? Is it appropriate to return 0 ?
	    return 0;
	}
    } elsif ($b1==1 && domain_belongs_to_cluster($d1, $cluster1) or 
	     $b2==1 && domain_belongs_to_cluster($d2, $cluster2)) {
	if ($b1== -1 or $b2== -1) {
	    return - $TERMINAL_EXTENTION_GAP_PENALTY;
	} else {
	    return - $EXTENTION_GAP_PENALTY;
	}
    } else {
	return 0;
    }
}

sub domain_belongs_to_cluster {
    my ($d, $cluster) = @_;

    if (! defined $cluster) {
	if ($d) {
	    return 1;
	} else {
	    return 0;
	}
    } else {
	if ($d eq $cluster) {
	    return 1;
	} else {
	    return 0;
	}
    }
}

sub calculate_score_local { # BUG. score is not consistent.
    my ($r_a, $r_b, $r_d, $from_j, $to_j, $head_cluster, $tail_cluster, $r_gene, $r_cluster) = @_;

    my $sum_of_pairs = 0;
    my $n_open_gap = 0;

    my $n = @{$r_b};
    my $m = @{${$r_b}[0]};
    for (my $i=0; $i<$n; $i++) {
	for (my $i2=0; $i2<$n; $i2++) {
	    if ($i2 != $i) {
# 		for (my $j=$from_j; $j<=$to_j; $j++) {
		for (my $j=0; $j<$m; $j++) {
# 		    $sum_of_pairs += pair_score(${$r_a}[$i][$j], ${$r_a}[$i2][$j], ${$r_b}[$i][$j], ${$r_b}[$i2][$j], ${$r_d}[$i][$j], ${$r_d}[$i2][$j]); # cannot treat multi-clusters
		    $sum_of_pairs += pair_score(${$r_a}[$i][$j], ${$r_a}[$i2][$j], ${$r_b}[$i][$j], ${$r_b}[$i2][$j], ${$r_d}[$i][$j], ${$r_d}[$i2][$j], $head_cluster, $head_cluster);
		    $sum_of_pairs += pair_score(${$r_a}[$i][$j], ${$r_a}[$i2][$j], ${$r_b}[$i][$j], ${$r_b}[$i2][$j], ${$r_d}[$i][$j], ${$r_d}[$i2][$j], $tail_cluster, $tail_cluster);
		}
	    }
	}
# 	$n_open_gap += count_open_gaps_local_light($r_b, $r_d, $from_j, $to_j, $i); # cannot treat multi-clusters
    }
    $n_open_gap += count_open_gaps($r_gene, $r_b, $r_d, $head_cluster, $r_cluster); # this is slow
    $n_open_gap += count_open_gaps($r_gene, $r_b, $r_d, $tail_cluster, $r_cluster); # this is slow

    my $score_local_sum = $sum_of_pairs / 2 - ($n_open_gap) * $OPEN_GAP_PENALTY;

    return $score_local_sum;
}

sub calculate_scores_local_sum_light {
    my ($r_a_a, $r_a_b, $r_a_d, $gene, $r_h_gene_idx, $from_pos, $to_pos, $r_a_get_j, @cluster_idx) = @_;

    my $sum_of_pairs = 0;
    my $n_open_gap = 0;

    for my $c (@cluster_idx) {
	my $i = ${${$r_h_gene_idx}[$c]}{$gene};
	my $from_j = ${${$r_a_get_j}[$c][$i]}{$from_pos};
	my $to_j = ${${$r_a_get_j}[$c][$i]}{$to_pos};
	my $r_a = ${$r_a_a}[$c];
	my $r_b = ${$r_a_b}[$c];
	my $r_d = ${$r_a_d}[$c];
	my $n = @{$r_b};
	for (my $i2=0; $i2<$n; $i2++) {
	    if ($i2 != $i) {
		for (my $j=$from_j; $j<=$to_j; $j++) {
		    $sum_of_pairs += pair_score(${$r_a}[$i][$j], ${$r_a}[$i2][$j], ${$r_b}[$i][$j], ${$r_b}[$i2][$j], ${$r_d}[$i][$j], ${$r_d}[$i2][$j]); # cannot treat multi-clusters
		}
	    }
	}
	$n_open_gap += count_open_gaps_local_light($r_b, $r_d, $from_j, $to_j, $i); # cannot treat multi-clusters
    }
    my $score_local_sum = $sum_of_pairs - ($n_open_gap) * $OPEN_GAP_PENALTY;

    return $score_local_sum;
}

sub calculate_scores_local {
    my ($r_a, $r_b, $r_d, $r_cluster, $gene, $r_gene, $i, $from_pos, $to_pos, $r_get_j) = @_;

    my ($from_j, $to_j) = (${${$r_get_j}[$i]}{$from_pos}, ${${$r_get_j}[$i]}{$to_pos});

    # my $score = 0;
    my $sum_of_pairs = 0;
    my $n_open_gap = 0;
    for my $cluster (keys %{$r_cluster}) {
	# my $sum_of_pairs = 0;
	unless (${$r_cluster}{$cluster}{$gene}) {
	    next;
	}
	for (my $j=$from_j; $j<=$to_j; $j++) {
	    for (my $i2=0; $i2<@{$r_b}; $i2++) {
		unless (${$r_cluster}{$cluster}{${$r_gene}[$i2]}) {
		    next;
		}
		if ($i2 != $i) {
		    my $pair_score = pair_score(${$r_a}[$i][$j], ${$r_a}[$i2][$j], ${$r_b}[$i][$j], ${$r_b}[$i2][$j], ${$r_d}[$i][$j], ${$r_d}[$i2][$j], $cluster, $cluster);
		    $sum_of_pairs += $pair_score;
		}
	    }
	}
	# my $n_open_gap = count_open_gaps_local($r_b, $r_d, $cluster, $from_j, $to_j, $i);
	# $score += $sum_of_pairs - $n_open_gap * $OPEN_GAP_PENALTY;
	$n_open_gap += count_open_gaps_local($r_b, $r_d, $cluster, $from_j, $to_j, $i);
    }

    # return $score;
    return ($sum_of_pairs - $n_open_gap * $OPEN_GAP_PENALTY, $sum_of_pairs, $n_open_gap);
}

sub calculate_scores_local_i {
    my ($r_a, $r_b, $r_d, $r_cluster, $gene, $r_gene, $i) = @_;

    my $m = @{${$r_a}[0]};

    my $score = 0;
    for my $cluster (keys %{$r_cluster}) {
	my $sum_of_pairs = 0;
	unless (${$r_cluster}{$cluster}{$gene}) {
	    next;
	}
	for (my $j=0; $j<$m; $j++) {
	    for (my $i2=0; $i2<@{$r_b}; $i2++) {
		unless (${$r_cluster}{$cluster}{${$r_gene}[$i2]}) {
		    next;
		}
		if ($i2 != $i) {
		    my $pair_score = pair_score(${$r_a}[$i][$j], ${$r_a}[$i2][$j], ${$r_b}[$i][$j], ${$r_b}[$i2][$j], ${$r_d}[$i][$j], ${$r_d}[$i2][$j], $cluster, $cluster);
		    $sum_of_pairs += $pair_score;
		}
	    }
	}
	my $n_open_gap = count_open_gaps_local($r_b, $r_d, $cluster, 0, $m-1, $i);
	$score += $sum_of_pairs - $n_open_gap * $OPEN_GAP_PENALTY;
    }

    return $score;
}

################################################################################
### Gap ########################################################################
################################################################################

sub count_open_gaps_for_line { # cannot treat multi-clusters
    my ($r_b, $r_d, $i) = @_;

    my $N = @{$r_b};

    my $n_open_gap = 0;
    my @gap_start = ();
    my @gap_end = ();
    extract_gaps(${$r_b}[$i], ${$r_d}[$i], \@gap_start, \@gap_end, 0);
    for (my $i2=0; $i2<$N; $i2++) {
	for (my $k=0; $k<@gap_end; $k++) {
	    if (aa_in_other_seq(${$r_b}[$i2], ${$r_d}[$i2], $gap_start[$k], $gap_end[$k])) { # cannot treat multi-clusters
		$n_open_gap ++;
	    }
	}
    }

    return $n_open_gap;
}

sub count_open_gaps_for_line_with_sp { # cannot treat multi-clusters
    my ($r_gene, $r_b, $r_d, $i) = @_;

    my $N = @{$r_b};

    my $n_open_gap = 0;
    my @gap_start = ();
    my @gap_end = ();
    extract_gaps(${$r_b}[$i], ${$r_d}[$i], \@gap_start, \@gap_end, 0);
    for (my $i2=0; $i2<$N; $i2++) {
	my ($sp1) = decompose_gene_id(${$r_gene}[$i]);
	my ($sp2) = decompose_gene_id(${$r_gene}[$i2]);
	if ($sp1 eq $sp2) {
	    next;
	}
	for (my $k=0; $k<@gap_end; $k++) {
	    if (aa_in_other_seq(${$r_b}[$i2], ${$r_d}[$i2], $gap_start[$k], $gap_end[$k])) { # cannot treat multi-clusters
		$n_open_gap ++;
	    }
	}
    }

    return $n_open_gap;
}

sub count_open_gaps {
    my ($r_gene, $r_b, $r_d, $cluster, $r_h_cluster) = @_;

    my $N = @{$r_b};

    my $n_open_gap = 0;
    for (my $i=0; $i<$N; $i++) {
	unless (${$r_h_cluster}{$cluster}{${$r_gene}[$i]}) {
	    next;
	}
	my @gap_start = ();
	my @gap_end = ();
	extract_gaps(${$r_b}[$i], ${$r_d}[$i], \@gap_start, \@gap_end, $cluster);
	for (my $i2=0; $i2<$N; $i2++) {
	    unless (${$r_h_cluster}{$cluster}{${$r_gene}[$i2]}) {
		next;
	    }
	    if ($IGNORE_SAME_SPECIES){
		my ($sp1) = decompose_gene_id(${$r_gene}[$i]);
		my ($sp2) = decompose_gene_id(${$r_gene}[$i2]);
		if ($sp1 eq $sp2) {
		    next;
		}
	    }
	    for (my $k=0; $k<@gap_end; $k++) {
		if (aa_in_other_seq(${$r_b}[$i2], ${$r_d}[$i2], $gap_start[$k], $gap_end[$k], $cluster)) {
		    $n_open_gap ++;
		}
	    }
	}
    }

    return $n_open_gap;
}

sub aa_in_other_seq {
    my ($r_b_line, $r_d_line, $gap_start, $gap_end, $cluster) = @_;
    
    for (my $j=$gap_start; $j<=$gap_end; $j++) {
	if (${$r_b_line}[$j]) {
	    if (defined $cluster) {
		if (${$r_d_line}[$j] eq $cluster) {
		    return 1;
		}
	    } else {
		if (${$r_d_line}[$j]) {
		    return 1;
		}
	    }
	}
    }
    
    return 0;
}

sub extract_gaps {
    my ($r_b_line, $r_d_line, $r_gap_start, $r_gap_end, $cluster, %opt) = @_;

    my $M = @{$r_b_line};

    my $inside_gap = 0;
    for (my $j=0; $j<$M; $j++) {
	if (domain_exists(${$r_d_line}[$j], $cluster)) {
	    # gap starts
	    if (${$r_b_line}[$j] != 1 and ! $inside_gap) {
		$inside_gap = 1;
		push @{$r_gap_start}, $j;
	    }
	    # gap ends
	    if (${$r_b_line}[$j] == 1 and $inside_gap) {
		$inside_gap = 0;
		push @{$r_gap_end}, $j - 1;
	    }
	}
	# first or last gap starts
	if (! domain_exists(${$r_d_line}[$j], $cluster) and ! $inside_gap) {
	    $inside_gap = 1;
	    push @{$r_gap_start}, $j;
	}
    }
    # last gap ends
    if ($inside_gap) {
	push @{$r_gap_end}, $M - 1;
    }

    if ($opt{ignore_terminal_gaps}) {
	if (${$r_gap_start}[0] and ${$r_gap_start}[0] == 0) {
	    shift @{$r_gap_start};
	    shift @{$r_gap_end};
	}
	if (${$r_gap_end}[-1] and ${$r_gap_end}[-1] == $M-1) {
	    pop @{$r_gap_start};
	    pop @{$r_gap_end};
	}
    }

    if (@{$r_gap_start} != @{$r_gap_end}) {
	die;
    }
}

sub domain_exists {
    my ($d, $cluster)  = @_;
    if ($cluster) {
	# d is clusterId or other clsuterId
	if ($d eq $cluster) {
	    return 1;
	} else {
# 	    return 2;
	    return 0;
	}
    } else {
	# d is 0 or 1
	if ($d) {
	    return 1;
	} else {
	    return 0;
	}
    }
}

sub count_open_gaps_local_light { # cannot treat multi-clusters
    my ($r_b, $r_d, $from_j, $to_j, $i) = @_;

    my $n_open_gap = 0;

    my @gap_start = ();
    my @gap_end = ();
    extract_gaps_local(${$r_b}[$i], ${$r_d}[$i], $from_j, $to_j, \@gap_start, \@gap_end);
    my $n = @{$r_b};
    for (my $i2=0; $i2<$n; $i2++) {
	if ($i != $i2) {
	    for (my $k=0; $k<@gap_end; $k++) {
		if (aa_in_other_seq(${$r_b}[$i2], ${$r_d}[$i2], $gap_start[$k], $gap_end[$k])) { # cannot treat multi-clusters
		    $n_open_gap ++;
		}
	    }
	    my @gap_start2 = ();
	    my @gap_end2 = ();
	    extract_gaps_local(${$r_b}[$i2], ${$r_d}[$i2], $from_j, $to_j, \@gap_start2, \@gap_end2);
	    for (my $k2=0; $k2<@gap_end2; $k2++) {
		if (aa_in_other_seq(${$r_b}[$i], ${$r_d}[$i], $gap_start2[$k2], $gap_end2[$k2])) { # cannot treat multi-clusters
		    $n_open_gap ++;
		}
	    }
	}
    }

    return $n_open_gap;
}

sub count_open_gaps_local {
    my ($r_b, $r_d, $cluster, $from_j, $to_j, $i) = @_;

    my $n_open_gap = 0;

    my @gap_start = ();
    my @gap_end = ();
    extract_gaps_of_cluster_local(${$r_b}[$i], ${$r_d}[$i], $cluster, $from_j, $to_j, \@gap_start, \@gap_end);
    for (my $i2=0; $i2<@{$r_b}; $i2++) {
	if ($i != $i2) {
	    for (my $k=0; $k<@gap_end; $k++) {
		if (aa_in_other_seq(${$r_b}[$i2], ${$r_d}[$i2], $gap_start[$k], $gap_end[$k], $cluster)) {
		    $n_open_gap ++;
		}
	    }
	}
    }
    for (my $i2=0; $i2<@{$r_b}; $i2++) {
	if ($i != $i2) {
	    @gap_start = ();
	    @gap_end = ();
	    extract_gaps_of_cluster_local(${$r_b}[$i2], ${$r_d}[$i2], $cluster, $from_j, $to_j, \@gap_start, \@gap_end);
	    for (my $k=0; $k<@gap_end; $k++) {
		if (aa_in_other_seq(${$r_b}[$i], ${$r_d}[$i], $gap_start[$k], $gap_end[$k], $cluster)) {
		    $n_open_gap ++;
		}
	    }
	}
    }

    return $n_open_gap;
}

### some bug ??
sub extract_gaps_local {
    my ($r_b_line, $r_d_line, $from_j, $to_j, $r_gap_start, $r_gap_end, %opt) = @_;

    $from_j = max(0, $from_j-1);

    my $M = @{$r_b_line};

    my $inside_gap = 0;
    # first gap starts
    if (! ${$r_d_line}[$from_j]) {
	$inside_gap = 1;
	push @{$r_gap_start}, 0;
    } elsif (${$r_b_line}[$from_j] != 1) {
	for (my $j=$from_j; $j>=0; $j--) {
	    if (${$r_b_line}[$from_j] == 1) {
		$inside_gap = 1;
		push @{$r_gap_start}, $j+1;
		last;
	    }
	}
    }
    for (my $j=$from_j+1; $j<=$to_j; $j++) {
	if (${$r_d_line}[$j]) {
	    # gap starts
	    if (${$r_b_line}[$j] != 1 and ! $inside_gap) {
		$inside_gap = 1;
		push @{$r_gap_start}, $j;
	    }
	    # gap ends
	    if (${$r_b_line}[$j] == 1 and $inside_gap) {
		$inside_gap = 0;
		push @{$r_gap_end}, $j - 1;
	    }
	}
	# last gap starts
	if (! ${$r_d_line}[$j] and ! $inside_gap) {
	    $inside_gap = 1;
	    push @{$r_gap_start}, $j;
	}
    }
    # last gap ends
    if ($inside_gap) {
	if (! ${$r_d_line}[$to_j]) {
	    push @{$r_gap_end}, $M;
	} elsif (${$r_b_line}[$to_j] != 1) {
	    for (my $j=$to_j; $j<$M; $j++) {
		if (${$r_b_line}[$j] == 1) {
		    push @{$r_gap_end}, $j-1;
		    last;
		}
	    }
	}
    }

    if ($opt{ignore_terminal_gaps}) {
	if (${$r_gap_start}[0] and ${$r_gap_start}[0] == 0) {
	    shift @{$r_gap_start};
	    shift @{$r_gap_end};
	}
	if (${$r_gap_end}[-1] and ${$r_gap_end}[-1] == $M-1) {
	    pop @{$r_gap_start};
	    pop @{$r_gap_end};
	}
    }

    if (@{$r_gap_start} != @{$r_gap_end}) {
	die;
    }
}

### some bug ?? (inherited from extract_gaps_local)
sub extract_gaps_of_cluster_local {
    my ($r_b_line, $r_d_line, $cluster, $from_j, $to_j, $r_gap_start, $r_gap_end, %opt) = @_;

    $from_j = max(0, $from_j-1);

    my $M = @{$r_b_line};

    my $inside_gap = 0;
    # first gap starts
    if (${$r_d_line}[$from_j] ne $cluster) {
	$inside_gap = 1;
	push @{$r_gap_start}, 0;
    } elsif (${$r_b_line}[$from_j] != 1) {
	for (my $j=$from_j; $j>=0; $j--) {
	    if (${$r_b_line}[$from_j] == 1) {
		$inside_gap = 1;
		push @{$r_gap_start}, $j+1;
		last;
	    }
	}
    }
    for (my $j=$from_j+1; $j<=$to_j; $j++) {
	if (${$r_d_line}[$j] eq $cluster) {
	    # gap starts
	    if (${$r_b_line}[$j] != 1 and ! $inside_gap) {
		$inside_gap = 1;
		push @{$r_gap_start}, $j;
	    }
	    # gap ends
	    if (${$r_b_line}[$j] == 1 and $inside_gap) {
		$inside_gap = 0;
		push @{$r_gap_end}, $j - 1;
	    }
	}
	# last gap starts
	if (${$r_d_line}[$j] ne $cluster and ! $inside_gap) {
	    $inside_gap = 1;
	    push @{$r_gap_start}, $j;
	}
    }
    # last gap ends
    if ($inside_gap) {
	if (${$r_d_line}[$to_j] ne $cluster) {
	    push @{$r_gap_end}, $M;
	} elsif (${$r_b_line}[$to_j] != 1) {
	    for (my $j=$to_j; $j<$M; $j++) {
		if (${$r_b_line}[$j] == 1) {
		    push @{$r_gap_end}, $j-1;
		    last;
		}
	    }
	}
    }

    if ($opt{ignore_terminal_gaps}) {
	if (${$r_gap_start}[0] and ${$r_gap_start}[0] == 0) {
	    shift @{$r_gap_start};
	    shift @{$r_gap_end};
	}
	if (${$r_gap_end}[-1] and ${$r_gap_end}[-1] == $M-1) {
	    pop @{$r_gap_start};
	    pop @{$r_gap_end};
	}
    }

    if (@{$r_gap_start} != @{$r_gap_end}) {
	die;
    }
}

################################################################################
### Between clusters ###########################################################
################################################################################

sub only_stats_between_two_domains {
    my ($r_gene, $r_a, $r_b, $r_d, $r_h_cluster, $cluster1, $cluster2) = @_;
    
    my $n = @{$r_a};
    my $m = @{${$r_a}[0]};

    my $total_n_pair = 0;
    for (my $i=0; $i<$n; $i++) {
	unless (${$r_h_cluster}{$cluster1}{${$r_gene}[$i]}) {
	    next;
	}
	for (my $i2=0; $i2<$n; $i2++) {
	    if ($i == $i2) {
		next;
	    }
	    unless (${$r_h_cluster}{$cluster2}{${$r_gene}[$i2]}) {
		next;
	    }
	    for (my $j=0; $j<$m; $j++) {
		$total_n_pair += check_pair_between_domains(${$r_d}[$i][$j], ${$r_d}[$i2][$j], $cluster1, $cluster2);
	    }
	}
    }

    return $total_n_pair;
}

sub check_pair_between_domains {
    my ($d1, $d2, $cluster1, $cluster2) = @_;
    
    if ($d1 eq $cluster1 and $d2 eq $cluster2) {
	return 1;
    } else {
	return 0;
    }
}

sub score_within_two_domains {
    my ($r_gene, $r_a, $r_b, $r_d, $r_h_cluster, $cluster1, $cluster2) = @_;
    
    my $n = @{$r_a};
    my $m = @{${$r_a}[0]};

    my $sum_of_pairs = 0;
    for (my $i=0; $i<$n; $i++) {
	for (my $i2=$i+1; $i2<$n; $i2++) {
	    for (my $j=0; $j<$m; $j++) {
		if (check_pair_within_domains(${$r_d}[$i][$j], ${$r_d}[$i2][$j], $cluster1, $cluster2)) {
		    $sum_of_pairs += pair_score(${$r_a}[$i][$j], ${$r_a}[$i2][$j], ${$r_b}[$i][$j], ${$r_b}[$i2][$j], ${$r_d}[$i][$j], ${$r_d}[$i2][$j]);
		}
		# do not include terminal extension gap penalty
	    }
	}
    }

    my $n_open_gap = count_open_gaps($r_gene, $r_b, $r_d, $cluster1, $r_h_cluster) + count_open_gaps($r_gene, $r_b, $r_d, $cluster2, $r_h_cluster);

    my $score = $sum_of_pairs - $n_open_gap * $OPEN_GAP_PENALTY;
    print STDERR "score = $score (= $sum_of_pairs - $n_open_gap x $OPEN_GAP_PENALTY)\n";

    return $score;
}

sub only_stats_within_two_domains {
    my ($r_gene, $r_a, $r_b, $r_d, $r_h_cluster, $cluster1, $cluster2) = @_;
    
    my $n = @{$r_a};
    my $m = @{${$r_a}[0]};

    my $total_n_pair = 0;
    for (my $i=0; $i<$n; $i++) {
	for (my $i2=$i+1; $i2<$n; $i2++) {
	    for (my $j=0; $j<$m; $j++) {
		$total_n_pair += check_pair_within_domains(${$r_d}[$i][$j], ${$r_d}[$i2][$j], $cluster1, $cluster2);
	    }
	}
    }

    return $total_n_pair;
}

sub check_pair_within_domains {
    my ($d1, $d2, $cluster1, $cluster2) = @_;
    
    if ($d1 eq $cluster1 and $d2 eq $cluster1 or
	$d1 eq $cluster2 and $d2 eq $cluster2) {
	return 1;
    } else {
	return 0;
    }
}

################################################################################
### Slow code written in Perl (Use C version instead) ##########################
################################################################################

# should be replaced by score_dclst_one_alignment
# bug in trim option ?
sub scores_dclst_of_a_cluster {
    my ($dclst, %opt) = @_;

    save_contents($dclst, $TMP_DCLST_TO_GET_SCORE_4);

    my @gene = ();
    my @a =();
    my @b = ();
    my %p = ();
    get_alignment_structure($TMP_DCLST_TO_GET_SCORE_4, \@gene, \@a, \@b, \%p, %opt);
    my %domain = ();
    my %cluster = ();
    get_dclst_structure($TMP_DCLST_TO_GET_SCORE_4, \%cluster, \%domain); # dclst should be trimmed ?
    my @d = (); 
    map_domains(\@gene, \%domain, \@b, \@d);

    my $n = @a;
    my $m = @{$a[0]};
    my $n_amino_acids_in_domain = 0;
    for (my $i=0; $i<$n; $i++) {
	for (my $j=0; $j<$m; $j++) {
	    if ($d[$i][$j] and $b[$i][$j] == 1) {
		$n_amino_acids_in_domain ++;
	    }
	}
    }

    my $score;
    if ($IGNORE_SAME_SPECIES) {
	$score = calculate_sp_score_with_sp(\@gene, \@a, \@b, \@d);
    } else {
	$score = calculate_sp_score(\@a, \@b, \@d);
    }

    return ($n, $m, $n_amino_acids_in_domain, $score, $score/($n*$n_amino_acids_in_domain));
}

sub stats_dclst_of_a_cluster {
    my ($dclst, %opt) = @_;

    save_contents($dclst, $TMP_DCLST_TO_GET_SCORE_3);

    my @gene = ();
    my @a =();
    my @b = ();
    my %p = ();
    get_alignment_structure($TMP_DCLST_TO_GET_SCORE_3, \@gene, \@a, \@b, \%p, %opt);
    my %domain = ();
    my %cluster = ();
    get_dclst_structure($TMP_DCLST_TO_GET_SCORE_3, \%cluster, \%domain);
    my @d = (); 
    map_domains(\@gene, \%domain, \@b, \@d);

    my $n = @a;
    my $m = @{$a[0]};
    my $n_amino_acids_in_domain = 0;
    for (my $i=0; $i<$n; $i++) {
	for (my $j=0; $j<$m; $j++) {
	    if ($d[$i][$j] and $b[$i][$j] == 1) {
		$n_amino_acids_in_domain ++;
	    }
	}
    }

    return ($n, $m, $n_amino_acids_in_domain);
}

sub score_dclst_one_alignment {
    my ($dclst, %opt) = @_;

    save_contents($dclst, $TMP_DCLST_TO_GET_SCORE_2);

    my $score_file = cache_file_path($ENV{DOMREFINE_CACHE}, md5_hex($dclst), ".score");
    my $score_cached = check_score_cache($dclst, $score_file);
    if (defined $score_cached) {
	return $score_cached;
    }

    my @gene = ();
    my @a =();
    my @b = ();
    my %p = ();
    my ($n_seq, $n_col) = get_alignment_structure($TMP_DCLST_TO_GET_SCORE_2, \@gene, \@a, \@b, \%p, %opt);
    my %domain = ();
    my %cluster = ();
    get_dclst_structure($TMP_DCLST_TO_GET_SCORE_2, \%cluster, \%domain);

    my $score = calculate_total_score(\@gene, \%domain, \@a, \@b, \%cluster);
    if ($USE_SCORE_CACHE) {
	save_contents("$score\n$dclst", $score_file);
    }

    return $score;
}

sub check_score_cache {
    my ($dclst_table, $score_file) = @_;

    if ($USE_SCORE_CACHE and -f $score_file) {
	my ($score, @lines) = `cat $score_file`;
	if (join("", @lines) eq $dclst_table) {
	    chomp($score);
	    return $score;
	} else {
	    print STDERR "WARNING: inconsistent score $score_file\n";
	    return;
	}
    } else {
	print STDERR "no score cache $score_file\n";
	return;
    }
}

sub score_dclst {
    my ($dclst, %opt) = @_;

    save_contents($dclst, $TMP_DCLST_TO_GET_SCORE_1);

    my $score_sum = 0;
    for my $cluster (get_clusters($TMP_DCLST_TO_GET_SCORE_1)) {
	$score_sum += score_dclst_one_alignment(extract_dclst($TMP_DCLST_TO_GET_SCORE_1, $cluster), %opt);
    }

    return $score_sum;
}

sub calculate_total_score {
    my ($r_gene, $r_domain, $r_a, $r_b, $r_cluster) = @_;

    my @d = (); 
    print STDERR "  mapping..\n";
    map_domains($r_gene, $r_domain, $r_b, \@d);
    
    my $total_score = 0;
    for my $cluster (keys %{$r_cluster}) { # slow to loop here ?
	print STDERR "  score..\n";
	$total_score += calculate_sp_score_restricted_cluster($r_gene, $r_a, $r_b, \@d, $r_cluster, $cluster);
    }

    return $total_score;
}

sub calculate_sp_score_restricted_cluster {
    my ($r_gene, $r_a, $r_b, $r_d, $r_h_cluster, $cluster) = @_;
    
    my $m = @{${$r_a}[0]};
    
    my $sum_of_pairs = 0;
    for (my $j=0; $j<$m; $j++) {
	$sum_of_pairs += calculate_sp_score_restricted_cluster_column($r_a, $r_b, $r_d, $cluster, $j, $r_h_cluster, $r_gene);
    }

    my $n_open_gap = count_open_gaps($r_gene, $r_b, $r_d, $cluster, $r_h_cluster);
    my $score = $sum_of_pairs - $n_open_gap * $OPEN_GAP_PENALTY;
    print STDERR "[$cluster] score = $score (= $sum_of_pairs - $n_open_gap x $OPEN_GAP_PENALTY)\n";

    return $score;
}

sub calculate_sp_score_restricted_cluster_column {
    my ($r_a, $r_b, $r_d, $cluster, $j, $r_h_cluster, $r_gene) = @_;

    my $n = @{$r_a};

    my $sum_of_pairs = 0;
    for (my $i=0; $i<$n; $i++) {
	unless (${$r_h_cluster}{$cluster}{${$r_gene}[$i]}) {
	    next;
	}
	for (my $i2=$i+1; $i2<$n; $i2++) {
	    unless (${$r_h_cluster}{$cluster}{${$r_gene}[$i2]}) {
		next;
	    }
	    if ($IGNORE_SAME_SPECIES) {
		$sum_of_pairs += pair_score_with_sp(${$r_gene}[$i], ${$r_gene}[$i2], ${$r_a}[$i][$j], ${$r_a}[$i2][$j], ${$r_b}[$i][$j], ${$r_b}[$i2][$j], ${$r_d}[$i][$j], ${$r_d}[$i2][$j], $cluster, $cluster);
	    } else {
		$sum_of_pairs += pair_score(${$r_a}[$i][$j], ${$r_a}[$i2][$j], ${$r_b}[$i][$j], ${$r_b}[$i2][$j], ${$r_d}[$i][$j], ${$r_d}[$i2][$j], $cluster, $cluster);
	    }
	}
    }

    return $sum_of_pairs;
}

sub score_between_two_domains {
    my ($r_gene, $r_a, $r_b, $r_d, $r_h_cluster, $cluster1, $cluster2) = @_;
    
    my $n = @{$r_a};
    my $m = @{${$r_a}[0]};

    my $pair_score_max;
    for (my $i=0; $i<$n; $i++) {
	unless (${$r_h_cluster}{$cluster1}{${$r_gene}[$i]}) {
	    next;
	}
	for (my $i2=0; $i2<$n; $i2++) {
	    if ($i == $i2) {
		next;
	    }
	    unless (${$r_h_cluster}{$cluster2}{${$r_gene}[$i2]}) {
		next;
	    }
	    my $pair_score = 0;
	    for (my $j=0; $j<$m; $j++) {
		$pair_score += pair_score_between(${$r_a}[$i][$j], ${$r_a}[$i2][$j], ${$r_b}[$i][$j], ${$r_b}[$i2][$j], ${$r_d}[$i][$j], ${$r_d}[$i2][$j], $cluster1, $cluster2);
	    }
	    # print "$i,$i2\t$pair_score\t";
	    my $n_open_gap = 0;
	    $n_open_gap += count_open_gaps_between_two_domains($r_gene, $r_b, $r_d, $i, $i2, $cluster1, $cluster2);
	    $n_open_gap += count_open_gaps_between_two_domains($r_gene, $r_b, $r_d, $i2, $i, $cluster2, $cluster1);
	    $pair_score -= $n_open_gap * $OPEN_GAP_PENALTY;
	    if (! defined $pair_score_max or $pair_score > $pair_score_max) {
		$pair_score_max = $pair_score;
	    }
	    # print "$n_open_gap\n";
	    # print "$pair_score\n";
	}
    }

    return ($pair_score_max);
}

sub pair_score_between {
    my ($a1, $a2, $b1, $b2, $d1, $d2, $cluster1, $cluster2) = @_;
    
    if (! %M) {
	read_score_matrix("$ENV{DOMREFINE_DIR}/lib/BLOSUM45");
    }

    if (($b1==1 && $d1 eq $cluster1) && ($b2==1 && $d2 eq $cluster2)) {
	if (defined $M{$a1}{$a2}){
	    return $M{$a1}{$a2};
	} else {
	    print STDERR "Warning: no score for $a1 and $a2, so return 0\n"; # BUG ? Is it appropriate to return 0 ?
	    return 0;
	}

    	# 0: internal gap, -1:terminal gap
    } elsif (($b1==1 && $d1 eq $cluster1)           && ($b2==0 && $d2 ne $cluster2 || $b2==1)
	  or ($b1==0 && $d1 ne $cluster1 || $b1==1) && ($b2==1 && $d2 eq $cluster2)
	) {
	return 0;

    } elsif (($b1==1 && $d1 eq $cluster1) || ($b2==1 && $d2 eq $cluster2)) {
	return - $EXTENTION_GAP_PENALTY;

    } else {
	return 0;
    }
}

sub count_open_gaps_between_two_domains {
    my ($r_gene, $r_b, $r_d, $i, $i2, $cluster1, $cluster2) = @_;

    my $n_open_gap = 0;
    my @gap_start = ();
    my @gap_end = ();

    extract_gaps(${$r_b}[$i], ${$r_d}[$i], \@gap_start, \@gap_end, $cluster1);
    for (my $k=0; $k<@gap_end; $k++) {
	# print STDERR " ($k)$gap_start[$k]-$gap_end[$k]:", $gap_end[$k]-$gap_start[$k]+1;
	if (aa_in_other_seq(${$r_b}[$i2], ${$r_d}[$i2], $gap_start[$k], $gap_end[$k], $cluster2)) {
	    # print STDERR " open";
	    $n_open_gap ++;
	}
    }
    # print STDERR "\n open:$n_open_gap\n";

    return $n_open_gap;
}

1;
