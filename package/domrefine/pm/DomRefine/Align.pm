package DomRefine::Align;
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(aligner create_alignment get_alignment_structure get_alignment_matrices
             incremental_alignment
             read_alignment map_domains create_get_j create_get_pos_from_a
             fasta_to_sorted_tsv
	     );

use strict;
use Digest::MD5 qw(md5_hex);
use DomRefine::General;
use DomRefine::Read;

my $TMP_SEQ_TO_CREATE_ALIGNMENT = define_tmp_file("seq_to_create_alignment");
my $TMP_SEQ_FOR_CLUSTALW = define_tmp_file("seq_for_clustalw");
my $TMP_ALIGNMENT_BY_ALIGNER = define_tmp_file("output_for_clustalw");
my $TMP_ALIGNMENT_TO_GET_ALIGNMENT_STRUCTURE = define_tmp_file("alignment_to_get_alignment_structure");
my $TMP_SEQ_COMMON = define_tmp_file("seq_common");
my $TMP_SEQ_ADD = define_tmp_file("seq_add");
my $TMP_SEQ_UNION = define_tmp_file("union_seq");
my $TMP_ALIGNMENT_COMMON = define_tmp_file("alignment_common");
my $TMP_ALIGNMENT_UNION = define_tmp_file("alignment_union");
END {
    remove_tmp_file($TMP_SEQ_FOR_CLUSTALW);
    remove_tmp_file($TMP_ALIGNMENT_BY_ALIGNER);
    remove_tmp_file($TMP_SEQ_TO_CREATE_ALIGNMENT);
    remove_tmp_file($TMP_ALIGNMENT_TO_GET_ALIGNMENT_STRUCTURE);
    remove_tmp_file("$TMP_ALIGNMENT_TO_GET_ALIGNMENT_STRUCTURE.err");
    remove_tmp_file($TMP_SEQ_COMMON);
    remove_tmp_file($TMP_SEQ_ADD);
    remove_tmp_file($TMP_SEQ_UNION);
    remove_tmp_file($TMP_ALIGNMENT_COMMON);
    remove_tmp_file($TMP_ALIGNMENT_UNION);
}

################################################################################
### Public functions ###########################################################
################################################################################

sub aligner {
    my ($seq_file, $alignment_file, %opt) = @_;

    my $aligner = $ENV{DOMREFINE_ALIGNER} || die "set environmental DOMREFINE_ALIGNER\n";

    # my $cache_file = cache_file_path($ENV{DOMREFINE_CACHE}, md5_hex(`cat $seq_file`), ".$aligner");
    my $cache_file = cache_file_path($ENV{DOMREFINE_CACHE}, md5_hex(fasta_to_sorted_tsv($seq_file)), ".$aligner");
    # sort before cache -> output of the aligner should be orderd ?
    my $cache_file_stepwise_merge = cache_file_path($ENV{DOMREFINE_CACHE_STEPWISE_MERGE}, md5_hex(fasta_to_sorted_tsv($seq_file)), ".$aligner");

    if ($ENV{DOMREFINE_READ_ALIGNMENT_CACHE} and found_alignment_cache($seq_file, $cache_file)) {
        print STDERR " found alignment cache $cache_file\n";
        system "cp $cache_file $alignment_file";
        return $cache_file;
    }
    if ($opt{stepwise_merge} and found_alignment_cache($seq_file, $cache_file_stepwise_merge)) {
        print STDERR " found alignment cache $cache_file_stepwise_merge\n";
        system "cp $cache_file_stepwise_merge $alignment_file";
        return $cache_file_stepwise_merge;
    }

    my ($n_seq, $n_aa, $mean_len) = seq_file_stat($seq_file);
    my $expected_time = $n_seq * log($n_seq) * $mean_len ** 2 * 4 / 10000000;
    printf STDERR "$aligner: n_seq=$n_seq n_aa=$n_aa mean_len=%.1f exp_time=%.1f ", $mean_len, $expected_time;
    my $start_time = time;

    if ($aligner eq "mafft6") {
        system "mafft --auto --anysymbol $seq_file > $alignment_file 2> $alignment_file.err";
        # system "mafft --localpair $seq_file > $alignment_file 2> $alignment_file.err";
    } elsif ($aligner eq "mafft-fast") {
        system "/bio/bin/mafft --retree 1 --maxiterate 0 --anysymbol $seq_file > $alignment_file 2> $alignment_file.err";
    } elsif ($aligner eq "mafft") {
        system "/bio/bin/mafft --anysymbol $seq_file > $alignment_file 2> $alignment_file.err";
    } elsif ($aligner eq "mafft-fftns") {
        system "/bio/bin/fftns --anysymbol $seq_file > $alignment_file 2> $alignment_file.err";
    } elsif ($aligner eq "mafft-auto") {
        system "/bio/bin/mafft --auto --anysymbol $seq_file > $alignment_file 2> $alignment_file.err";
    } elsif ($aligner eq "mafft-fftnsi") {
        system "/bio/bin/fftnsi --anysymbol $seq_file > $alignment_file 2> $alignment_file.err";
    } elsif ($aligner eq "mafft-linsi") {
        system "/bio/bin/linsi --anysymbol $seq_file > $alignment_file 2> $alignment_file.err";
    } elsif ($aligner eq "muscle") {
        system "muscle -in $seq_file -out $alignment_file 2> $alignment_file.err";
        # system "muscle -stable -in $seq_file -out $alignment_file 2> $alignment_file.err";
    } elsif ($aligner eq "clustalo1.1") {
        system "clustalo-1.1.0-linux-64 --auto -i $seq_file -o $alignment_file -l $alignment_file.err --threads=4 --force -v -v";
    } elsif ($aligner eq "clustalo1.2") {
        system "clustalo-1.2.0-Ubuntu-x86_64 --auto -i $seq_file -o $alignment_file -l $alignment_file.err --threads=4 --force -v -v";
    } elsif ($aligner eq "clustalo1.2kimura") {
        system "clustalo-1.2.0-Ubuntu-x86_64 --auto -i $seq_file -o $alignment_file -l $alignment_file.err --threads=4 --force -v -v --use-kimura";
    } elsif ($aligner eq "clustalo1.2order") {
        system "clustalo-1.2.0-Ubuntu-x86_64 --auto -i $seq_file -o $alignment_file -l $alignment_file.err --threads=4 --force -v -v --output-order=tree-order";
    } elsif ($aligner eq "clustalo1.2orderkimura") {
        system "clustalo-1.2.0-Ubuntu-x86_64 --auto -i $seq_file -o $alignment_file -l $alignment_file.err --threads=4 --force -v -v --output-order=tree-order --use-kimura";
    } elsif ($aligner eq "clustalw") {
        convert_seq_file($seq_file, $TMP_SEQ_FOR_CLUSTALW);
        system "clustalw2 -endgaps -infile=$TMP_SEQ_FOR_CLUSTALW -outfile=$TMP_ALIGNMENT_BY_ALIGNER -output=gde -outorder=input 2>&1 > $alignment_file.err";
    } elsif ($aligner eq "t_coffee") {
        system "t_coffee $seq_file -mode regular -outfile=$TMP_ALIGNMENT_BY_ALIGNER -output=fasta_aln > /dev/null 2> $alignment_file.err";
    } elsif ($aligner eq "t_coffee_quick") {
        system "t_coffee $seq_file -mode quickaln -outfile=$alignment_file -output=fasta_aln > /dev/null 2> $alignment_file.err";
    } elsif ($aligner eq "famsa") {
        system "famsa-1.6.1-linux-static $seq_file $alignment_file 2> $alignment_file.err";
    }

    if ($aligner =~ /^(clustalw|t_coffee|t_coffee_quick)$/) {
        convert_alignment_file($TMP_ALIGNMENT_BY_ALIGNER, $alignment_file);
    }

    my $end_time = time;
    printf STDERR "%d sec\n", $end_time - $start_time;

    if ($opt{stepwise_merge}) {
        print STDERR " saving alignment cache $cache_file_stepwise_merge\n";
        print STDERR " saving alignment cache $cache_file_stepwise_merge.err\n";
        system "cp -a $alignment_file $cache_file_stepwise_merge";
        system "cp -a $alignment_file.err $cache_file_stepwise_merge.err";
        return $cache_file_stepwise_merge;
    }
    if ($ENV{DOMREFINE_WRITE_ALIGNMENT_CACHE}) {
        print STDERR " saving alignment cache $cache_file\n";
        print STDERR " saving alignment cache $cache_file.err\n";
        system "cp -a $alignment_file $cache_file";
        system "cp -a $alignment_file.err $cache_file.err";
        return $cache_file;
    }
}

sub incremental_alignment {
    my ($dclst, @cluster) = @_;

    print STDERR "\n";
    print STDERR "Incremental merge\n";

    if (@cluster != 2) {
        die;
    }
    my $domtbl_1 = extract_domtbl($dclst, $cluster[0]);
    my $domtbl_2 = extract_domtbl($dclst, $cluster[1]);
    my @genes_1 = extract_genes($domtbl_1);
    my @genes_2 = extract_genes($domtbl_2);
    # print STDERR "cluster $cluster[0] (", scalar(@genes_1), " genes): @genes_1\n";
    # print STDERR "cluster $cluster[1] (", scalar(@genes_2), " genes): @genes_2\n";

    my @genes_main = ();
    my @genes_add = ();
    if (@genes_1 > @genes_2) {
        @genes_main = @genes_1;
        @genes_add = subtract_geneset(\@genes_1, \@genes_2);
    } else {
        @genes_main = @genes_2;
        @genes_add = subtract_geneset(\@genes_2, \@genes_1);
    }
    # print STDERR "diff clusters (", scalar(@genes_add), " genes): @genes_add\n";

    my @genes_union = uniq(@genes_1, @genes_2);
    my %seq = ();
    read_seq(\%seq, \@genes_union); # does not care geneset seq (fused seq). c.f. get_seq command

    open(TMP_SEQ_UNION, ">$TMP_SEQ_UNION") || die;
    open(TMP_SEQ_COMMON, ">$TMP_SEQ_COMMON") || die;
    for my $gene (@genes_main) {
        my $seq = get_geneset_seq($gene, \%seq);
        if ($seq) {
            print TMP_SEQ_COMMON $seq;
            print TMP_SEQ_UNION $seq;
        } else {
            print STDERR "WARNING: incremental_alignment: no seq for $gene\n";
        }
    }
    close(TMP_SEQ_COMMON);
    open(TMP_SEQ_ADD, ">$TMP_SEQ_ADD") || die;
    for my $gene (@genes_add) {
        my $seq = get_geneset_seq($gene, \%seq);
        if ($seq) {
            print TMP_SEQ_ADD $seq;
            print TMP_SEQ_UNION $seq;
        } else {
            print STDERR "WARNING: incremental_alignment: no seq for $gene\n";
        }
    }
    close(TMP_SEQ_ADD);
    close(TMP_SEQ_UNION);

    my $aligner = $ENV{DOMREFINE_ALIGNER};
    my $cache_file = cache_file_path($ENV{DOMREFINE_CACHE}, md5_hex(fasta_to_sorted_tsv($TMP_SEQ_UNION)), ".$aligner");
    my $cache_file_stepwise_merge = cache_file_path($ENV{DOMREFINE_CACHE_STEPWISE_MERGE}, md5_hex(fasta_to_sorted_tsv($TMP_SEQ_UNION)), ".$aligner");

    if ($ENV{DOMREFINE_READ_ALIGNMENT_CACHE} and found_alignment_cache($TMP_SEQ_UNION, $cache_file)) {
        print STDERR " found alignment cache $cache_file\n";
        return;
    }
    if (found_alignment_cache($TMP_SEQ_UNION, $cache_file_stepwise_merge)) {
        print STDERR " found alignment cache $cache_file_stepwise_merge\n";
        return;
    }

    print STDERR "Aligning common seq (clusters @cluster)\n";
    aligner($TMP_SEQ_COMMON, $TMP_ALIGNMENT_COMMON, stepwise_merge => 1);
    print STDERR "Aligning additional seq (clusters @cluster)\n";
    my $start_time = time;
    my $aligner_command = "clustalo-1.2.0-Ubuntu-x86_64 --auto --p1 $TMP_ALIGNMENT_COMMON -o $TMP_ALIGNMENT_UNION -l $TMP_ALIGNMENT_UNION.err --force -v -v";
    if ($aligner =~ /kimura$/) {
        $aligner_command .= " --use-kimura";
    }
    if (@genes_add == 0) {
        print STDERR "WARNING: incremental_alignment: no seq diff for cluster pair @cluster\n";
        system "cp -a $TMP_ALIGNMENT_COMMON $TMP_ALIGNMENT_UNION";
    } elsif (@genes_add == 1) {
        print STDERR " add 1 seq: ";
        system "$aligner_command --p2 $TMP_SEQ_ADD";
    } else {
        print STDERR " align ", scalar(@genes_main), " seq and ", scalar(@genes_add), " seq: ";
        system "$aligner_command -i $TMP_SEQ_ADD";
        # using --output-order=tree-order results in segmentation fault
    }
    my $end_time = time;
    printf STDERR " %d sec\n", $end_time - $start_time;

    print STDERR " saving alignment cache $cache_file_stepwise_merge\n";
    print STDERR " saving alignment cache $cache_file_stepwise_merge.err\n";
    system "cp -a $TMP_ALIGNMENT_UNION $cache_file_stepwise_merge";
    system "cp -a $TMP_ALIGNMENT_UNION.err $cache_file_stepwise_merge.err";
}

sub create_alignment {
    my ($dclst_file, $r_gene, $alignment_file, %opt) = @_;

    # prepare seq file
    my %before_seq = ();
    my %after_seq = ();
    my @before_seq = ();
    my @after_seq = ();
    if ($opt{region} or $opt{REGION}) {
        print STDERR "read region seq\n";
        read_dclst_to_seq_file_region($dclst_file, $TMP_SEQ_TO_CREATE_ALIGNMENT, \%before_seq, \%after_seq);
    } elsif ($opt{domain} or $opt{DOMAIN}) {
        print STDERR "read domain seq\n";
        read_dclst_to_seq_file_domain($dclst_file, $TMP_SEQ_TO_CREATE_ALIGNMENT, \@before_seq, \@after_seq);
    } else {
        print STDERR "read gene seq\n";
        read_dclst_to_seq_file($dclst_file, $TMP_SEQ_TO_CREATE_ALIGNMENT);
    }

    # create alignment file
    my @line = `cat $TMP_SEQ_TO_CREATE_ALIGNMENT | grep '^>'`;
    if (@line == 1) {
        system "cp $TMP_SEQ_TO_CREATE_ALIGNMENT $alignment_file";
    } elsif (@line >= 2) {
        aligner($TMP_SEQ_TO_CREATE_ALIGNMENT, $alignment_file, %opt);
    } else {
        die "no seq to align";
    }

    # post process
    if ($opt{region}) {
        append_seq_to_alignments($alignment_file, \%before_seq, \%after_seq);
    } elsif ($opt{domain}) {
        # This option is only used in dom_draw.
        # WARNING: do NOT change the order of sequences in the alignment.
        #          And do not use CACHE (different order of input could be cached as the same calculation).
        append_seq_to_alignments_array($alignment_file, \@before_seq, \@after_seq);
    }
    read_fasta_entries($alignment_file, $r_gene);


    # if ($opt{DOMAIN}) {
    # 	add_domain_number($alignment_file);
    # }
}

sub get_alignment_structure {
    my ($dclst_file, $r_gene, $r_a, $r_b, $r_p, %opt) = @_;

    create_alignment($dclst_file, $r_gene, $TMP_ALIGNMENT_TO_GET_ALIGNMENT_STRUCTURE, %opt);

    my $n_seq = read_alignment($TMP_ALIGNMENT_TO_GET_ALIGNMENT_STRUCTURE, $r_a);
    my $n_pos = check_alignment_length(\@{$r_a});
    if (@{$r_gene} != @{$r_a}) {
        die scalar(@{$r_gene}) . " != " . scalar(@{$r_a});
    }
    print STDERR " n_seq=$n_seq, n_pos=$n_pos\n";

    summarize_amino_acid_frequency($r_a, $r_b, $r_p);

    return (scalar(@{$r_a}), scalar(@{${$r_a}[0]}));
}

sub get_alignment_matrices {
    my ($dclst_file, $tmp_alignment_file, $r_gene, $r_a, $r_b, $r_d, $r_p, $h_cluster, $h_domain, $r_get_j, %opt) = @_;

    my ($n, $m) = get_alignment_structure_from_file($tmp_alignment_file, $r_gene, $r_a, $r_b, $r_p, %opt);
    print STDERR " n_seq=$n n_pos=$m\n";

    get_dclst_structure($dclst_file, $h_cluster, $h_domain);
    if ($opt{domain}) {
        map_domains_one_by_one($r_gene, $h_domain, $r_b, $r_d, $dclst_file);
    } else {
        map_domains($r_gene, $h_domain, $r_b, $r_d);
    }
    create_get_j($r_b, $r_get_j);

    return ($n, $m);
}

sub read_alignment {
    my ($file, $r_a) = @_;

    my $alignment = "";
    open(F, $file) || die;
    while (<F>) {
        chomp;
        if (/^>(.*)/) {
            if ($alignment) {
                $alignment =~ s/\s//g;
                my @a = split //, $alignment;
                push @{$r_a}, \@a;
                $alignment = "";
            }
        } else {
            $alignment .= $_;
        }
    }
    $alignment =~ s/\s//g;
    my @a = split //, $alignment;
    push @{$r_a}, \@a;
    return @{$r_a};
}

sub map_domains {
    my ($r_gene, $r_h_domain, $r_b, $r_d) = @_;

    my $n = @{$r_b};
    my $m = @{${$r_b}[0]};

    initialize_matrix($r_d, $n, $m, 0);
    for (my $i=0; $i<$n; $i++) {
        my $gene = ${$r_gene}[$i];
        my @domains = keys %{${$r_h_domain}{$gene}};
        for my $domain (@domains) {
            my $begin_pos = ${$r_h_domain}{$gene}{$domain}{begin};
            my $end_pos = ${$r_h_domain}{$gene}{$domain}{end};
            my $cluster = ${$r_h_domain}{$gene}{$domain}{cluster};
            my ($begin_j, $end_j) = get_j_of_domain(${$r_b}[$i], $begin_pos, $end_pos);
            if (defined($begin_j) and defined($end_j)) {
                for (my $j=$begin_j; $j<=$end_j; $j++) {
                    ${$r_d}[$i][$j] = $cluster;
                }
            }
        }
    }
}

sub create_get_j {
    my ($r_b, $r_get_j) = @_;

    my $n = @{$r_b};
    my $m = @{${$r_b}[0]};

    for (my $i=0; $i<$n; $i++) {
        my $count = 0;
# 	print STDERR "$i\n";
        for (my $j=0; $j<$m; $j++) {
            if (${$r_b}[$i][$j] == 1) {
                $count ++;
                ${${$r_get_j}[$i]}{$count} = $j;
# 		print STDERR "$count,$j ";
            }
        }
# 	print STDERR "\n";
    }
}

sub create_get_pos_from_a {
    my ($r_a, $r_get_pos, $r_get_j) = @_;

    my $n = @{$r_a};
    my $m = @{${$r_a}[0]};

    for (my $i=0; $i<$n; $i++) {
        my $count = 0;
        for (my $j=0; $j<$m; $j++) {
            if (${$r_a}[$i][$j] ne '-') {
                $count ++;
                ${$r_get_pos}[$i][$j] = $count;
                ${${$r_get_j}[$i]}{$count} = $j;
            }
        }
    }
}

sub fasta_to_sorted_tsv {
    my ($fasta) = @_;

    my @tsv = ();
    my $tmp = "";
    open(FASTA_TO_TSV, $fasta) || die;
    while (my $line = <FASTA_TO_TSV>) {
        chomp($line);
        if ($line =~ /^>/) {
            if ($tmp ne "") {
                push @tsv, $tmp;
            }
            $tmp = "$line\t";
        } else {
            $line =~ s/[-\s]//g;
            $tmp .= $line;
        }
    }
    push @tsv, $tmp;
    close(FASTA_TO_TSV);

    @tsv = sort { $a cmp $b } @tsv;
    return join("\n", @tsv);
}

################################################################################
### Private functions ##########################################################
################################################################################

sub create_get_pos {
    my ($r_b, $r_get_pos, $r_get_len) = @_;

    my $n = @{$r_b};
    my $m = @{${$r_b}[0]};

    initialize_matrix($r_get_pos, $n, $m, ""); # Is this necessary?
    for (my $i=0; $i<$n; $i++) {
        my $count = 0;
        for (my $j=0; $j<$m; $j++) {
            if (${$r_b}[$i][$j] == 1) {
                $count ++;
                ${$r_get_pos}[$i][$j] = $count;
            } else {
                ${$r_get_pos}[$i][$j] = 0; # if it is gap
            }
        }
        ${$r_get_len}[$i] = $count;
    }
}

sub found_alignment_cache {
    my ($seq_file, $alignment_cache) = @_;

    if (! -f $alignment_cache) {
        return 0;
    }

    my $seq = fasta_to_sorted_tsv($seq_file);
    my $alignment = fasta_to_sorted_tsv($alignment_cache);
    if ($seq ne $alignment) {
        print STDERR "WARNING: inconsistent alignment_cache $alignment_cache\n";
        return 0;
    }
    
    return 1;
}

sub convert_seq_file {
    my ($seq_file1, $seq_file2) = @_;

    my $seq = `cat $seq_file1`;
    $seq =~ s/^>.*\s(\S+:\S+)\s.*/>$1/mg;

    open(SEQ_FILE2, "> $seq_file2") || die;
    print SEQ_FILE2 $seq;
    close(SEQ_FILE2);
}

sub convert_alignment_file {
    my ($file1, $file2) = @_;

    my $alignment = `cat $file1`;
    $alignment =~ s/^>(\w+?)_(\S+)/>$1:$2/mg;
    $alignment =~ s/^%(\w+?)_(\S+)/>$1:$2/mg;

    open(ALIGNMENT_FILE2, "> $file2") || die;
    print ALIGNMENT_FILE2 $alignment;
    close(ALIGNMENT_FILE2);
}

sub append_seq_to_alignments {
    my ($tmp_alignment, $r_before_seq, $r_after_seq) = @_;

    my @alignment = ();
    read_alignment_list($tmp_alignment, \@alignment);

    my @gene = ();
    read_fasta_entries($tmp_alignment, \@gene);

    # additional lengths needed (before and after)
    my $max_before_len = 0;
    for (my $i=0; $i<@alignment; $i++) {
        if ($alignment[$i] =~ /^(-*)/) {
            my $start_gap_len = length($1);
            my $before_seq_len = length(${$r_before_seq}{$gene[$i]});
            if ($before_seq_len > $start_gap_len) {
                my $required_len = $before_seq_len - $start_gap_len;
                if ($required_len > $max_before_len) {
                    $max_before_len = $required_len;
                }
            }
        }
    }
    my $max_after_len = 0;
    for (my $i=0; $i<@alignment; $i++) {
        if ($alignment[$i] =~ /[^-](-*)$/) {
            my $end_gap_len = length($1);
            my $after_seq_len = length(${$r_after_seq}{$gene[$i]});
            if ($after_seq_len > $end_gap_len) {
                my $required_len = $after_seq_len - $end_gap_len;
                if ($required_len > $max_after_len) {
                    $max_after_len = $required_len;
                }
            }
        }
    }

    # append before_seq and after_seq
    open(TMP_ALIGNMENT, ">$tmp_alignment") || die;
    for (my $i=0; $i<@alignment; $i++) {
        my $alignment = "-" x $max_before_len . $alignment[$i] . "-" x $max_after_len;
        if ($alignment =~ /^(-*)/) {
            my $start_gap_len = length($1);
            $alignment =~ s/^-*//;
            my $start_gap_len_required = $start_gap_len - length(${$r_before_seq}{$gene[$i]});
            $alignment = "-" x $start_gap_len_required . ${$r_before_seq}{$gene[$i]} . $alignment;
        }
        if ($alignment =~ /[^-](-*)$/) {
            my $end_gap_len = length($1);
            $alignment =~ s/([^-])-*$/$1/;
            my $end_gap_len_required = $end_gap_len - length(${$r_after_seq}{$gene[$i]});
            $alignment = $alignment . ${$r_after_seq}{$gene[$i]} . "-" x $end_gap_len_required;
        }
        print TMP_ALIGNMENT ">", $gene[$i], "\n";
        print TMP_ALIGNMENT $alignment, "\n";
    }
    close(TMP_ALIGNMENT);
}

sub append_seq_to_alignments_array {
    my ($tmp_alignment, $r_before_seq, $r_after_seq) = @_;

    my @alignment = ();
    read_alignment_list($tmp_alignment, \@alignment);

    my @gene = ();
    read_fasta_entries($tmp_alignment, \@gene);

    # additional lengths needed (before and after)
    my $max_before_len = 0;
    for (my $i=0; $i<@alignment; $i++) {
        if ($alignment[$i] =~ /^(-*)/) {
            my $start_gap_len = length($1);
            my $before_seq_len = length(${$r_before_seq}[$i]);
            if ($before_seq_len > $start_gap_len) {
                my $required_len = $before_seq_len - $start_gap_len;
                if ($required_len > $max_before_len) {
                    $max_before_len = $required_len;
                }
            }
        }
    }
    my $max_after_len = 0;
    for (my $i=0; $i<@alignment; $i++) {
        if ($alignment[$i] =~ /[^-](-*)$/) {
            my $end_gap_len = length($1);
            my $after_seq_len = length(${$r_after_seq}[$i]);
            if ($after_seq_len > $end_gap_len) {
                my $required_len = $after_seq_len - $end_gap_len;
                if ($required_len > $max_after_len) {
                    $max_after_len = $required_len;
                }
            }
        }
    }

    # append before_seq and after_seq
    open(TMP_ALIGNMENT, ">$tmp_alignment") || die;
    for (my $i=0; $i<@alignment; $i++) {
        my $alignment = "-" x $max_before_len . $alignment[$i] . "-" x $max_after_len;
        if ($alignment =~ /^(-*)/) {
            my $start_gap_len = length($1);
            $alignment =~ s/^-*//;
            my $start_gap_len_required = $start_gap_len - length(${$r_before_seq}[$i]);
            $alignment = "-" x $start_gap_len_required . ${$r_before_seq}[$i] . $alignment;
        }
        if ($alignment =~ /[^-](-*)$/) {
            my $end_gap_len = length($1);
            $alignment =~ s/([^-])-*$/$1/;
            my $end_gap_len_required = $end_gap_len - length(${$r_after_seq}[$i]);
            $alignment = $alignment . ${$r_after_seq}[$i] . "-" x $end_gap_len_required;
        }
        print TMP_ALIGNMENT ">", $gene[$i], "\n";
        print TMP_ALIGNMENT $alignment, "\n";
    }
    close(TMP_ALIGNMENT);
}

sub read_alignment_list {
    my ($file, $r_a) = @_;

    my $alignment = "";
    open(TMP_ALIGNMENT, $file) || die;
    while (<TMP_ALIGNMENT>) {
        chomp;
        if (/^>(.*)/) {
            if ($alignment) {
                $alignment =~ s/\s//g;
                push @{$r_a}, $alignment;
                $alignment = "";
            }
        } else {
            $alignment .= $_;
        }
    }
    $alignment =~ s/\s//g;
    push @{$r_a}, $alignment;
    close(TMP_ALIGNMENT);
    return @{$r_a};
}

sub get_alignment_structure_from_file {
    my ($tmp_alignment_file, $r_gene, $r_a, $r_b, $r_p, %opt) = @_;

    my $r_gene_list = $opt{gene_list};
    my $n_seq;
    if ($r_gene_list and %{$r_gene_list}) {
        $n_seq = read_alignment_for_gene($tmp_alignment_file, $r_a, $r_gene, $r_gene_list);
    } else {
        $n_seq = read_alignment($tmp_alignment_file, $r_a);
    }
    my $n_pos = check_alignment_length(\@{$r_a});
    if (@{$r_gene} != @{$r_a}) {
        die "n_gene=" . scalar(@{$r_gene}) . ", n_row=" . scalar(@{$r_a});
    }
    print STDERR " n_seq=$n_seq, n_pos=$n_pos\n";

    summarize_amino_acid_frequency($r_a, $r_b, $r_p);

    return (scalar(@{$r_a}), scalar(@{${$r_a}[0]}));
}

sub map_domains_one_by_one {
    my ($r_gene, $r_h_domain, $r_b, $r_d, $dclst_file) = @_;

    my $n = @{$r_b};
    my $m = @{${$r_b}[0]};

    my @domain = `cat $dclst_file | cut.sh 3`;
    chomp(@domain);

    initialize_matrix($r_d, $n, $m, 0);
    for (my $i=0; $i<$n; $i++) {
        my $gene = ${$r_gene}[$i];
        my $begin_pos = ${$r_h_domain}{$gene}{$domain[$i]}{begin};
        my $end_pos = ${$r_h_domain}{$gene}{$domain[$i]}{end};
        my $cluster = ${$r_h_domain}{$gene}{$domain[$i]}{cluster};
        my ($begin_j, $end_j) = get_j_of_domain(${$r_b}[$i], $begin_pos, $end_pos);
        if (defined($begin_j) and defined($end_j)) {
            for (my $j=$begin_j; $j<=$end_j; $j++) {
                ${$r_d}[$i][$j] = $cluster;
            }
        }
    }
}

sub get_j_of_domain {
    my ($r_b, $begin_pos, $end_pos) = @_;
    
    my $m = @{$r_b};

    my ($begin_j, $end_j);
    my ($seq_begin_j, $seq_end_j);
    
    my $count = 0;
    for (my $j=0; $j<$m; $j++) {
        if (${$r_b}[$j] == 1) {
            $count ++;

            if ($count == 1) {
                $seq_begin_j = $j;
            }
            $seq_end_j = $j;
            
            if ($begin_pos and $begin_pos == $count) {
                $begin_j = $j;
            }
            if ($end_pos and $end_pos == $count) {
                $end_j = $j;
            }
        }
    }

    if (defined($begin_j) and defined($end_j)) {
        return ($begin_j, $end_j);
    } elsif (defined($begin_j) and ! defined($end_j)) {
        return ($begin_j, $seq_end_j);
    } elsif (! defined($begin_j) and defined($end_j)) {
        return ($seq_begin_j, $end_j);
    } elsif (! defined($begin_j) and ! defined($end_j)) {
        return;
    } else {
        die;
    }
}

sub read_alignment_for_gene {
    my ($file, $r_a, $r_gene, $r_gene_list) = @_;

    @{$r_gene} = ();
    my $alignment = "";
    my $entry_id = "";
    open(F, $file) || die;
    while (<F>) {
        chomp;
        if (/^>(\S+)/) {
            $entry_id = $1;
            if ($alignment) {
                $alignment =~ s/\s//g;
                my @a = split //, $alignment;
                push @{$r_a}, \@a;
                $alignment = "";
            }
            if (${$r_gene_list}{$entry_id}) {
                push @{$r_gene}, $entry_id;
            }
        } elsif (/^>/) {
            die;
        } else {
            if ($entry_id ne "") {
                if (${$r_gene_list}{$entry_id}) {
                    $alignment .= $_;
                }
            } else {
                die;
            }
        }
    }
    if ($alignment) {
        $alignment =~ s/\s//g;
        my @a = split //, $alignment;
        push @{$r_a}, \@a;
    }
    return @{$r_a};
}

sub check_alignment_length {
    my ($r_a) = @_;

    my $length = @{${$r_a}[0]};
    for (my $i=0; $i<@{$r_a}; $i++) {
        my $length2 = @{${$r_a}[$i]};
        if ($length2 != $length) {
            die "$length != $length2";
        }
    }

    return $length;
}

sub summarize_amino_acid_frequency {
    my ($r_a, $r_b, $r_p) = @_;

    my $n = @{$r_a};
    my $m = @{${$r_a}[0]};

    my %f = ();
    for (my $j=0; $j<$m; $j++) {
# 	my $count_b = 0;
        for (my $i=0; $i<$n; $i++) {
            if (${$r_a}[$i][$j] ne "-") {
                ${$r_b}[$i][$j] = 1;
# 		$count_b ++;
            } else {
                ${$r_b}[$i][$j] = 0;
            }
            $f{$j}{${$r_a}[$i][$j]}++;
        }
        if (! $f{$j}{"-"}) {
            $f{$j}{"-"} = 0;
        }
        for my $a (keys %{$f{$j}}) {
            ${$r_p}{$j}{$a} = $f{$j}{$a}/$n;
# 	    ${$r_p}{$j}{$a} = $f{$j}{$a}/$count_b;
        }
    }

    # check terminal gaps
    for (my $i=0; $i<$n; $i++) {
        my $j = 0;
        while (${$r_a}[$i][$j] eq "-" and $j<$m) {
            ${$r_b}[$i][$j] = -1;
            $j ++;
        }
        $j = $m - 1;
        while (${$r_a}[$i][$j] eq "-" and $j>=0) {
            ${$r_b}[$i][$j] = -1;
            $j --;
        }
    }

}

sub add_domain_number {
    my ($tmp_alignment) = @_;

    my @alignment = `cat $tmp_alignment`;
    chomp(@alignment);

    for (my $i=0; $i<@alignment; $i++) {
        if ($alignment[$i] =~ /^>\S+/) {
            if ($alignment[$i] =~ /^>(\S+) (\S+)/) {
                my ($gene, $domain) = ($1, $2);
                if ($domain eq "0") {
                    $alignment[$i] = ">$gene";
                } else {
                    $alignment[$i] = ">$gene($domain)";
                }
            } else {
                die $alignment[$i];
            }
        }
    }
    
    # overwite tmp_alignment
    open(TMP_ALIGNMENT, ">$tmp_alignment") || die;
    print TMP_ALIGNMENT join("\n", @alignment), "\n";
    close(TMP_ALIGNMENT);
}

sub subtract_geneset {
    my ($r_genes_1, $r_genes_2) = @_;

    my %checked = ();
    for my $gene (@{$r_genes_1}) {
        $checked{$gene} = 1;
    }

    my @diff = ();
    for my $gene (@{$r_genes_2}) {
        if ($checked{$gene}) {
        } else {
            push @diff, $gene;
        }
    }

    return @diff;
}

1;
