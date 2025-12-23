package DomRefine::Read;
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(parse_cluster_members get_dclst_structure decompose_dclst_line parse_dclst_to_structure
             output_clusters output_domains get_dclst_of_domain
             extract_domtbl extract_dclst extract_dclst_compl extract_dclst_by_genes
             read_dclst_geneset extract_genes read_dclst_seq
             read_dclst_to_seq_file read_dclst_to_seq_file_domain read_dclst_to_seq_file_region
             decompose_gene_id get_rep_cluster choose_rep_cluster
             get_annotation get_clusters
             cache_file_path define_tmp_file remove_tmp_file
             save_stdin save_contents
             get_gene_idx get_seq_fast read_seq get_dbh
             get_geneset_descr get_genes_fused get_geneset_seq
             read_cluster_members read_cluster_members_from_file
             output_cluster_members get_cluster_count read_fasta_entries
             get_gene_length calc_total_aa
             read_homology_data put_homology_in_hash
             );

use strict;
use Digest::MD5 qw(md5_hex);
use DomRefine::General;

$ENV{DOMREFINE_CACHE} || die "set environmental DOMREFINE_CACHE\n";
$ENV{DOMREFINE_TMP} || die "set environmental DOMREFINE_TMP\n";

my $TMP_GET_A_CLUSTER = define_tmp_file("get_a_cluster");
my $TMP_DCLST_TO_GET_STRUCTURE = define_tmp_file("dclst_to_get_structure");
END {
    remove_tmp_file($TMP_GET_A_CLUSTER);
}

################################################################################
### Public functions ###########################################################
################################################################################

sub parse_cluster_members {
    my ($fh, $r_domain, $r_cluster_member) = @_;

    while (my $line = <$fh>) {
        chomp($line);
        if ($line eq "") {
            print STDERR "read blank line $.\n";
            next;
        }
        my ($cluster, $gene, @domain_info) = decompose_dclst_line($line);
        unless (@domain_info and @domain_info % 3 == 0) {
            die $line;
        }
        for (my $i=0; $i<@domain_info; $i+=3) {
            my ($domain, $begin_pos, $end_pos) = ($domain_info[$i], $domain_info[$i+1], $domain_info[$i+2]);
            ${$r_domain}{$gene}{$domain}{cluster} = $cluster;
            ${$r_domain}{$gene}{$domain}{begin} = $begin_pos;
            ${$r_domain}{$gene}{$domain}{end} = $end_pos;
            ${$r_cluster_member}{$cluster}{"$gene($domain)"} = 1; ### this is different from get_dclst_structure
        }
    }
}

sub get_dclst_structure {
    my ($dclst_table_file, $h_cluster, $h_domain) = @_;

    open(DCLST, $dclst_table_file) || die "$dclst_table_file: $!";
    while (my $line = <DCLST>) {
        chomp($line);
        if ($line eq "") {
            print STDERR "read blank line $.\n";
            next;
        }
        my ($cluster, $gene, @domain_info) = decompose_dclst_line($line);
        unless (@domain_info and @domain_info % 3 == 0) {
            die $line;
        }
        for (my $i=0; $i<@domain_info; $i+=3) {
            my ($domain, $begin_pos, $end_pos) = ($domain_info[$i], $domain_info[$i+1], $domain_info[$i+2]);
            ${$h_domain}{$gene}{$domain}{cluster} = $cluster;
            ${$h_domain}{$gene}{$domain}{begin} = $begin_pos;
            ${$h_domain}{$gene}{$domain}{end} = $end_pos;
            ${$h_cluster}{$cluster}{$gene} ++;
        }
    }
    close(DCLST);
}

sub decompose_dclst_line {
    my ($dclst_line) = @_;

    chomp($dclst_line);
    if ($dclst_line =~ /^(\S+)\s(\S+)\s(\S+)\s(\d+)\s(\d+)/) {
        my @f = split(/\s+/, $dclst_line);
        if ($ENV{DOMREFINE_PRECLUST_INFO}) {
            my ($cluster, $gene, $domain, $start, $end) = @f;
            return ($cluster, $gene, $domain, $start, $end);
        } else {
            return @f;
        }
    } elsif ($dclst_line =~ /^(\S+)\s(\S+)\s(\d+)\s(\d+)$/) { # eggNOG member format
        my ($cluster, $gene, $start, $end) = ($1, $2, $3, $4);
        my $domain = 0;
        return ($cluster, $gene, $domain, $start, $end);
    }
}

sub parse_dclst_to_structure { # obsolete
    my ($dclst_table, $h_cluster, $h_domain) = @_;

    save_contents($dclst_table, $TMP_DCLST_TO_GET_STRUCTURE);
    get_dclst_structure($TMP_DCLST_TO_GET_STRUCTURE, $h_cluster, $h_domain);
}

sub output_clusters {
    my ($r_cluster, $r_domain) = @_;

    for my $cluster (sort {$a cmp $b} keys %{$r_cluster}) {
        for my $gene (sort {$a cmp $b} keys %{${$r_cluster}{$cluster}}) {
            for my $domain (sort {$a <=> $b} keys %{${$r_domain}{$gene}}) {
                my $cluster_for_this_domain = ${$r_domain}{$gene}{$domain}{cluster};
                my $begin = ${$r_domain}{$gene}{$domain}{begin};
                my $end = ${$r_domain}{$gene}{$domain}{end};
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

sub output_domains {
    my ($r_domain) = @_;

    for my $gene (sort {$a cmp $b} keys %{$r_domain}) {
        for my $domain (sort {$a <=> $b} keys %{${$r_domain}{$gene}}) {
            my $cluster = ${$r_domain}{$gene}{$domain}{cluster};
            my $begin = ${$r_domain}{$gene}{$domain}{begin};
            my $end = ${$r_domain}{$gene}{$domain}{end};
            print "$cluster $gene $domain $begin $end\n";
        }
    }
}

sub get_dclst_of_domain {
    my ($r_domain) = @_;

    my $dclst = "";
    for my $gene (sort {$a cmp $b} keys %{$r_domain}) {
        for my $domain (sort {$a <=> $b} keys %{${$r_domain}{$gene}}) {
            my $cluster = ${$r_domain}{$gene}{$domain}{cluster};
            my $begin = ${$r_domain}{$gene}{$domain}{begin};
            my $end = ${$r_domain}{$gene}{$domain}{end};
            $dclst .= "$cluster $gene $domain $begin $end\n";
        }
    }
    return $dclst;
}

sub extract_domtbl {
    my ($domtbl, @cluster) = @_;

    my %target = ();
    for my $cluster (@cluster) {
        $target{$cluster} = 1;
    }
    my %extracted = ();
    
    chomp($domtbl);
    my @line = split("\n", $domtbl);
    my @extracted = ();
    for my $line (@line) {
        if ($line =~ /./) {
            my ($cluster_id) = split(/\s+/, $line);
            if (! defined $cluster_id) {
                die;
            }
            if ($target{$cluster_id}) {
                push @extracted, $line;
            }
            $extracted{$cluster_id} = 1;
        }
    }

    for my $cluster (@cluster) {
        if (! $extracted{$cluster}) {
            print STDERR "[$cluster] is not containded.\n";
        }
    }

    return join("\n", @extracted) . "\n";
}

sub extract_dclst {
    my ($dclst_table_file, @cluster) = @_;
    
    my %hash = ();
    for my $cluster (@cluster) {
        $hash{$cluster} = 1;
    }

    my $extracted_dclst = "";
    open(DCLST_TABLE_FILE, $dclst_table_file) || die;
    while (<DCLST_TABLE_FILE>) {
        unless (/./) {
            next;
        }
        my ($cluster_id) = split;
        if (! defined $cluster_id) {
            die;
        }
        if ($hash{$cluster_id}) {
            $extracted_dclst .= $_;
        }
    }
    close(DCLST_TABLE_FILE);

    return $extracted_dclst;
}

sub extract_dclst_compl {
    my ($dclst_file, @cluster) = @_;
    
    my %hash = ();
    for my $cluster (@cluster) {
        $hash{$cluster} = 1;
    }

    my $dclst = "";
    open(EXTRACT_DCLST_COMPL, $dclst_file) || die;
    while (<EXTRACT_DCLST_COMPL>) {
        my ($cluster_id) = split;
        if ($hash{$cluster_id}) {
        } else {
            $dclst .= $_;
        }
    }
    close(EXTRACT_DCLST_COMPL);

    # my $pattern =  join("|", @cluster);
    # my $dclst = `cat $dclst_file | grep -P -v '^($pattern) '`;

    return $dclst;
}

sub extract_dclst_by_genes {
    my ($dclst_table_file, $r_select) = @_;

    my $extracted_dclst = "";

    my @line = `cat $dclst_table_file`;
    for my $line (@line) {
        my ($cluster, $gene) = split(/\s+/, $line);
        if (${$r_select}{$gene}) {
            $extracted_dclst .= $line;
        }
    }

    return $extracted_dclst;
}

sub read_dclst_geneset {
    my ($dclst_file) = @_;

    my @geneset = `cat $dclst_file | cut.sh 2`;
    chomp(@geneset);

    return uniq(@geneset);
}

sub extract_genes { # cannot tread fused genes
    my ($dclst) = @_;

    my @line = split("\n", $dclst);
    chomp(@line);

    my %gene = ();
    for my $line (@line) {
        if ($line =~ /^\S+\s([A-Za-z0-9]+[:\.]\S+)\s/) {
            $gene{$1} = 1;
        } else {
            die;
        }
    }

    return keys %gene;
}

sub read_dclst_seq {
    my ($dclst_file, $r_seq) = @_;

    my @geneset = read_dclst_geneset($dclst_file);
    my @gene = get_genes_fused(\@geneset);
    read_seq($r_seq, \@gene);
}

sub read_dclst_to_seq_file {
    my ($dclst_file, $tmp_seq) = @_;

    my %seq = ();
    my @geneset = read_dclst_geneset($dclst_file);
    my @gene = get_genes_fused(\@geneset);
    read_seq(\%seq, \@gene);

    open(TMP_SEQ, ">$tmp_seq") || die;
    for my $geneset (@geneset) {
        my $seq = get_geneset_seq($geneset, \%seq);
        if ($seq) {
            print TMP_SEQ $seq;
        }
    }
    close(TMP_SEQ);
}

sub read_dclst_to_seq_file_region {
    my ($dclst_file, $tmp_seq, $r_before_seq, $r_after_seq) = @_;

    my %seq = ();
    my @geneset = read_dclst_geneset($dclst_file);
    my @gene = get_genes_fused(\@geneset);
    read_seq(\%seq, \@gene);

    my @dclst_line = `cat $dclst_file`;
    open(TMP_SEQ, ">$tmp_seq") || die "cannot create $tmp_seq";
    for (my $i=0; $i<@geneset; $i++) {
        my $seq = get_geneset_seq($geneset[$i], \%seq);
        if ($seq) {
            my ($start, $end) = extract_region($dclst_file, $geneset[$i]);
            my ($header, $region_seq, $before_seq, $after_seq) = trim_to_region($seq, $start, $end);
            ${$r_before_seq}{$geneset[$i]} = $before_seq;
            ${$r_after_seq}{$geneset[$i]} = $after_seq;
            if ($region_seq) {
                print TMP_SEQ "$header\n$region_seq\n";
            }
        } else {
            print STDERR "WARNING: cant get seq for $geneset[$i]\n";
        }
    }
    close(TMP_SEQ);
}

sub read_dclst_to_seq_file_domain {
    my ($dclst_file, $tmp_seq, $r_before_seq, $r_after_seq) = @_;
    
    my %seq = ();
    my @geneset = `cat $dclst_file | cut.sh 2`;
    chomp(@geneset);
    my @gene = get_genes_fused(\@geneset);
    read_seq(\%seq, \@gene);

    my @dclst_line = `cat $dclst_file`;
    open(TMP_SEQ, ">$tmp_seq") || die;
    for (my $i=0; $i<@geneset; $i++) {
        my $seq = get_geneset_seq($geneset[$i], \%seq);
        if ($seq) {
            my ($cluster, $gene, $domain, $start, $end) = decompose_dclst_line($dclst_line[$i]);
            my ($header, $region_seq, $before_seq, $after_seq) = trim_to_region($seq, $start, $end);
            ${$r_before_seq}[$i] = $before_seq;
            ${$r_after_seq}[$i] = $after_seq;
            if ($region_seq) {
                # print TMP_SEQ "$header $domain $start $end\n";
                if ($domain eq "0") {
                    # print TMP_SEQ "$header $start $end\n";
                    print TMP_SEQ "$header\n";
                } else {
                    # print TMP_SEQ "$header($domain) $start $end\n";
                    print TMP_SEQ "$header($domain)\n";
                }
                print TMP_SEQ "$region_seq\n";
            } else {
                print STDERR "ERROR: could not get sequence of $gene($domain), thus it is omitted.\n";
            }
        } else {
            print STDERR "WARNING: cant get seq for $geneset[$i]\n";
        }
    }
    close(TMP_SEQ);
}

sub get_genes_fused {
    my ($r_geneset) = @_;

    my @gene = ();
    for my $geneset (@{$r_geneset}) {
        for my $gene (split(/\|/, $geneset)) {
            push @gene, $gene;
        }
    }

    return @gene;
}

sub get_clusters { ### Maybe slow; should be replaced with hash
    my ($dclst_file, $sort_option) = @_;
    $sort_option ||= "";
    
    my @cluster = `cat $dclst_file | cut.sh 1 | sort -n $sort_option | uniq`; # inappropreate for cluster_id of string
    chomp(@cluster);

    return @cluster;
}

sub get_a_cluster_with_check {
    my ($dclst) = @_;

    save_contents($dclst, $TMP_GET_A_CLUSTER);

    my @cluster = get_clusters($TMP_GET_A_CLUSTER);
    if (@cluster != 1) {
        die;
    }

    return $cluster[0];
}

sub wrap_lines_for_same_cluster_and_gene {
    my ($input, $output) = @_;
    
    open(INPUT, $input) || die;
    my %hash = ();
    while (<INPUT>) {
        chomp;
        my ($cluster_id, $gene, @domain_info) = split;
        unless (@domain_info and @domain_info % 3 == 0) {
            die;
        }
        for (my $i=0; $i<@domain_info; $i+=3) {
            my ($domain, $begin_pos, $end_pos) = ($domain_info[$i], $domain_info[$i+1], $domain_info[$i+2]);
            $hash{$cluster_id}{$gene}{$domain}{begin} = $begin_pos;
            $hash{$cluster_id}{$gene}{$domain}{end} = $end_pos;
        }
    }
    close(INPUT);

    open(OUTPUT, ">$output") || die;;
    for my $cluster_id (sort {$a cmp $b} keys %hash) {
        for my $gene (sort {$a cmp $b} keys %{$hash{$cluster_id}}) {
            my @domain_info = ();
            for my $domain (sort {$a<=>$b} keys %{$hash{$cluster_id}{$gene}}) {
                my $begin = $hash{$cluster_id}{$gene}{$domain}{begin};
                my $end = $hash{$cluster_id}{$gene}{$domain}{end};
                push @domain_info, $domain, $begin, $end;
            }
            print OUTPUT "$cluster_id $gene @domain_info\n";
        }
    }
    close(OUTPUT);
}

sub get_rep_cluster {
    my ($dclst) = @_;

    my @cluster = ();
    my $flg_string = 0;
    while ($dclst =~ /^(\S+) /gm) {
        my $cluster = $1;
        push @cluster, $cluster;
        if ($cluster !~ /^\d+$/) {
            $flg_string = 1;
        }
    }

    my @cluster_sorted = ();
    if ($flg_string) {
        @cluster_sorted = sort {$a cmp $b} @cluster;
    } else {
        @cluster_sorted = sort {$a <=> $b} @cluster;
    }

    return $cluster_sorted[0];
}

sub choose_rep_cluster {
    my @cluster = @_;

    my $flg_string = 0;
    for my $cluster (@cluster) {
        if ($cluster !~ /^\d+$/) {
            $flg_string = 1;
        }
    }

    my @cluster_sorted = ();
    if ($flg_string) {
        @cluster_sorted = sort {$a cmp $b} @cluster;
    } else {
        @cluster_sorted = sort {$a <=> $b} @cluster;
    }

    return $cluster_sorted[0];
}

sub read_cluster_members {
    my ($dclst, $r_member) = @_;

    for my $line (split("\n", $dclst)) {
        my ($cluster, $gene, $domain, $start, $end) = split(/\s+/, $line);
        if (defined $end) {
            ${$r_member}{$cluster}{$gene}{$domain}{start} = $start;
            ${$r_member}{$cluster}{$gene}{$domain}{end} = $end;
        } else {
            print STDERR "invalid line: $line\n";
        }
    }
}

sub read_cluster_members_from_file {
    my ($dclst_file, $r_member_txt, $r_member) = @_;
    
    open(DCLST_TO_READ_MEMBERS, $dclst_file) || die;
    while (my $line = <DCLST_TO_READ_MEMBERS>) {
        chomp($line);
        if ($line =~ /^(\S+)\s(\S.*)/) {
            my ($cluster, $member) = ($1, $2);
            if (! defined ${$r_member_txt}{$cluster}) {
                ${$r_member_txt}{$cluster} = "$member\n";
            } else {
                ${$r_member_txt}{$cluster} .= "$member\n";
            }
            if (defined $r_member) {
                my ($cluster, $gene, $domain, $start, $end) = split(/\s+/, $line);
                if (defined $end) {
                    ${$r_member}{$cluster}{$gene}{$domain}{start} = $start;
                    ${$r_member}{$cluster}{$gene}{$domain}{end} = $end;
                } else {
                    print STDERR "invalid line: $line\n";
                }
            }
        } else {
            print STDERR "invalid line: $line\n";
        }
    }
    close(DCLST_TO_READ_MEMBERS);
}

sub output_cluster_members {
    my ($r_member) = @_;

    for my $cluster (sort { $a cmp $b } keys(%{$r_member})) {
        my $members = ${$r_member}{$cluster};
        chomp($members);
        my @member = split("\n", $members);
        for my $member (@member) {
            print "$cluster $member\n";
        }
    }
}

sub get_cluster_count {
    my ($r_cluster_count, $r_cluster) = @_;

    for my $cluster (keys %{$r_cluster}) {
        for my $gene (keys %{${$r_cluster}{$cluster}}) {
            ${$r_cluster_count}{$cluster} += ${$r_cluster}{$cluster}{$gene};
        }
    }
}

sub calc_total_aa {
    my ($r_domain) = @_;

    my $total_aa = 0;
    for my $gene (%{$r_domain}) {
        for my $domain (keys %{${$r_domain}{$gene}}) {
            my $begin = ${$r_domain}{$gene}{$domain}{begin};
            my $end = ${$r_domain}{$gene}{$domain}{end};
            my $len = ($end - $begin + 1);
            $total_aa += $len;
        }
    }

    return $total_aa;
}

sub read_seq {
    my ($r_seq, $r_gene) = @_;

    my %seq_to_read = ();
    my $gene_cnt = 0 ;
    for my $gene (@{$r_gene}) {
        my ($sp, $name) = decompose_gene_id($gene);
        $seq_to_read{$sp}{$gene} = 1;
        $gene_cnt++;
    }

    unless ($ENV{DOMREFINE_SEQ_DB} and $ENV{DOMREFINE_SEQ_DB} ne "") {
        print STDERR "ERROR: set DOMREFINE_SEQ_DB\n";
        exit;
    }

    my $start_time = time;
    if (-f $ENV{DOMREFINE_SEQ_DB}) {
        print STDERR " read $ENV{DOMREFINE_SEQ_DB} ..\n";
        if ($ENV{DOMREFINE_BLASTDBCMD} && $gene_cnt < 5000){
            for my $gene (@{$r_gene}) {
                my ($sp, $name) = decompose_gene_id($gene);
                #my $cmd="/bio/bin/blastdbcmd -outfmt '%s' -entry $gene -db $ENV{DOMREFINE_SEQ_DB}";
                #2.11.0では処理速度が遅くなるので、versionを指定
                my $cmd="/bio/package/blast+/2.10.1/bin/blastdbcmd -outfmt '%s' -entry $gene -db $ENV{DOMREFINE_SEQ_DB}";
                my $tmpout=`$cmd`;
                $tmpout=~s/\r\n|\r|\n//g;
                ${$r_seq}{$sp}{$gene} = $tmpout;
            }
        } else {
            open(BLDB, "$ENV{DOMREFINE_SEQ_DB}") || die;
            my $read_flg = 0;
            my $gene;
            my $organism_to_read;
            while (<BLDB>) {
                # get gene and set read_flg
                if (/^>(\S+)/) {
                    $gene = $1;
                    $read_flg = 0;
                    for my $organism (keys %seq_to_read) {
                        if ($seq_to_read{$organism}{$gene}) {
                            # duplicated sequence will be error
                            if (${$r_seq}{$organism}{$gene}) {
                                die;
                            }
                            $read_flg ++;
                            $organism_to_read = $organism;
                            ${$r_seq}{$organism_to_read}{$gene} = "";
                        }
                    }
                    # error handling
                    if ($read_flg >= 2) {
                        die;
                    }
                }
                if ($read_flg) {
                    ${$r_seq}{$organism_to_read}{$gene} .= $_;
                }
            }
            close(BLDB);
        }
    } elsif (-d $ENV{DOMREFINE_SEQ_DB}) {
        print STDERR " read directory $ENV{DOMREFINE_SEQ_DB} ..\n";
        for my $organism (keys %seq_to_read) {
            my $seq_file = "";
            if (-f "$ENV{DOMREFINE_SEQ_DB}/$organism") {
                $seq_file = "$ENV{DOMREFINE_SEQ_DB}/$organism";
            } elsif (-f "$ENV{DOMREFINE_SEQ_DB}/$organism.fa") {
                $seq_file = "$ENV{DOMREFINE_SEQ_DB}/$organism.fa";
            } else {
                die "can't open $organism";
            }
            open(BLDB, "$seq_file") || die "can't open $seq_file";
            my $read_flg = 0;
            my $gene_to_read;
            while (<BLDB>) {
                if (/^>(\S+)/) {
                    my $gene = $1;
                    if ($gene =~ /^(\S+):\S+:(\S+)$/) {
                        $gene = "$1:$2";
                    }
                    $read_flg = 0;
                    if ($seq_to_read{$organism}{$gene}) {
                        # duplicated sequence will be error
                        if (${$r_seq}{$organism}{$gene}) {
                            die;
                        }
                        $read_flg = 1;
                        $gene_to_read = $gene;
                        ${$r_seq}{$organism}{$gene_to_read} = "";
                    }
                }
                if ($read_flg) {
                    ${$r_seq}{$organism}{$gene_to_read} .= $_;
                }
            }
            close(BLDB);
        }
    } else {
        die "\nInvalid value for DOMREFINE_SEQ_DB: $ENV{DOMREFINE_SEQ_DB}";
    }
    my $end_time = time;
    printf STDERR " %d sec\n", $end_time - $start_time;
}

sub get_gene_length {
    my ($h_seq, $gene) = @_;

    my ($sp, $name) = decompose_gene_id($gene);
    my $seq = ${$h_seq}{$sp}{$gene};
    $seq =~ s/^>.*//;
    $seq =~ s/\s+//g;

    return length($seq);
}

sub get_seq_fast {
    my ($r_gene, $tmp_seq, $dbh) = @_;

    open(TMP_SEQ, ">$tmp_seq") || die;
    for (my $i=0; $i<@{$r_gene}; $i++) {
        my $seq = get_seq_mysql($dbh, ${$r_gene}[$i]);
        if ($seq) {
            print TMP_SEQ $seq;
        }
    }
    close(TMP_SEQ);
}

sub get_geneset_seq {
    my ($geneset, $r_seq) = @_;

    my $seq_fused = "";
    my @genes_fused = ();
    for my $gene (split(/\|/, $geneset)) {
        my $seq = get_gene_seq($gene, $r_seq);
        if ($seq) {
            push @genes_fused, $gene;
            $seq_fused .= $seq;
        }
    }
    $seq_fused =~ s/^>.*//gm;
    $seq_fused =~ s/^\s*\n//gm;

    if ($seq_fused) {
        my $genes_fused = join("|", @genes_fused);
        return ">$genes_fused\n$seq_fused\n";
    } else {
        return "";
    }
}

sub get_gene_seq {
    my ($gene, $r_seq) = @_;

    my ($sp, $name) = decompose_gene_id($gene);

    my $seq = ${$r_seq}{$sp}{$gene};
    if ($seq) {
        return $seq;
    } else {
        print STDERR "Warning: no seq for $gene, thus drop the seq.\n";
        return "";
    }
}

sub get_annotation {
    my ($annotation_file, $r_annotation) = @_;

    if (-f $annotation_file) {
        open(ANNOT, $annotation_file) || die;
        while (<ANNOT>) {
            chomp;
            my ($key, $annotation) = split("\t", $_);
            if ($key and $annotation) {
                if (${$r_annotation}{$key}) {
                    print STDERR "Warning: duplicated annotation for $key\n";
                    ${$r_annotation}{$key} .= "; $annotation";
                }
                ${$r_annotation}{$key} = $annotation;
            }
        }
        close(ANNOT);
    } else {
        # print STDERR "Warning: cannot read annotation file $annotation_file\n";
    }
}

sub md5_of_file_contents {
    my ($file) = @_;

    return md5_hex(`cat $file`);
}

sub md5_of_fasta {
    my ($fasta) = @_;

    my $sequences = `cat $fasta`;
    # $sequences =~ s/^>.*/>/gm;
    # $sequences =~ s/[-\s]//g;
    my $md5_value = md5_hex($sequences);
    
    return $md5_value;
}

sub get_gene_idx {
    my ($r_gene, $r_gene_idx) = @_;
    
    for (my $i=0; $i<@{$r_gene}; $i++) {
        my $gene = ${$r_gene}[$i];
        ${$r_gene_idx}{$gene} = $i;
    }
}

sub read_fasta_entries {
    my ($fasta_file, $r_gene) = @_;
    
    my @line = `cat $fasta_file | grep '^>'`;
    chomp(@line);
    for my $line (@line) {
        if ($line =~ /^>(\S+)/) {
            my $gene = $1;
            push @{$r_gene}, $gene;
        } else {
            die;
        }
    }
}

sub decompose_gene_id {
    my ($gene_id) = @_;

    if (! defined $gene_id) {
        die;
    }
    
    my ($sp, $name);
    if ($gene_id =~ /^[A-Za-z0-9]+:/) {
        my @x = split(":", $gene_id);
        if (@x == 2) {
            ($sp, $name) = @x;
        } elsif (@x == 3) {
            ($sp, $name) = @x[0,2]; # bug ?
        } else {
            die;
        }
    } elsif ($gene_id =~ /^(\d+)\.(\S+)/) { # eggNOG
        ($sp, $name) = ($1, $2);
    } elsif ($gene_id =~ /^(\w+?)_(\S+)/) {
        ($sp, $name) = ($1, $2);
    } else {
        die $gene_id;
    }

    return ($sp, $name);
}

sub extract_region {
    my ($tmp_dclst_file, $geneset) = @_;

    my @start = ();
    my @end = ();
    open(EXTRACT_REGION_DCLST, "$tmp_dclst_file") || die;
    while (<EXTRACT_REGION_DCLST>) {
        chomp;
        my ($cluster, $gene, @domain_info) = split;
        unless (@domain_info and @domain_info % 3 == 0) {
            die $_;
        }
        if ($gene eq $geneset) {
            for (my $i=0; $i<@domain_info; $i+=3) {
                my ($domain, $begin_pos, $end_pos) = ($domain_info[$i], $domain_info[$i+1], $domain_info[$i+2]);
                push @start, $begin_pos;
                push @end, $end_pos;
            }
        }
    }
    close(EXTRACT_REGION_DCLST);
    my $start = min(@start);
    my $end = max(@end);
    if ($start > $end) {
        die;
    }

    return ($start, $end);
}

sub trim_to_region {
    my ($seq, $start, $end) = @_;

    if ($seq =~ /^(>.*?)\n/) {
        my $header = $1;
        $seq =~ s/^.*?\n//;
        $seq =~ s/\s//g;
        my $seq_len = length($seq);
        # if ($seq_len < $start || $seq_len < $end) {
        # print STDERR "$start-$end\n$header\n$seq";
        # }
        if ($seq_len < $end) {
            print STDERR "WARNING: $header $start-$end len=$seq_len\n";
        }
        if ($seq_len < $start) {
            print STDERR "ERROR: $header $start-$end len=$seq_len\n";
        }
        my $before_seq = substr($seq, 0, $start-1);
        my $region_seq = substr($seq, $start-1, $end-$start+1);
        my $after_seq = substr($seq, $end);
        return ($header, $region_seq, $before_seq, $after_seq);
    } else {
        die;
    }
}

sub get_dbh {
    my ($database) = @_;

    if (! defined $database) {
        $database = "mbgd";
    }

    use DBI;
    my $dbh = DBI->connect("DBI:mysql:$database:localhost;mysql_socket=/tmp/mysql.sock", "chiba", "chiba", {'RaiseError' => 1});

    return $dbh;
}

sub define_tmp_file {
    my ($prefix) = @_;
    
    my $tmp_file = "$ENV{DOMREFINE_TMP}/$ENV{HOSTNAME}.$$.$prefix";
    
    return $tmp_file;
}

sub remove_tmp_file {
    my ($file) = @_;
    
    if ($file and -f $file) {
        unlink $file;
    }
}

sub cache_file_path {
    my ($cache_dir, $md5_string, $suffix) = @_;

    my $sub_dir = substr($md5_string, 0, 2);

    my $cache_file = "$cache_dir/$sub_dir/${md5_string}${suffix}";

    return $cache_file;
}

sub save_stdin {
    my ($file) = @_;

    my $contents;
    if (defined $file) {
        system "cat > $file";
        $contents = `cat $file`;
    } else {
        $contents = `cat`;
    }

    return $contents;
}

sub save_contents {
    my ($contents, $file) = @_;
    
    open(FILE, ">$file") || die;
    print FILE $contents;
    close(FILE);
}

sub get_geneset_descr {
    my ($geneset, %opt) = @_;
    my $r_gene_descr = $opt{r_gene_descr};

    my @gene = split(/\|/, $geneset);

    my @descr = ();
    for (my $i=0; $i<@gene; $i++) {
        if ($opt{mysql}) {
            $descr[$i] = mysql_gene_descr($gene[$i]);
        } elsif (%{$r_gene_descr}) {
            my $r_gene_descr = $opt{r_gene_descr};
            if (${$r_gene_descr}{$gene[$i]}) {
                $descr[$i] = ${$r_gene_descr}{$gene[$i]};
            } elsif ($gene[$i] =~ /^([A-Za-z0-9]+):(\S+)$/) {
                my ($sp_lc, $name_uc) = (lc($1), uc($2));
                if (${$r_gene_descr}{"$sp_lc:$name_uc"}) {
                    $descr[$i] = ${$r_gene_descr}{"$sp_lc:$name_uc"};
                } else {
                    $descr[$i] = "NOT_FOUND";
                }
            } else {
                die;
            }
        }
    }

    return join(" | ", @descr);
}

sub get_seq_mysql {
    my ($dbh, $gene) = @_;

    my ($sp, $name) = decompose_gene_id($gene);

    my $r_r_id = $dbh->selectall_arrayref("select aaseq from gene where sp='$sp' and name='$name'");
    if (@{$r_r_id} != 1 or @{${$r_r_id}[0]} != 1) {
        die;
    }
    my $id = ${$r_r_id}[0][0];

    my $r_r_seq = $dbh->selectall_arrayref("select seq from proteinseq where id=$id");
    if (@{$r_r_seq} != 1 or @{${$r_r_seq}[0]} != 1) {
        die;
    }
    my $seq = ${$r_r_seq}[0][0];

    $seq =~ s/(.{1,60})/$1\n/g;
    return ">$gene\n$seq\n";
}

sub read_homology_data {
    my ($tabout, $r_hom) = @_;

    open(TABOUT, $tabout) || die;
    while (<TABOUT>) {
        my @f = split;
        if (@f != 8) {
            next;
        }
        my ($gene1, $gene2, $start1, $end1, $start2, $end2, $pam, $score) = @f;
        if (defined $score) {
            put_homology_in_hash($r_hom, $gene1, $gene2, $start1, $end1, $start2, $end2, $pam, $score);
        } else {
            print STDERR;
        }
    }
    close(TABOUT);
}

sub put_homology_in_hash {
    my ($r_hom, $gene1, $gene2, $start1, $end1, $start2, $end2, $pam, $score) = @_;

    # if (defined ${$r_hom}{$gene1}{$gene2}) {
    if (defined ${$r_hom}{"$gene1 $gene2"}) {
        die;
    }

    # ${$r_hom}{$gene1}{$gene2}{start1} = $start1;
    # ${$r_hom}{$gene1}{$gene2}{end1} = $end1;
    # ${$r_hom}{$gene1}{$gene2}{start2} = $start2;
    # ${$r_hom}{$gene1}{$gene2}{end2} = $end2;
    # ${$r_hom}{$gene1}{$gene2}{pam} = $pam;
    # ${$r_hom}{$gene1}{$gene2}{score} = $score;

    # ${$r_hom}{$gene1}{$gene2} = "start1:$start1 end1:$end1 start2:$start2 end2:$end2";
    ${$r_hom}{"$gene1 $gene2"} = "s1:$start1 e1:$end1 s2:$start2 e2:$end2";
}

################################################################################
### Private functions ##########################################################
################################################################################

sub extract_dclst_with_check {
    my ($dclst_table_file, @cluster) = @_;

    my $dclst = "";
    for my $cluster (@cluster) {
        my $dclst_part = extract_dclst($dclst_table_file, $cluster);
        if ($dclst_part eq "") {
            print STDERR "[$cluster] is not containded.\n";
        }
        $dclst .= $dclst_part;
    }
    
    return $dclst;
}

sub mysql_gene_descr {
    my ($gene) = @_;

    my ($sp, $name) = decompose_gene_id($gene);

    my $dbh = get_dbh();

    my $r_r_descr = $dbh->selectall_arrayref("select descr from gene where sp='$sp' and name='$name'");
    if (@{$r_r_descr} != 1 or @{${$r_r_descr}[0]} != 1) {
        print STDERR "Cannot find annotation for $gene\n";
        # die;
    }

    my $descr = ${$r_r_descr}[0][0];
    if (! $descr) {
        $descr = "NOT_FOUND";
    }

    return $descr;
}

1;
