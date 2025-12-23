#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM [OPTION] CLUSTER_FILE 2>> domrefine.log &
-t THRESHOLD
-e GAP_EXT
-o GAP_OPEN
";

my %OPT;
getopts('t:e:o:', \%OPT);

### Parameters ###
my $T = "";
my $E = "";
my $O = "";
if ($OPT{t}) {
    $T = "-t $OPT{t}";
}
if ($OPT{e}) {
    $E = "-e $OPT{e}";
}
if ($OPT{o}) {
    $O = "-o $OPT{o}";
}

if (@ARGV != 1) {
    print STDERR $USAGE;
    exit 1;
}
my ($INPUT_FILE) = @ARGV;

my $START_TIME = time;

### Preparation ###
if ($ENV{DOMREFINE_BLASTDBCMD}){
    if (!-e "$ENV{DOMREFINE_SEQ_DB}.phd"){
        system "/bio/bin/makeblastdb -dbtype prot -hash_index -parse_seqids -in $ENV{DOMREFINE_SEQ_DB}";
    }
}

system "domrefine.prepare_cache.pl $ENV{DOMREFINE_CACHE}";
if (! $ENV{DOMREFINE_FAST_MERGE}) {
    system "domrefine.prepare_cache.pl $ENV{DOMREFINE_CACHE_STEPWISE_MERGE}";
}

if (! -e "cluster") {
    if ($INPUT_FILE eq "cluster") {
    } elsif ($INPUT_FILE =~ /cluster$/) {
        system "ln -s $INPUT_FILE cluster";
    } elsif ($INPUT_FILE =~ /\.o0$/) {
        system "cat $INPUT_FILE | domclust_to_tsv.pl > cluster";
    } elsif ($INPUT_FILE =~ /\.o11$/) {
        system "cat $INPUT_FILE | o11_to_files.pl";
    } elsif ($INPUT_FILE =~ /\.dclst$/) {
        system "cat $INPUT_FILE | domrefine.dclst_to_files.pl";
    } elsif ($INPUT_FILE =~ /\.dclst.gz$/) {
        system "cat $INPUT_FILE | gunzip | domrefine.dclst_to_files.pl";
    } else {
        system "cat $INPUT_FILE | domclust_to_tsv.pl > cluster";
    }
}
if (-s "cluster" and (! -s "cluster.checked" || ! -f "cluster.checked.log")) {
    system "cat cluster | dom_check -lt > cluster.checked 2> cluster.checked.log";
    system "diff_domtbl cluster cluster.checked > cluster.checked.diff";
}

### Main ###
if (! -s "cluster.merge") {
    my $start_time = time;
    if (-s "cluster.checked" and -s "cluster.checked.diff") {
        system "domrefine.merge.pl $T $E $O cluster.checked";
        system "ln -s cluster.checked.merge cluster.merge";
    } else {
        system "domrefine.merge.pl $T $E $O cluster";
    }
    my $end_time = time;
    printf STDERR "merge:\t%.2f\tmin\n", ($end_time - $start_time)/60;
}
if (-s "cluster.merge" and ! -s "cluster.merge.merge_divide") {
    my $start_time = time;
    system "domrefine.merge_divide.pl $E $O cluster.merge";
    my $end_time = time;
    printf STDERR "merge_divide:\t%.2f\tmin\n", ($end_time - $start_time)/60;
}
if (-s "cluster.merge.merge_divide" and ! -s "cluster.merge.merge_divide.move.create") {
    my $start_time = time;
    system "domrefine.boundary.pl $E $O cluster.merge.merge_divide";
    my $end_time = time;
    printf STDERR "boundary:\t%.2f\tmin\n", ($end_time - $start_time)/60;
}
if (-s "cluster.merge.merge_divide.move.create" and ! -s "cluster.merge.merge_divide.move.create.paste.divide") {
    my $start_time = time;
    if ($ENV{DOMREFINE_PRECLUST_INFO}) {
        system "cat cluster.merge.merge_divide.move.create | dom_renumber -d | dom_paste 2> cluster.merge.merge_divide.move.create.paste.log | dom_add_preclust.pl $ENV{DOMREFINE_PRECLUST_INFO} > cluster.merge.merge_divide.move.create.paste";
    } else {
        system "cat cluster.merge.merge_divide.move.create | dom_renumber -d | dom_paste > cluster.merge.merge_divide.move.create.paste 2>  cluster.merge.merge_divide.move.create.paste.log";
    }
    system "domrefine.divide.pl -1 cluster.merge.merge_divide.move.create.paste";
    my $end_time = time;
    printf STDERR "divide:\t%.2f\tmin\n", ($end_time - $start_time)/60;
}

if (-s "cluster.merge.merge_divide.move.create.paste.divide" and ! -s "cluster.domrefine") {
    my $start_time = time;
    system "cat cluster.merge.merge_divide.move.create.paste.divide | dom_renumber -d | dom_renumber -c > cluster.domrefine";
    my $end_time = time;
    print STDERR "\n";
    printf STDERR "formatting:\t%.2f\tmin\n", ($end_time - $start_time)/60;
}
if (-s "cluster.domrefine" and ! -s "cluster.domrefine.o0") {
    system "cat cluster.domrefine | tsv_to_domclust.pl > cluster.domrefine.o0";
}
if (-s "cluster.domrefine" and ! -s "cluster.domrefine.size") {
    system "cat cluster.domrefine | dom_cluster_size > cluster.domrefine.size";
}

### Post-process ###
print STDERR "\n";
print STDERR '@' . `date`;
if (-s "cluster.domrefine" and ($ENV{DOMREFINE_SAVE_ALIGNMENT} || $ENV{DOMREFINE_SAVE_TREE} || $ENV{DOMREFINE_SAVE_O11} || $ENV{DOMREFINE_SAVE_NEWICK})) {
    if (! -s "cluster.domrefine.to_check") {
        my $start_time = time;
        if ($ENV{DOMREFINE_IGNORE_LARGEST_CLUSTER}) {
            system "cat cluster.domrefine | cut.sh 1 | sort | uniq | ignore_clusters.pl cluster.to_ignore > cluster.domrefine.to_check";
        } else {
            system "cat cluster.domrefine | cut.sh 1 | sort | uniq > cluster.domrefine.to_check";
        }
        my $end_time = time;
        printf STDERR "pre:\t%.2f\tmin\n", ($end_time - $start_time)/60;
    }
    if (-s "homcluster" and ! -s "cluster.domrefine.homcluster_cluster") {
        system "cat cluster.domrefine | dom_overlap -h homcluster > cluster.domrefine.homcluster_cluster";
    }
    if ($ENV{DOMREFINE_SAVE_ALIGNMENT}) {
        my $start_time = time;
        system "cat cluster.domrefine.to_check | sge_script.pl -i cluster.domrefine -o cluster.domrefine.alignment 'dom_align -D'";
        my $end_time = time;
        printf STDERR "alignment:\t%.2f\tmin\n", ($end_time - $start_time)/60;
    }
    if ($ENV{DOMREFINE_SAVE_TREE}) {
        my $start_time = time;
        system "cat cluster.domrefine.to_check | sge_script.pl -i cluster.domrefine -o cluster.domrefine.tree 'dom_tree_score'";
        my $end_time = time;
        printf STDERR "tree:\t%.2f\tmin\n", ($end_time - $start_time)/60;
        if (! -s "cluster.domrefine.o1") {
            my $start_time = time;
            system "trees_to_o1.pl -h cluster.domrefine.homcluster_cluster cluster.domrefine.tree > cluster.domrefine.o1 2> cluster.domrefine.o1.err";
            my $end_time = time;
            printf STDERR "o0_to_o1:\t%.2f\tmin\n", ($end_time - $start_time)/60;
        }
    }
    if ($ENV{DOMREFINE_SAVE_O11}) {
        my $start_time = time;
        system "cat cluster.domrefine.to_check | sge_script.pl -i cluster.domrefine -o cluster.domrefine.tree.o11 'dom_tree_score -o 11'";
        my $end_time = time;
        printf STDERR "o11:\t%.2f\tmin\n", ($end_time - $start_time)/60;
        if (! -s "cluster.domrefine.o11") {
            my $start_time = time;
            system "trees_to_o11.pl -h cluster.domrefine.homcluster_cluster -s homcluster.score -S cluster.domrefine.clusterScore cluster.domrefine.tree.o11 > cluster.domrefine.o11 2> cluster.domrefine.o11.err";
            system "cat cluster.domrefine.o11 | grep '^Cluster ' | cut.sh 2 > cluster.domrefine.o11.list";
            my $end_time = time;
            printf STDERR "o0_to_o11:\t%.2f\tmin\n", ($end_time - $start_time)/60;
        }
    }
    if ($ENV{DOMREFINE_SAVE_NEWICK}) {
        my $start_time = time;
        system "cat cluster.domrefine.to_check | sge_script.pl -i cluster.domrefine -o cluster.domrefine.newick 'dom_tree_score -N'";
        my $end_time = time;
        printf STDERR "newick:\t%.2f\tmin\n", ($end_time - $start_time)/60;
    }
}
if ($ENV{DOMREFINE_SAVE_ARCHIVE}) {
    if (-s "cluster.domrefine.o11.list" and -s "cluster.domrefine.size") {
        my $archive_prefix = $ENV{DOMREFINE_SAVE_ARCHIVE};
        my $start_time = time;
        if (-d "cluster.domrefine.alignment" and ! -s "${archive_prefix}.alignment.tar.gz") {
            system "cat cluster.domrefine.o11.list | domrefine.archive.pl -a -s cluster.domrefine.size -n 2 cluster.domrefine.alignment ${archive_prefix}.alignment";
            system "tar zcf ${archive_prefix}.alignment.tar.gz ${archive_prefix}.alignment";
        }
        if (-d "cluster.domrefine.newick" and ! -s "${archive_prefix}.tree.tar.gz") {
            system "cat cluster.domrefine.o11.list | domrefine.archive.pl -N -s cluster.domrefine.size -n 2 cluster.domrefine.newick ${archive_prefix}.tree";
            system "tar zcf ${archive_prefix}.tree.tar.gz ${archive_prefix}.tree";
        }
        my $end_time = time;
        printf STDERR "archive:\t%.2f\tmin\n", ($end_time - $start_time)/60;
    }
}

my $END_TIME = time;
printf STDERR "\ndomrefine:\t%.2f\tmin\n", ($END_TIME - $START_TIME)/60;
