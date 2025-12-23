#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM -S clusterScore cluster.domrefine.tree.o11
-c CLUSTER_FILE
-h HOMCLUSTER_CLUSTER
-s HOMCLUSTER_SCORE
";

use DomRefine::Read;

### Settings ###
my %OPT;
getopts('c:h:s:S:', \%OPT);

if (@ARGV != 1) {
    print STDERR $USAGE;
    exit 1;
}
my ($TREE_DIR) =  @ARGV;
if (! -d $TREE_DIR) {
    die $USAGE;
}

my $CLUSTER_FILE;
if ($TREE_DIR =~ /^(\S+).tree.o11$/) {
    $CLUSTER_FILE = $1;
}
if ($OPT{c}) {
    $CLUSTER_FILE = $1;
}
unless ($CLUSTER_FILE and -s $CLUSTER_FILE) {
    die $USAGE;
}

my $SCORE_FILE = $OPT{S} || die $USAGE;

my %HOMCLUSTER_SCORE = ();
my $HOMCLUSTER_ID_NEXT;
if ($OPT{s} and -f $OPT{s}) {
    get_annotation($OPT{s}, \%HOMCLUSTER_SCORE);

    my $homcluster_id_max = 0;
    open(HOMCLUSTER_SCORE, $OPT{s}) || die;
    while (<HOMCLUSTER_SCORE>) {
        my ($homcluster_id) = split;
        if ($homcluster_id > $homcluster_id_max) {
            $homcluster_id_max = $homcluster_id;
        }
    }
    close(HOMCLUSTER_SCORE);
    $HOMCLUSTER_ID_NEXT = $homcluster_id_max + 1;
}

### Main ###
my %CLUSTER = ();
my %DOMAIN = ();
my %HOMCLUSTER_CLUSTER = ();

get_dclst_structure($CLUSTER_FILE, \%CLUSTER, \%DOMAIN);

### Print ###
print "# simtype=0\n";
print "# cutoff=300.000000\n";
print "# missdist=320.000000\n";
print "# phylocutratio=0.500000\n";
print "# distscale=40\n";

my $NODE_COUNT = 0;
my %NODE_ID = ();

open(SCORE, ">$SCORE_FILE") || die;
if ($OPT{h} and -s $OPT{h}) {
    my %printed_cluster = ();
    open(HOMCLUSTER_CLUSTER, "$OPT{h}") || die;
    while (<HOMCLUSTER_CLUSTER>) {
        my ($homcluster, $cluster) = split;
        if ($HOMCLUSTER_CLUSTER{$homcluster}) {
            push @{$HOMCLUSTER_CLUSTER{$homcluster}}, $cluster;
        } else {
            $HOMCLUSTER_CLUSTER{$homcluster} = [$cluster];
        }
        $printed_cluster{$cluster} = 1;
    }
    close(HOMCLUSTER_CLUSTER);

    for my $homcluster (sort {$a <=> $b} keys %HOMCLUSTER_CLUSTER) {
        my @out = ();
        my @score_txt = ();
        for my $cluster (sort {$a <=> $b} @{$HOMCLUSTER_CLUSTER{$homcluster}}) {
            output_cluster_tree($TREE_DIR, $cluster, \@out, \@score_txt);
        }
        print_homcluster(\@out, $homcluster, $HOMCLUSTER_SCORE{$homcluster});
        if (@score_txt) {
            print SCORE join("\n", @score_txt), "\n";
        }
    }

    my @out = ();
    my @score_txt = ();
    for my $cluster (sort {$a cmp $b} keys %CLUSTER) {
        if ($printed_cluster{$cluster}) {
            next;
        }
        output_cluster_tree($TREE_DIR, $cluster, \@out, \@score_txt);
    }
    print_homcluster(\@out, $HOMCLUSTER_ID_NEXT++);
    if (@score_txt) {
        print SCORE join("\n", @score_txt), "\n";
    }
} else {
    for my $cluster (sort {$a cmp $b} keys %CLUSTER) {
        my @out = ();
        my @score_txt = ();
        output_cluster_tree($TREE_DIR, $cluster, \@out, \@score_txt);
        if (@out) {
            print @out;
        }
        if (@score_txt) {
            print SCORE join("\n", @score_txt), "\n";
        }
    }
}
close(SCORE);

################################################################################
### Functions ##################################################################
################################################################################
sub print_homcluster {
    my ($r_out, $homcluster, $homcluster_score) = @_;

    if (@{$r_out}) {
        if (defined $homcluster_score) {
            print "HomCluster $homcluster : $homcluster_score\n";
        } else {
            print "HomCluster $homcluster : 0.0 0.0\n";
            print STDERR "WARNING: no score for HomCluster $homcluster\n";
        }
        print @{$r_out};
    }
}

sub output_cluster_tree {
    my ($tree_dir, $cluster, $r_out, $r_score_txt) = @_;

    if (! -s "$tree_dir/$cluster.out") {
        print STDERR "WARNING: no output for Cluster $cluster\n";
        return;
    }

    my $tree = "";
    my $score_txt;
    open(CLUSTER_TREE, "$tree_dir/$cluster.out") || die;
    my @cluster_tree = <CLUSTER_TREE>;
    close(CLUSTER_TREE);
    for (my $i=0; $i<@cluster_tree; $i++) {
        if ($cluster_tree[$i] =~ /Cluster \S:(score=\S+:dist=\S+)/) {
            $score_txt = $1;
        } elsif ($cluster_tree[$i] =~ /Cluster .*/) {
        } elsif ($cluster_tree[$i] =~ /^([IL]) (\S+) (\S+) (.*)/) {
            # $tree .= $cluster_tree[$i];
            my ($node_type, $node_id, $parent_id, $node_info) = ($1, $2, $3, $4);
            save_node_id($cluster, $node_id);
            save_node_id($cluster, $parent_id);
            # $tree .= "$node_type $node_id $parent_id $node_info\n";
            $tree .= "$node_type $NODE_ID{$cluster}{$node_id} $NODE_ID{$cluster}{$parent_id} $node_info\n";
        } elsif (/^#/) {
        } else {
            print STDERR "WARNING: Cluster $cluster: $cluster_tree[$i]";
        }
    }

    if ($tree eq "") {
        print STDERR "ERROR: no tree for Cluster $cluster\n";
        return;
    }

    if ($score_txt) {
        # $score_txt = "Cluster $cluster:$score_txt";
        push @{$r_score_txt}, "Cluster $cluster:$score_txt"
    }

    push @{$r_out}, "Cluster $cluster\n$tree\n//\n";
    # return ("Cluster $cluster\n$tree\n//\n", $score_txt);
}

sub save_node_id {
    my ($cluster, $node_id) = @_;

    if (defined $NODE_ID{$cluster}{$node_id}) {
        return;
    }

    if ($node_id == "-1") {
        $NODE_ID{$cluster}{$node_id} = "-1";
        return;
    }

    $NODE_COUNT ++;
    $NODE_ID{$cluster}{$node_id} = $NODE_COUNT;
}
