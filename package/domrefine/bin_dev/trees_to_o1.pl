#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM cluster.domrefine.tree
-c CLUSTER_FILE
-h HOMCLUSTER_CLUSTER
";

use DomRefine::Read;

### Settings ###
my %OPT;
getopts('c:h:', \%OPT);

if (@ARGV != 1) {
    print STDERR $USAGE;
    exit 1;
}
my ($TREE_DIR) =  @ARGV;
if (! -d $TREE_DIR) {
    die $USAGE;
}

my $CLUSTER_FILE;
if ($TREE_DIR =~ /^(\S+).tree$/) {
    $CLUSTER_FILE = $1;
}
if ($OPT{c}) {
    $CLUSTER_FILE = $1;
}
unless ($CLUSTER_FILE and -s $CLUSTER_FILE) {
    die $USAGE;
}

### Main ###
my %CLUSTER = ();
my %DOMAIN = ();
my %HOMCLUSTER_CLUSTER = ();

get_dclst_structure($CLUSTER_FILE, \%CLUSTER, \%DOMAIN);

if ($OPT{h} and -s $OPT{h}) {
    open(HOMCLUSTER_CLUSTER, "$OPT{h}") || die;
    while (<HOMCLUSTER_CLUSTER>) {
        my ($homcluster, $cluster) = split;
        if ($HOMCLUSTER_CLUSTER{$homcluster}) {
            push @{$HOMCLUSTER_CLUSTER{$homcluster}}, $cluster;
        } else {
            $HOMCLUSTER_CLUSTER{$homcluster} = [$cluster];
        }
    }
    close(HOMCLUSTER_CLUSTER);

    for my $homcluster (sort {$a <=> $b} keys %HOMCLUSTER_CLUSTER) {
        my @out = ();
        for my $cluster (sort {$a <=> $b} @{$HOMCLUSTER_CLUSTER{$homcluster}}) {
            my $out = output_cluster_tree($TREE_DIR, $cluster);
            if ($out) {
                push @out, $out;
            }
        }
        if (@out) {
            print "HomCluster $homcluster\n";
            print @out;
        }
    }
} else {
    for my $cluster (sort {$a cmp $b} keys %CLUSTER) {
        my $out = output_cluster_tree($TREE_DIR, $cluster);
        if ($out) {
            print $out;
        }
    }
}

################################################################################
### Functions ##################################################################
################################################################################
sub output_cluster_tree {
    my ($tree_dir, $cluster) = @_;

    if (! -s "$tree_dir/$cluster.out") {
        print STDERR "WARNING: no output for Cluster $cluster\n";
        return;
    }

    my ($score, $dist, $score_txt, $dist_txt);
    my $tree = "";
    open(CLUSTER_TREE, "$tree_dir/$cluster.out") || die;
    my @cluster_tree = <CLUSTER_TREE>;
    close(CLUSTER_TREE);
    for (my $i=0; $i<@cluster_tree; $i++) {
        if ($cluster_tree[$i] =~ /^#score\t(\S+)/) {
            $score = $1;
            $score_txt = "score=$score";
        } elsif ($cluster_tree[$i] =~ /^#dist\t(\S+)/) {
            $dist = $1;
            $dist_txt = "dist=$dist";
        } elsif ($cluster_tree[$i] =~ /^#/) {
        } else {
            $tree .= $cluster_tree[$i];
        }
    }

    if ($tree eq "") {
        print STDERR "ERROR: no tree for Cluster $cluster\n";
        return;
    }

    my $out = "Cluster $cluster";
    my @score_dist_txt = ();
    if ($score_txt) {
        push @score_dist_txt, $score_txt;
    }
    if ($dist_txt) {
        push @score_dist_txt, $dist_txt;
    }
    if (@score_dist_txt) {
        $out .= " [@score_dist_txt]";
    }
    $out .= "\n";
    $out .= "$tree";
    $out .= "\n";
    return $out;
}
