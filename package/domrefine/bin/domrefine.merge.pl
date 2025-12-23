#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM [OPTION] CLUSTER_FILE [CLUSTESET_FILE]
-r: recursive merge
-t THRESHOLD
-e GAP_EXT
-o GAP_OPEN

Oversplit domains are merged iteratively using the DSP score change.
Exceptionally, the largest cluster is not merged.

Intermediate output files:
CLUSTER_FILE*.merge_test(|.jobs|.jobs.watch|.summary)
CLUSTER_FILE*.merge(|.log|.renumber)
CLUSTER_FILE*.link.(to_check|to_merge)
CLUSTER_FILE*.clusterset.merged
CLUSTER_FILE.to_ignore

Final output files:
CLUSTER_FILE.merge
CLUSTER_FILE.clusterset.merged.total
";

print STDERR "\n";
print STDERR '@' . `date`;
print STDERR ">$PROGRAM @ARGV\n";

my %OPT;
getopts('rt:e:o:q:', \%OPT);

if ($OPT{r}) {
    $ENV{DOMREFINE_RECURSIVE_MERGE} = 1;
}

### Settings ###
my $THRESHOLD = -0.05;
if (defined $OPT{t}) {
    $THRESHOLD = $OPT{t};
}

my $MERGE_TEST = "dom_merge_test_c";
if($ENV{DOMREFINE_ALIGN_REGION}){
   $MERGE_TEST .= " -r"; # align by region
}
if (defined $OPT{e}) {
    $MERGE_TEST .= " -e $OPT{e}";
}
if (defined $OPT{o}) {
    $MERGE_TEST .= " -o $OPT{o}";
}

my $QUEUE = "";
if (defined $OPT{q}) {
    $QUEUE = "-q $OPT{q}";
}

my $N_IGNORE = 1;

### Main ###
if (@ARGV == 0) {
    print STDERR $USAGE;
    exit 1;
}
my ($INPUT_CLUSTER, $CLUSTERSET_FILE) = @ARGV;

if ($ENV{DOMREFINE_RECURSIVE_MERGE}) {
    my ($prefix, $n);
    if ($INPUT_CLUSTER =~ /^(\S*\D).merge(\d+)$/) {
	($prefix, $n) = ($1, $2);
    } elsif ($INPUT_CLUSTER =~ /^(\S*\D)(\d+)$/) {
	($prefix, $n) = ($1, $2);
    } else {
	($prefix, $n) = ($INPUT_CLUSTER, 0);
    }
    my $tmp_cluster = $INPUT_CLUSTER;
    do {
	($CLUSTERSET_FILE) = merge($tmp_cluster, $CLUSTERSET_FILE);
	$n ++ ;
	if (! -e "$prefix.merge$n") {
	    system "ln -fs ${tmp_cluster}.merge.out.renumber $prefix.merge$n"; # overwrite when resuming
	}
	$tmp_cluster = "$prefix.merge$n";
    } until (-z $CLUSTERSET_FILE);
    system "ln -fs $tmp_cluster $prefix.merge";
    system "cat *.link.to_merge | links_to_clustersets.pl -n | sort -k1,1nr > $prefix.clusterset.merged.total";
} else {
    merge($INPUT_CLUSTER, $CLUSTERSET_FILE);
    system "ln -fs ${INPUT_CLUSTER}.merge.out.renumber ${INPUT_CLUSTER}.merge";
}
    
################################################################################
### Function ###################################################################
################################################################################
sub merge {
    my (${tmp_cluster}, $clusterset) = @_;

    my $outdir = "${tmp_cluster}.merge_test";

    my $start_time = time;
    if (! -s "${tmp_cluster}.size") {
	system "cat $tmp_cluster | dom_cluster_size > ${tmp_cluster}.size";
    }
    # restrict links
    if (! defined $clusterset) {
	$clusterset = "";
    }
    my @cluster_to_ignore = `head -n ${N_IGNORE} ${tmp_cluster}.size | cut -f1`;
    chomp(@cluster_to_ignore);
    add_cluster_to_ignore("cluster.to_ignore", @cluster_to_ignore);
    if (! -s "${tmp_cluster}.link") {
	system "cat $tmp_cluster | dom_network -l > ${tmp_cluster}.link";
    }
    if (! -s "${tmp_cluster}.link.to_check") {
	system "cat ${tmp_cluster}.link | ignore_clusters.pl cluster.to_ignore | restrict_links.pl $clusterset > ${tmp_cluster}.link.to_check";
    }
    my $end_time = time;
    printf STDERR "pre:\t%.2f\tmin\n", ($end_time - $start_time)/60;

    # dom_merge_test
    if (! $ENV{DOMREFINE_SKIP_MERGE_TEST}) {
	system "cat ${tmp_cluster}.link.to_check | shuffle.pl | sge_script.pl $QUEUE -i ${tmp_cluster} -o $outdir '$MERGE_TEST'";
    }

    $start_time = time;

    # summarize
    if (! $ENV{DOMREFINE_SKIP_MERGE_TEST}) {
	my @link = `cat ${tmp_cluster}.link.to_check`;
	chomp(@link);
	open(OUTDIR_SUMMARY, ">$outdir.summary") || die;
	for my $link (@link) {
	    if (-f "$outdir/$link.out" and ! -z "$outdir/$link.out") {
		my $out = `cat $outdir/$link.out`;
		print OUTDIR_SUMMARY "$link\t$out";
	    }
	}
	close(OUTDIR_SUMMARY);
    }

    # execute
    if (! -s "${tmp_cluster}.link.to_merge") {
	system "cat $outdir.summary | perl -lane '\$F[1]>=$THRESHOLD and \$F[7]>0 and print' | sort -k2,2gr | cut -f1 > ${tmp_cluster}.link.to_merge";
    }
    if ($ENV{DOMREFINE_FAST_MERGE}) {
	system "cat ${tmp_cluster}.link.to_merge | links_to_clustersets.pl -n | sort -k1,1nr > ${tmp_cluster}.clusterset.merged";
	system "cat ${tmp_cluster}.clusterset.merged | cut -f2 | clusterset_merge -i ${tmp_cluster} > ${tmp_cluster}.merge.out 2> ${tmp_cluster}.merge.log";
    } else {
	# system "cat ${tmp_cluster} | dom_merge_stepwise ${tmp_cluster}.link.to_merge > ${tmp_cluster}.merge.out 2> ${tmp_cluster}.merge.log";
	system "cat ${tmp_cluster} | cut.sh 2 | cut.sh 1 : | sort | uniq > ${tmp_cluster}.sp";
	my @sp = `cat ${tmp_cluster}.sp`;
	chomp(@sp);
	my $sp_list = join(",", @sp);
	my $pwd = `pwd`;
	chomp($pwd);
	my $large_queue = $ENV{DOMREFINE_LARGE_QUEUE} || "smpl";
	my $stepwise_merge_command = "dom_merge_stepwise -r 0 -R 0.2 -l $pwd/${tmp_cluster}.link.to_merge $pwd/${tmp_cluster}";
	my $command;
	if ($ENV{DOMREFINE_HOMOLOGY_INFO}) {
	    my $homology_info = $ENV{DOMREFINE_HOMOLOGY_INFO};
	    if ($homology_info !~ /^\//) {
		$homology_info = "$pwd/$homology_info";
	    }
            if ($large_queue eq "local") {
                system "cat $homology_info | $stepwise_merge_command > ${tmp_cluster}.merge.out 2> ${tmp_cluster}.merge.log";
            } else {
                $command = "sge.pl -N merge -q $large_queue 'cat $homology_info | $stepwise_merge_command > $pwd/${tmp_cluster}.merge.out 2> $pwd/${tmp_cluster}.merge.log'";
            }
	} elsif ($ENV{DOMREFINE_HOMOLOGY_DIR}) {
	    my $homology_dir = $ENV{DOMREFINE_HOMOLOGY_DIR};
	    if ($homology_dir !~ /^\//) {
		$homology_dir = "$pwd/$homology_dir";
	    }
	    my $select_command = "source /db/project/MBGD/etc/profile.mbgd && /db/project/MBGD/WWW/bin/select.pl -DIR=$homology_dir -SPEC=$sp_list -EVAL=0.001 -SCORE=60 -tabout";
            if ($large_queue eq "local") {
                system "$select_command 2> $pwd/${tmp_cluster}.tabout.log | $stepwise_merge_command > ${tmp_cluster}.merge.out 2> ${tmp_cluster}.merge.log";
            } else {
                $command = "sge.pl -N merge -q $large_queue '$select_command 2> $pwd/${tmp_cluster}.tabout.log | $stepwise_merge_command > $pwd/${tmp_cluster}.merge.out 2> $pwd/${tmp_cluster}.merge.log'";
            }
	} else {
	    die;
	}
        if ($large_queue ne "local") {
            open(MERGE_COMMAND, ">$pwd/${tmp_cluster}.merge.command") || die;
            print MERGE_COMMAND "$command\n";
            close(MERGE_COMMAND);
            my $submitted_job = `$command`;
            open(SUBMITTED_JOB, ">$pwd/${tmp_cluster}.merge.job") || die;
            print SUBMITTED_JOB $submitted_job;
            close(SUBMITTED_JOB);
            system "sge_check_jobs.pl -N -n0 $pwd/${tmp_cluster}.merge.job > $pwd/${tmp_cluster}.merge.job.check";
        }
	system "cat ${tmp_cluster}.merge.log | grep '^MERGED: ' | perl -pe 's/^MERGED: //' | links_to_clustersets.pl -n | sort -k1,1nr > ${tmp_cluster}.clusterset.merged";
    }
    system "cat ${tmp_cluster}.merge.out | dom_renumber | dom_renumber -c > ${tmp_cluster}.merge.out.renumber";

    $end_time = time;
    printf STDERR "post:\t%.2f\tmin\n", ($end_time - $start_time)/60;

    return ("${tmp_cluster}.clusterset.merged");
}

sub add_cluster_to_ignore {
    my ($file, @cluster_to_ignore) = @_;

    my %hash = ();
    if (-f $file) {
	my @cluster = `cat $file`;
	chomp(@cluster);
	for my $cluster (@cluster) {
	    $hash{$cluster} = 1;
	}
    }

    my @new_cluster = ();
    for my $cluster_to_ignore (@cluster_to_ignore) {
	if (! $hash{$cluster_to_ignore}) {
	    push @new_cluster, $cluster_to_ignore;
	}
    }

    open(CLUSTER_TO_IGNORE, ">>$file") || die;
    for my $new_cluster (@new_cluster) {
	print CLUSTER_TO_IGNORE "$new_cluster\n";
    }
    close($file);
}
