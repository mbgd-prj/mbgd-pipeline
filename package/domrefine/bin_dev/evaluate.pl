#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM -d REF_CLUSTER [-i INITIAL_CLUSTER | -l FILE_LIST | FILE...]
-1: only 1 to 1 relationship
-g GENE_LIST: extract genes before evaluation
-r: resume
-s: save summary as a file (default: stdout)
";
# -w: well defined ortholog groups
# -o r_over
# -O r_over2
# -t r_over_t
# -T r_over_T

use DomRefine::General;

my %OPT;
getopts('d:i:l:1wg:ro:O:t:T:s', \%OPT);

my $PWD = `pwd`;
chomp($PWD);
$ENV{DOMREFINE_TMP} = $PWD;

### Input files ###
my @FILE = ();
if ($OPT{i}) {
    my $initial = $OPT{i};
    @FILE = ("$initial",
	     "$initial.merge",
	     "$initial.merge.merge_divide",
	     "$initial.merge.merge_divide.move",
	     "$initial.merge.merge_divide.move.create.paste",
	     "$initial.merge.merge_divide.move.create.paste.divide"
	);
} elsif ($OPT{l}) {
    @FILE = `cat $OPT{l}`;
    chomp(@FILE);
}
push @FILE, @ARGV;

my $REF_CLUSTER;
if ($OPT{d}) {
    $REF_CLUSTER = $OPT{d};
} else {
    die $USAGE;
}

### Get overlaps with reference ###
my $EXTRACT_GENE = "dom_renumber -d";
if ($OPT{g}) {
    $EXTRACT_GENE .= " | dom_extract_gene $OPT{g}";
}
if ($OPT{w}) {
    $EXTRACT_GENE .= " | dom_wdog";
}

my $OVERLAP = "dom_overlap";
if ($OPT{o}) {
    $OVERLAP .= " -o $OPT{o}";
}
if ($OPT{O}) {
    $OVERLAP .= " -O $OPT{O}";
}
if ($OPT{t}) {
    $OVERLAP .= " -t $OPT{t}";
}
if ($OPT{T}) {
    $OVERLAP .= " -T $OPT{T}";
}
if ($OPT{1}) {
    $OVERLAP .= " -1";
}

my $SUFFIX = "";
if ($OPT{w}) {
    $SUFFIX = "wdog";
}
if ($OPT{1}) {
    $SUFFIX .= "1to1";
} else {
    $SUFFIX .= "corr";
}
$SUFFIX .= "." . basename $REF_CLUSTER;

my @EVALUATED_FILE = ();
for my $file (@FILE) {
    if ($file !~ /^\s*#/ && $file !~ /\s/ && -f $file) {
	if ($OPT{r} && -f "$file.$SUFFIX") {
	    print STDERR "$file.$SUFFIX exists.\n";
	} else {
	    system "cat $file | $EXTRACT_GENE | $OVERLAP $REF_CLUSTER > $file.$SUFFIX";
	}
	push @EVALUATED_FILE, "$file.$SUFFIX";
    }
}

### Summary ###
if (@EVALUATED_FILE) {
    my @summary;
    if ($OPT{1}) {
	@summary = evaluation_summary_1to1(\@EVALUATED_FILE);
    } else {
	@summary = evaluation_summary_corr(\@EVALUATED_FILE);
    }
    if ($OPT{s}) {
	open(RESULT, ">evaluate.$SUFFIX") || die;
	print RESULT join("\n", @summary), "\n";
	close(RESULT);
    } else {
	print join("\n", @summary), "\n";
    }
}

################################################################################
### Function ###################################################################
################################################################################

sub evaluation_summary_corr {
    my ($r_evaluated_file) = @_;

    my @summary = ();

    push @summary, join("\t", "equiv", "sub", "super", "others", "either");
    for my $evaluated_file (@{$r_evaluated_file}) {
	my %type = ();
	my %either = ();
	open(FILE, $evaluated_file) || die;
	while (<FILE>) {
	    my @x = split;
	    my ($reference_cluster, $type) = ($x[3], $x[11]);
	    $type{$type}{$reference_cluster} ++;
	    if ($type ne "others") {
		$either{$reference_cluster}++;
	    }
	}
	close(FILE);
	my $n_equiav = keys %{$type{"equiv"}};
	my $n_sub = keys %{$type{"sub"}};
	my $n_super = keys %{$type{"super"}};
	my $n_others = keys %{$type{"others"}};
	my $n_either = keys %either;
	push @summary, join("\t", $n_equiav, $n_sub, $n_super, $n_others, $n_either, $evaluated_file);
    }    

    return @summary;
}

sub evaluation_summary_1to1 {
    my ($r_evaluated_file) = @_;

    my @num = ();
    for my $evaluated_file (@{$r_evaluated_file}) {
	my $result = `wc -l ${evaluated_file}`;
	chomp($result);
	if ($result =~ /^(\d+) (\S+)$/) {
	    my ($num, $file) = ($1, $2);
	    if ($file ne $evaluated_file) {
		die;
	    }
	    push @num, $num;
	} else {
	    die;
	}
    }

    my @summary = ();

    push @summary, join("\t", "Num", "Change", "Gain", "Loss", "Full List", "Gain List", "Loss List");
    push @summary, join("\t", $num[0], "NA", "NA", "NA", ${$r_evaluated_file}[0], "NA", "NA");
    for (my $i=1; $i<@num; $i++) {
	my $diff = $num[$i] - $num[$i-1];
	my @gain = ();
	my @loss = ();
	get_gain_loss(${$r_evaluated_file}[$i-1], ${$r_evaluated_file}[$i], \@gain, \@loss);
	my $gain = @gain;
	my $loss = @loss;
	push @summary, join("\t", $num[$i], $diff, $gain, -$loss, ${$r_evaluated_file}[$i], "'@gain'", "'@loss'");
	check_diff_gain_loss($diff, $gain, $loss);
    }
    my $diff_total = $num[-1] - $num[0];
    my @gain_total = ();
    my @loss_total = ();
    get_gain_loss(${$r_evaluated_file}[0], ${$r_evaluated_file}[-1], \@gain_total, \@loss_total);
    my $gain_total = @gain_total;
    my $loss_total = @loss_total;
    push @summary, join("\t", $num[-1], $diff_total, $gain_total, -$loss_total, "(Total Change)");
    check_diff_gain_loss($diff_total, $gain_total, $loss_total);

    return @summary;
}

sub check_diff_gain_loss {
    my ($diff, $gain, $loss) = @_;
    
    if ($diff == $gain - $loss) {
    } else {
	print STDERR "ERROR: $diff != $gain - $loss\n";
    }
}

sub get_gain_loss {
    my ($file1, $file2, $r_gain, $r_loss) = @_;

    my @motif1 = `cat $file1 | cut -f4`;
    chomp(@motif1);
    my %motif1;
    for my $motif (@motif1) {
	$motif1{$motif}++;
    }

    my @motif2 = `cat $file2 | cut -f4`;
    chomp(@motif2);
    my %motif2;
    for my $motif (@motif2) {
	$motif2{$motif}++;
    }
    
    my @motif = uniq(@motif1, @motif2);
    for my $motif (@motif) {
	if (! $motif1{$motif} and $motif2{$motif}) {
	    push @{$r_gain}, $motif;
	}
	if ($motif1{$motif} and ! $motif2{$motif}) {
	    push @{$r_loss}, $motif;
	}
    }
}
