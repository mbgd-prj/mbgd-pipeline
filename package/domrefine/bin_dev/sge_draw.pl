#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: cat CLUSTERSET | $PROGRAM -i FILE_OR_DIR -o DIR [OPTION] SCRIPT
-p PREFIX: prefix of output files in DIR
-q QUEUE (by default, small)
-r : resume

-m MOTIF_FILE_PREFIX
-a GENE_ANNOTATION
";

use DomRefine::General;

### Settings ###
my %OPT;
getopts('i:o:p:q:rm:a:1', \%OPT);

if (@ARGV != 1) {
    print STDERR $USAGE;
    exit 1;
}
my ($SCRIPT) = @ARGV;

my $DIR = $ENV{PWD};
my $OUT_DIR;
if ($OPT{o}) {
    $OUT_DIR = "$DIR/$OPT{o}";
} else {
    die $USAGE;
}
my $QUEUE = "";
if ($OPT{q}) {
    $QUEUE = "-q $OPT{q}";
}
my $PREFIX = "";
if ($OPT{p}) {
    $PREFIX = "$OPT{p}.";
}

if (! -e $OUT_DIR) {
    system "mkdir $OUT_DIR";
}

my $OPTION = "";
if ($OPT{m}) {
    $OPTION .= " -m $DIR/$OPT{m}";
}
if ($OPT{a}) {
    $OPTION .= " -a $DIR/$OPT{a}";
}

### Main ###
-t and die $USAGE;

my $SUBMITTED_JOBS = "$OUT_DIR.${PREFIX}jobs";
open(SUBMITTED_JOBS, ">$SUBMITTED_JOBS") || die;
my $count = 0;
while (my $cluster_set = <STDIN>) {
    $cluster_set =~ s/^(\S+?)\s.*/$1/;
    my $get_dclst;
    if ($OPT{i}) {
	my $input = "$DIR/$OPT{i}";
	if (-f $input) {
	    $get_dclst = "cat $input | dom_extract $cluster_set";
	} elsif (-d $input) {
	    if (! -f "$input/$cluster_set.out" or -z "$input/$cluster_set.out") {
		next;
	    }
	    $get_dclst = "cat $input/$cluster_set.out";
	} else {
	    die;
	}
    } else {
	die;
    }
    my $out_file_prefix = "$OUT_DIR/$PREFIX$cluster_set";
    my $out_file = "$out_file_prefix.png";
    if ($SCRIPT =~ /^dom_tree/) {
	$out_file = "$out_file_prefix.pdf";
    }
    if ($OPT{r}) {
	if (-f $out_file) {
	    next;
	}
    }
    if ($count % $ENV{DOMREFINE_QSUB_UNIT} == 0) {
	check_queue();
    }
    # my $submitted_job = `$get_dclst | $SCRIPT $OPTION $out_file_prefix`;
    my $submitted_job = `sge.pl $QUEUE -N d$cluster_set '$get_dclst | $SCRIPT $OPTION $out_file_prefix 2> $out_file_prefix.log'`;
    # my $submitted_job = `sge.pl $QUEUE -N d$cluster_set 'touch $out_file'`;
    print SUBMITTED_JOBS $submitted_job;
    $count ++;
}
close(SUBMITTED_JOBS);

system "sge_check_jobs.pl -N -n0 $SUBMITTED_JOBS > $SUBMITTED_JOBS.check";
