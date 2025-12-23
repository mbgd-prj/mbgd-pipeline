#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM
";

my %OPT;
getopts('', \%OPT);

my $ALIGNER = $ENV{DOMREFINE_ALIGNER};

my $START_TIME = time;

### merge ###
my $CLUSTER = "cluster";
if (-s 'cluster.checked.link.to_check') {
    $CLUSTER = "cluster.checked";
}
if (-s "$CLUSTER.link.to_check") {
    my $total = extract_execution_time_with_file_list("$CLUSTER.link.to_check", "$CLUSTER.merge_test/", '.log', '^merge_test: \d+ sec$', 1);
    my $aligner = extract_execution_time_with_file_list("$CLUSTER.link.to_check", "$CLUSTER.merge_test/", '.log', '^'.$ALIGNER.':.* sec$', -2);
    my $dsp_score = extract_execution_time_with_file_list("$CLUSTER.link.to_check", "$CLUSTER.merge_test/", '.log', '^dsp_score: .* sec$', 3);
    my $others = $total - $aligner - $dsp_score;
    printf "merge_test\t$total\t%.2f\n", $total/60;
    printf "$ALIGNER\t$aligner\t%.2f\n", $aligner/60;
    printf "dsp_score\t$dsp_score\t%.2f\n", $dsp_score/60;
    printf "others(Perl)\t$others\t%.2f\n", $others/60;
    print "\n";
}

### merge_divide ###
if (-s 'cluster.merge.link') {
    my $total = extract_execution_time_with_file_list('cluster.merge.link', 'cluster.merge.merge_divide_test/score.before.', '.log', '^dom_score_c: \d+ sec$', 1);
    my $aligner = extract_execution_time_with_file_list('cluster.merge.link', 'cluster.merge.merge_divide_test/score.before.', '.log', '^'.$ALIGNER.':.* sec$', -2);
    my $dsp_score = extract_execution_time_with_file_list('cluster.merge.link', 'cluster.merge.merge_divide_test/score.before.', '.log', '^dsp_score: .* sec$', 3);
    my $others = $total - $aligner - $dsp_score;
    printf "score_before\t$total\t%.2f\n", $total/60;
    printf "$ALIGNER\t$aligner\t%.2f\n", $aligner/60;
    printf "dsp_score\t$dsp_score\t%.2f\n", $dsp_score/60;
    printf "others(Perl)\t$others\t%.2f\n", $others/60;
    print "\n";
}
if (-s 'cluster.merge.link') {
    my $total = extract_execution_time_with_file_list('cluster.merge.link', 'cluster.merge.merge_divide_test/', '.log', '^merge_divide: \d+ sec$', 1);
    my $aligner = extract_execution_time_with_file_list('cluster.merge.link', 'cluster.merge.merge_divide_test/', '.log', '^'.$ALIGNER.':.* sec$', -2);
    my $FastTree = extract_execution_time_with_file_list('cluster.merge.link', 'cluster.merge.merge_divide_test/', '.log', '^FastTree: .* sec$', 3);
    my $others = $total - $aligner - $FastTree;
    printf "merge_divide\t$total\t%.2f\n", $total/60;
    printf "$ALIGNER\t$aligner\t%.2f\n", $aligner/60;
    printf "FastTree\t$FastTree\t%.2f\n", $FastTree/60;
    printf "others(Perl)\t$others\t%.2f\n", $others/60;
    print "\n";
}
if (-s 'cluster.merge.link') {
    my $total = extract_execution_time_with_file_list('cluster.merge.link', 'cluster.merge.merge_divide_test/score.after.', '.log', '^dom_score_c: \d+ sec$', 1);
    my $aligner = extract_execution_time_with_file_list('cluster.merge.link', 'cluster.merge.merge_divide_test/score.after.', '.log', '^'.$ALIGNER.':.* sec$', -2);
    my $dsp_score = extract_execution_time_with_file_list('cluster.merge.link', 'cluster.merge.merge_divide_test/score.after.', '.log', '^dsp_score: .* sec$', 3);
    my $others = $total - $aligner - $dsp_score;
    printf "score_after\t$total\t%.2f\n", $total/60;
    printf "$ALIGNER\t$aligner\t%.2f\n", $aligner/60;
    printf "dsp_score\t$dsp_score\t%.2f\n", $dsp_score/60;
    printf "others(Perl)\t$others\t%.2f\n", $others/60;
    print "\n";
}

### move boundary ###
if (-s 'cluster.merge.merge_divide.link') {
    my $total = extract_execution_time_with_file_list('cluster.merge.merge_divide.link', 'cluster.merge.merge_divide.move_test/score.before.', '.log', '^dom_score_c: \d+ sec$', 1);
    my $aligner = extract_execution_time_with_file_list('cluster.merge.merge_divide.link', 'cluster.merge.merge_divide.move_test/score.before.', '.log', '^'.$ALIGNER.':.* sec$', -2);
    my $dsp_score = extract_execution_time_with_file_list('cluster.merge.merge_divide.link', 'cluster.merge.merge_divide.move_test/score.before.', '.log', '^dsp_score: .* sec$', 3);
    my $others = $total - $aligner - $dsp_score;
    printf "score_before\t$total\t%.2f\n", $total/60;
    printf "$ALIGNER\t$aligner\t%.2f\n", $aligner/60;
    printf "dsp_score\t$dsp_score\t%.2f\n", $dsp_score/60;
    printf "others(Perl)\t$others\t%.2f\n", $others/60;
    print "\n";
}
if (-s 'cluster.merge.merge_divide.link') {
    my $total = extract_execution_time_with_file_list('cluster.merge.merge_divide.link', 'cluster.merge.merge_divide.move_test/', '.log', '^move_boundary: \d+ sec$', 1);
    my $aligner = extract_execution_time_with_file_list('cluster.merge.merge_divide.link', 'cluster.merge.merge_divide.move_test/', '.log', '^'.$ALIGNER.':.* sec$', -2);
    my $dsp_score = extract_execution_time_with_file_list('cluster.merge.merge_divide.link', 'cluster.merge.merge_divide.move_test/', '.log', '^dsp_score: .* sec$', 3);
    my $dsp_score2 = extract_execution_time_with_file_list('cluster.merge.merge_divide.link', 'cluster.merge.merge_divide.move_test/', '.log', '^\d+ sec$', 0);
    my $others = $total - $aligner - $dsp_score - $dsp_score2;
    printf "move_boundary\t$total\t%.2f\n", $total/60;
    printf "$ALIGNER\t$aligner\t%.2f\n", $aligner/60;
    printf "dsp_score\t$dsp_score\t%.2f\n", $dsp_score/60;
    printf "dsp_score2\t$dsp_score2\t%.2f\n", $dsp_score2/60;
    printf "others(Perl)\t$others\t%.2f\n", $others/60;
    print "\n";
}
if (-s 'cluster.merge.merge_divide.link') {
    my $total = extract_execution_time_with_file_list('cluster.merge.merge_divide.link', 'cluster.merge.merge_divide.move_test/score.after.', '.log', '^dom_score_c: \d+ sec$', 1);
    my $aligner = extract_execution_time_with_file_list('cluster.merge.merge_divide.link', 'cluster.merge.merge_divide.move_test/score.after.', '.log', '^'.$ALIGNER.':.* sec$', -2) || 0;
    my $dsp_score = extract_execution_time_with_file_list('cluster.merge.merge_divide.link', 'cluster.merge.merge_divide.move_test/score.after.', '.log', '^dsp_score: .* sec$', 3);
    my $others = $total - $aligner - $dsp_score;
    printf "score_after\t$total\t%.2f\n", $total/60;
    printf "$ALIGNER\t$aligner\t%.2f\n", $aligner/60;
    printf "dsp_score\t$dsp_score\t%.2f\n", $dsp_score/60;
    printf "others(Perl)\t$others\t%.2f\n", $others/60;
    print "\n";
}

### create boundary ###
my $MERGE_DIVIDE = "cluster.merge.merge_divide";
if (-s "${MERGE_DIVIDE}.link") {
    my $total = extract_execution_time_with_file_list("${MERGE_DIVIDE}.link", "${MERGE_DIVIDE}.move.create_test/", '.log', '^create_boundary: \d+ sec$', 1);
    my $aligner = extract_execution_time_with_file_list("${MERGE_DIVIDE}.link", "${MERGE_DIVIDE}.move.create_test/", '.log', '^'.$ALIGNER.':.* sec$', -2) || 0;
    my $dsp_score = extract_execution_time_with_file_list("${MERGE_DIVIDE}.link", "${MERGE_DIVIDE}.move.create_test/", '.log', '^dsp_score: .* sec$', 3);
    my $others = $total - $aligner - $dsp_score;
    printf "create_boundary\t$total\t%.2f\n", $total/60;
    printf "$ALIGNER\t$aligner\t%.2f\n", $aligner/60;
    printf "dsp_score\t$dsp_score\t%.2f\n", $dsp_score/60;
    printf "others(Perl)\t$others\t%.2f\n", $others/60;
    print "\n";
}
if (-s "${MERGE_DIVIDE}.link") {
    my $total = extract_execution_time_with_file_list("${MERGE_DIVIDE}.link", "${MERGE_DIVIDE}.move.create_test/score.after.", '.log', '^dom_score_c: \d+ sec$', 1);
    my $aligner = extract_execution_time_with_file_list("${MERGE_DIVIDE}.link", "${MERGE_DIVIDE}.move.create_test/score.after.", '.log', '^'.$ALIGNER.':.* sec$', -2);
    my $dsp_score = extract_execution_time_with_file_list("${MERGE_DIVIDE}.link", "${MERGE_DIVIDE}.move.create_test/score.after.", '.log', '^dsp_score: .* sec$', 3);
    my $others = $total - $aligner - $dsp_score;
    printf "score_after\t$total\t%.2f\n", $total/60;
    printf "$ALIGNER\t$aligner\t%.2f\n", $aligner/60;
    printf "dsp_score\t$dsp_score\t%.2f\n", $dsp_score/60;
    printf "others(Perl)\t$others\t%.2f\n", $others/60;
    print "\n";
}

### divide ###
my $DIVIDE_TEST = "cluster.merge.merge_divide.move.create.paste.divide_test";
if (-s "${DIVIDE_TEST}.to_check") {
    my $total = extract_execution_time_with_file_list("${DIVIDE_TEST}.to_check", "${DIVIDE_TEST}/", '.log', '^divide: \d+ sec$', 1);
    my $aligner = extract_execution_time_with_file_list("${DIVIDE_TEST}.to_check", "${DIVIDE_TEST}/", '.log', '^'.$ALIGNER.':.* sec$', -2);
    my $FastTree = extract_execution_time_with_file_list("${DIVIDE_TEST}.to_check", "${DIVIDE_TEST}/", '.log', '^FastTree: .* sec$', 3);
    my $others = $total - $aligner - $FastTree;
    printf "divide_test\t$total\t%.2f\n", $total/60;
    printf "$ALIGNER\t$aligner\t%.2f\n", $aligner/60;
    printf "FastTree\t$FastTree\t%.2f\n", $FastTree/60;
    printf "others(Perl)\t$others\t%.2f\n", $others/60;
    print "\n";
}

### alignment ###
if (-s 'cluster.domrefine.to_check') {
    my $total = extract_execution_time_with_file_list('cluster.domrefine.to_check', 'cluster.domrefine.alignment/', '.log', '^dom_align: \d+ sec$', 1);
    my $aligner = extract_execution_time_with_file_list('cluster.domrefine.to_check', 'cluster.domrefine.alignment/', '.log', '^'.$ALIGNER.':.* sec$', -2);
    my $others = $total - $aligner;
    printf "alignment\t$total\t%.2f\n", $total/60;
    printf "$ALIGNER\t$aligner\t%.2f\n", $aligner/60;
    printf "others(Perl)\t$others\t%.2f\n", $others/60;
    print "\n";
}

### tree ###
if (-s 'cluster.domrefine.to_check') {
    my $total = extract_execution_time_with_file_list('cluster.domrefine.to_check', 'cluster.domrefine.tree/', '.log', '^dom_tree_score: \d+ sec$', 1);
    my $FastTree = extract_execution_time_with_file_list('cluster.domrefine.to_check', 'cluster.domrefine.tree/', '.log', '^FastTree: .* sec$', 3);
    my $dsp_score = extract_execution_time_with_file_list('cluster.domrefine.to_check', 'cluster.domrefine.tree/', '.log', '^dsp_score: .* sec$', 3);
    my $others = $total - $FastTree - $dsp_score;
    printf "tree\t$total\t%.2f\n", $total/60;
    printf "FastTree\t$FastTree\t%.2f\n", $FastTree/60;
    printf "dsp_score\t$dsp_score\t%.2f\n", $dsp_score/60;
    printf "others(Perl)\t$others\t%.2f\n", $others/60;
    print "\n";
}

my $END_TIME = time;
printf STDERR "summarize time:\t%.2f\tmin\n", ($END_TIME - $START_TIME)/60;

################################################################################
### Function ###################################################################
################################################################################

sub extract_execution_time_with_file_list {
    my ($list, $prefix, $suffix, $pattern, $ind) = @_;

    my @list = `cat $list`;
    chomp(@list);
    
    my $total = 0;
    for my $name (@list) {
	my $file = "$prefix$name$suffix";
	if (-s $file) {
	    open(FILE, $file) || die;
	    while (<FILE>) {
		chomp;
		if (/$pattern/) {
		    my @x = split;
		    my $sec = $x[$ind];
		    if ($sec =~ /^\d+$/) {
			$total += $sec;
		    } else {
			print STDERR "WARNING: $file '$_' $pattern [$ind]";
		    }
		}
	    }
	    close(FILE);
	} elsif (-z $file) {
	    # print STDERR "$file is empty.\n";
	} else {
	    # print STDERR "$file does not exist.\n";
	}
    }

    return $total;
}
