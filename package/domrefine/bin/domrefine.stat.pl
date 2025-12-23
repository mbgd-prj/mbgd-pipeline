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

print STDERR "\n";
print STDERR '@' . `date`;
print STDERR ">$PROGRAM\n";

my $START_TIME = time;

### Main ###
if (-s "cluster") {
    my $cluster = `cat cluster | cut.sh 1 | sort | uniq | wc -l`;
    chomp($cluster);
    my $cluster_split = `cat cluster | grep -v ' 0 ' | cut.sh 1 | sort | uniq | wc -l`;
    chomp($cluster_split);
    print "$cluster clusters\n";
    print "$cluster_split/$cluster =",format_percent($cluster_split/$cluster)," split\n";
}
if (-s "cluster.link.to_merge") {
    print "\nmerge\n";
    domrefine_stat(to_test => "cluster.link.to_check",
		   result => "cluster.merge_test.summary",
		   to_modify => "cluster.link.to_merge");
}
if (-s "cluster.checked.link.to_merge") {
    print "\nmerge\n";
    domrefine_stat(to_test => "cluster.checked.link.to_check",
		   result => "cluster.checked.merge_test.summary",
		   to_modify => "cluster.checked.link.to_merge");
}
if (-s "cluster.merge.merge_divide_test.to_divide") {
    print "\nmerge_divide\n";
    domrefine_stat(to_test => "cluster.merge.link",
		   result => "cluster.merge.merge_divide_test.summary",
		   to_modify => "cluster.merge.merge_divide_test.to_divide");
}
if (-s "cluster.merge.merge_divide.moved") {
    print "\nmove_boundary\n";
    domrefine_stat(to_test => "cluster.merge.merge_divide.link",
		   result => "cluster.merge.merge_divide.move_test.summary",
		   to_modify => "cluster.merge.merge_divide.move_test.to_move",
		   modified => "cluster.merge.merge_divide.moved");
}
if (-s "cluster.merge.merge_divide.move.create_test.to_create" and -s "cluster.merge.merge_divide.link") {
    print "\ncreate_boundary\n";
    domrefine_stat(to_test => "cluster.merge.merge_divide.link",
		   result => "cluster.merge.merge_divide.move.create_test.summary",
		   to_modify => "cluster.merge.merge_divide.move.create_test.to_create");
}
if (-s "cluster.merge.merge_divide.move.create.paste.divide_test.to_divide") {
    print "\ndivide\n";
    domrefine_stat(to_test => "cluster.merge.merge_divide.move.create | cut.sh 1 | sort | uniq",
		   result => "cluster.merge.merge_divide.move.create.paste.divide_test.summary",
		   valid_results => "cluster.merge.merge_divide.move.create.paste.divide_test.summary | grep .",
		   to_modify => "cluster.merge.merge_divide.move.create.paste.divide_test.to_divide");

    my $cluster = `cat cluster.merge.merge_divide.move.create.paste.divide | cut.sh 1 | sort | uniq | wc -l`;
    chomp($cluster);
    my $cluster_split = `cat cluster.merge.merge_divide.move.create.paste.divide | dom_renumber -d | dom_renumber | grep -v ' 0 ' | cut.sh 1 | sort | uniq | wc -l`;
    chomp($cluster_split);
    print "\n";
    print "$cluster clusters\n";
    printf "$cluster_split/$cluster =%.1f%% split\n", $cluster_split/$cluster*100;

}

my $END_TIME = time;
printf STDERR "summarize stats:\t%.2f\tmin\n", ($END_TIME - $START_TIME)/60;

################################################################################
### Functions ##################################################################
################################################################################

sub domrefine_stat {
    my (%opt) = @_;
    
    my $to_test;
    if ($opt{to_test}) {
	$to_test = `cat $opt{to_test} | ignore_clusters.pl cluster.to_ignore | wc -l`;
	chomp($to_test);
    } else {
	die;
    }

    unless ($to_test) {
	return;
    }

    if ($opt{result}) {
	my $result = `cat $opt{result} | wc -l`;
	chomp($result);
	print "$result/$to_test =", format_percent($result/$to_test), " tested\n";
    }
    if ($opt{valid_results}) {
	my $valid_results = `cat $opt{valid_results} | wc -l`;
	chomp($valid_results);
	print "$valid_results/$to_test =", format_percent($valid_results/$to_test), " obtained\n";
    }
    if ($opt{to_modify}) {
	my $to_modify = `cat $opt{to_modify} | wc -l`;
	chomp($to_modify);
	print "$to_modify/$to_test =", format_percent($to_modify/$to_test), " targeted\n";
    }
    if ($opt{modified}) {
	my $modified = `cat $opt{modified} | wc -l`;
	chomp($modified);
	print "$modified/$to_test =", format_percent($modified/$to_test), " modified\n";
    }
}

sub format_percent {
    my ($ratio) = @_;

    if ($ratio == 1) {
	return "100%";
    } else {
	return sprintf("%.1f%%", $ratio * 100);
    }
}
