#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM
";

my %OPT;
getopts('d:', \%OPT);

if ($OPT{d}) {
    
} else {
    !@ARGV and die $USAGE;
    for my $log_file (@ARGV) {
	analyze_log_file($log_file)
    }
}

################################################################################
### Functions ##################################################################
################################################################################
sub analyze_log_file {
    my ($log_file) = @_;
    
    my @log = `cat $log_file`;
    chomp(@log);

    for (my $i=0; $i<@log; $i++) {
	if ($log[$i] =~ /clustalo.*:\sn_seq=(\d+) n_aa=(\d+) mean_len=(\S+) (\d+) sec/) {
	    my ($n_seq, $n_aa, $mean_len, $align_sec) = ($1, $2, $3, $4);
	    if ($i+2 < @log and $log[$i+2] =~ /dsp_score:\sn_seq=(\d+) n_pos=(\d+) (\d+) sec/) {
		my ($n_seq2, $n_pos, $score_sec) = ($1, $2, $3, $4);
		if ($n_seq != $n_seq2) {
		    print STDERR "error\n";
		} else {
		    print "$n_seq\t$n_pos\t$mean_len\t$align_sec\t$score_sec\n";
		}
	    }
	}
    }
}
