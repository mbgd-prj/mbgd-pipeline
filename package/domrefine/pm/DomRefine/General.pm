package DomRefine::General;
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(max max_i log2 initialize_matrix sum_sq sum match overlap_len uniq
	     min min_i mean
	     check_redundancy add_group_member
             check_queue diff_time convert_time get_time get_time_obj
             seq_file_stat
	     );

use strict;
# use Time::Piece;

sub log2 {
    my ($x) = @_;

    return (log($x)/log(2));
}

sub sum_sq {
    my @x = @_;
    my $sum = 0;
    for (my $i=0; $i<@x; $i++) {
	$sum += $x[$i]**2;
    }
    return $sum;
}

sub max {
    my @x = @_;
    my $max = $x[0];
    for (my $i=1; $i<@x; $i++) {
	if ($x[$i] > $max) {
	    $max = $x[$i];
	}
    }
    return $max;
}

sub max_i {
    my @x = @_;
    my $max = $x[0];
    my $max_i = 0;
    for (my $i=1; $i<@x; $i++) {
	if ($x[$i] > $max) {
	    $max = $x[$i];
	    $max_i = $i;
	}
    }
    return $max_i;
}

sub min {
    my @x = @_;
    my $min = $x[0];
    for (my $i=1; $i<@x; $i++) {
	if ($x[$i] < $min) {
	    $min = $x[$i];
	}
    }
    return $min;
}

sub min_i {
    my @x = @_;
    my $min;
    my $min_i;
    for (my $i=1; $i<@x; $i++) {
	if (! defined $x[$i]) {
	    next;
	}
	if (! defined $min) {
	    $min = $x[$i];
	    $min_i = $i;
	}
	if ($x[$i] < $min) {
	    $min = $x[$i];
	    $min_i = $i;
	}
    }
    return $min_i;
}

sub sum {
    my @x = @_;
    my $sum = 0;
    for (my $i=0; $i<@x; $i++) {
	$sum += $x[$i];
    }
    return $sum;
}

sub mean {
    my @x = @_;
    my $mean;
    if (@x != 0) {
	$mean = sum(@x)/@x;
    }
    return $mean;
}

sub uniq {
    my @x = @_;

    my %hash = ();
    for my $x (@x) {
	$hash{$x} = 1;
    }

    return sort {$a cmp $b} keys(%hash);
}

sub check_redundancy {
    my @x = @_;

    my %check = ();
    my @redundant = ();
    for my $x (@x) {
	if ($check{$x}) {
	    push @redundant, $x;
	} else {
	    $check{$x} = 1;
	}
    }

    return uniq(@redundant);
}

sub add_group_member {
    my ($r_hash, $group, @cluster) = @_;
    
    for my $cluster (@cluster) {
	if (${$r_hash}{$group}) {
	    push @{${$r_hash}{$group}}, $cluster;
	} else {
	    ${$r_hash}{$group} = [$cluster];
	}
    }
}

sub paste {
    my (@r_a) = @_;

    my @out = ();
    my $m = @r_a;
    my $n = @{$r_a[0]};
    for (my $j=0; $j<$m; $j++) {
	if (@{$r_a[$j]} > $n) {
	    $n = @{$r_a[$j]};
	}
    }
    
    for (my $i=0; $i<$n; $i++) {
	my @tmp = ();
	for (my $j=0; $j<$m; $j++) {
	    push @tmp, ${$r_a[$j]}[$i];
	}
	$out[$i] = join(" ", @tmp);
    }

    return @out;
}

sub deviation {
    my ($begin, $end, $begin2, $end2) = @_;

    my $deviation = abs($begin2 - $begin) + abs($end2 - $end);

    return $deviation/2;
}

sub seq_length {
    my ($seq) = @_;

    $seq =~ s/^>.*//m;
    $seq =~ s/\s//g;

    return length($seq);
}

sub seq_file_stat {
    my ($seq_file) = @_;

    my $n_seq = 0;
    my $n_aa = 0;
    open(SEQ_FILE_STAT, $seq_file) || die;
    while (my $line = <SEQ_FILE_STAT>) {
	if ($line =~ /^>/) {
	    $n_seq ++;
	} else {
	    $line =~ s/\s//g;
	    $n_aa += length($line);
	}
    }
    close(SEQ_FILE_STAT);
    my $mean_len;
    if ($n_seq) {
	$mean_len = $n_aa / $n_seq;
    } else {
	print STDERR "ERROR: no sequence in seq_file $seq_file\n";
    }

    return ($n_seq, $n_aa, $mean_len);
}

sub overlap_len {
    my ($begin, $end, $begin2, $end2) = @_;

    my $overlap = min($end, $end2) - max($begin, $begin2) + 1;
    if ($overlap < 0) {
	$overlap = 0;
    }

    return $overlap;
}

sub initialize_matrix {
    my ($r_d, $n, $m, $val) = @_;
    
    for (my $i=0; $i<$n; $i++) {
	for (my $j=0; $j<$m; $j++) {
	    ${$r_d}[$i][$j] = $val;
	}
    }
}

sub match {
    my ($query, @db) = @_;
    
    for my $db (@db) {
	if ($query eq $db) {
	    return 1;
	}
    }

    return 0;
}

sub check_queue {
    my ($n_jobs, $interval) = @_;

    while (1) {
	my @jobs;
	if ($ENV{DOMREFINE_QUEUE}) {
	    @jobs = `qstat | grep -Pw '$ENV{DOMREFINE_QUEUE}|qw'`;
	} else {
	    @jobs = `qstat | grep -w '$ENV{USER}'`;
	}
	if (@jobs + $ENV{DOMREFINE_QSUB_UNIT} > $ENV{DOMREFINE_QUEUE_SIZE}) {
	    sleep $ENV{DOMREFINE_QUEUE_CHECK_INTERVAL};
	} else {
	    last;
	}
    }
}

sub diff_time {
    my ($start_time, $end_time) = @_;

    my $sec = convert_time($end_time) - convert_time($start_time);
    if ($sec < 0) {
        $sec += 24 * 3600
    }

    return $sec;
}

sub convert_time {
    my ($time) = @_;

    if ($time =~ /^(\d+):(\d+)[:\.](\d+)$/) {
        my ($h, $m, $s) = ($1, $2, $3);
        my $sec = $h * 3600 + $m * 60 + $s;
        return $sec;
    } else {
        die;
    }
}

sub get_time {
    # my $time = `date '+%T'`;
    # chomp($time);
    # return $time;

    my ($s, $m, $h) = localtime;
    return "$h:$m.$s";
}

sub get_time_obj {
    my $time = localtime; # If Time::Piece is used, it overwrites the default localtime().
}

1;
