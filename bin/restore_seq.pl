#!/usr/bin/perl -s

use List::Util qw(max min);
use File::Path qw(make_path);

$domclust_file = $ARGV[0];
$preclust_file = $ARGV[1];
$seq_file = $ARGV[2];

if ($domclust_file =~ /\.o11$/) {
	$fmt='o11';
}

if (! $domclust_file || ! $preclust_file) {
	die "Usage: $0 domclust_file preclust_file [seq_file]\n";
}
if ($seqout && ! $seq_file) {
	die "Usage: $0 domclust_file preclust_file [seq_file]\n"
		. "seq_file is required when -seqout is specieid\n";
}

# read redundant sequences to be added
print STDERR "Reading preclust file\n";
open(F, $preclust_file) || die;
while(<F>){
	chomp;
	if (/^\* *(\S+)/) {
		$repr_seq = $1;
	} elsif (/^ /) {
		s/^ +//;
		($addseq, $repr_beg, $repr_end, $beg, $end, $ident) = split(/\t/);
		push(@{$AddSeq{$repr_seq}},
			{seq=>$addseq, tgt_begin=>$beg, tgt_end=>$end, repr_begin=>$repr_beg, repr_end=>$repr_end, ident=>$ident});
	}
###	last if (++$cnt > 100000);
}
close(F);
print STDERR "Done\n";

if ($seqout) {
	print STDERR "Reading sequence file\n";
	open(F, $seq_file) || die;
	while(<F>) {
		chomp;
		if (/>\s*(\S+)/) {
			$name = $1;
		} else {
			$Seq{$name} .= $_;
		}
	}
	print STDERR "Done\n";
	if ($seqout eq '1') {
		## defualt output directory
		$seqoutDir = "addseqs";
	} else {
		$seqoutDir = $seqout;
	}
	make_path( $seqoutDir );
}

print STDERR "Reading domclust file\n";
&read_domclust_file($domclust_file);

print STDERR "Assigning domain numbers to the added sequences\n";
foreach $gene (keys %GeneDom) {
	## assign domain ID for added sequence
	if (@{$GeneDom{$gene}} > 1) {
		my($domi);
		my(@posSortedList) = sort {$a->{begin}<=>$b->{begin}} @{$GeneDom{$gene}};
		foreach $dom (@posSortedList) {
#print "$gene, $dom->{begin}, $dom->{end}\n";
			$dom->{dom} = ++$domi;
		}
		$GeneDom{$gene} = \@posSortedList;
	}
}
print STDERR "Output results\n";
if ($fmt eq 'o11') {
	&outputAllClusters_o11;
} else {
	&outputAllClusters;
}

sub read_domclust_file {
	my($domclust_file) = @_;
	my($test_cnt);

	open(F, $domclust_file) || die;
	while(<F>){
		my($read_member_flag);
		my($name, $begin, $end);
		if (/^HomCluster (\S+)/) {
			$homclustid = $1;
			push(@AllHomClusters, $homclustid);
#			print;
		} elsif (/^Cluster (\S+)/) {
			$clustid = $1;
#			print;
			if ($homclustid) {
				push(@{$HomClust{$homclustid}}, $clustid);
			} else {
				push(@AllClusters, $clustid);
			}
#last if (++$test_cnt > 1000);
		} elsif ($fmt eq 'o11' && /^L (\d+) (\d+)/) {
			($dmy,$dmy1,$dmy2,$name, $begin, $end) = split(/ /);
			$read_member_flag = 1;
		} elsif ($fmt ne 'o11' && /^\w/) {
			chomp;
			($name, $begin, $end) = split(/ /);
			$read_member_flag = 1;
		}
		if ($read_member_flag) {
			my($dom);
#			print "$_\n";
			if ($name =~ /\((\d+)\)/) {
				$dom = $1;
				$name =~ s/\(\d+\)//;
			}
			$mem = {name=>$name, begin=>$begin, end=>$end, dom=>$dom};
			push(@{$Members{$clustid}}, $mem);
			if ($AddSeq{$name}) {
				# %adsInfo: information of the sequences to be added
				foreach $adsqInfo (@{$AddSeq{$name}}) {
					@ov = &ovlp($begin, $end, $adsqInfo->{repr_begin}, $adsqInfo->{repr_end});
					if (@ov) {
						my($obegin, $oend) = @ov;
#print ">$adsqInfo->{repr_begin} $adsqInfo->{repr_end}; $adsqInfo->{tgt_begin} $adsqInfo->{tgt_end}\n";
#print "ov>($obegin, $oend)\n";
						my($adsq_begin) = &conv_coord($obegin, $adsqInfo);
						my($adsq_end) = &conv_coord($oend, $adsqInfo);
						my($mem) = {name=>$adsqInfo->{seq}, begin=>$adsq_begin, end=>$adsq_end, add=>1};
						if ($fmt eq 'o11') {
							push(@{$AddGeneInfo{$clustid, $name}}, $mem);
						} else {
							push(@{$Members{$clustid}}, $mem);
#							print "+$adsqInfo->{seq} $adsq_begin $adsq_end\n";
						}
						push(@{$GeneDom{$adsqInfo->{seq}}}, $mem);
					} else {
						print "#### $adsqInfo->{seq} $adsqInfo->{repr_begin},$adsqInfo->{repr_end}\n" if ($verbose);
					}
				}
			}
		}
	}
}

sub outputAllClusters {
	if (@AllHomClusters) {
		foreach $homClust (@AllHomClusters) {
			print "HomCluster $homClust\n";
			foreach $clust (@{$HomClust{$homClust}}) {
				&outputCluster($clust);
			}
		}
	} else {
		foreach $clust (@AllClusters) {
			&outputCluster($clust);
			if ($seqout) {
				&outputSequence($clust);
			}
		}
	}
}
sub outputCluster {
	my($clust) = @_;
	print "Cluster $clust\n";
	foreach $mem (sort {$a->{name} cmp $b->{name}} @{$Members{$clust}}) {
		my($name) = $mem->{name};
		$name .= "($mem->{dom})" if ($mem->{dom});
		if ($DEBUG) {
			$name = "+ $name" if ($mem->{add});
		}
		print join(" ", $name, $mem->{begin}, $mem->{end}), "\n";
	}
	print "\n";
}
sub outputAllClusters_o11 {
	open(F, $domclust_file) || die;
	while(<F>){
		if (/^Cluster (\S+)/) {
			$clustid = $1;
			print;
		} elsif (/^L /) {
			print;
			chomp;
			my($stat, $nodenum, $parentnum, $name, $from, $to) = split(/ /);
			foreach $mem (@{$AddGeneInfo{$clustid, $name}}) {
				print join(" ", 'L2', $nodenum, $parentnum, $mem->{name}, $mem->{begin}, $mem->{end}, $mem->{dom}), "\n";
			}
		} else {
			print;
		}
	}
}
sub outputSequence {
	my($clust) = @_;
	open(O, ">$seqoutDir/$clust.add.fas");
	foreach my $mem (sort {$a->{name} cmp $b->{name}} @{$Members{$clust}}) {
		my($seq) = $Seq{$mem->{name}};
		if ($mem->{dom}) {
			$seq = substr($seq, $mem->{begin}, $mem->{end} - $mem->{begin} + 1);
		}
		if ($mem->{add}) {
			print O ">$mem->{name}\n";
			print O "$seq\n";
		}
	}
	close(O);
}

sub ovlp {
	my($begin1, $end1, $begin2, $end2) = @_;
	my($minend) = min($end1, $end2);
	my($maxbegin) = max($begin1, $begin2);
	my($ovlen) = ($minend - $maxbegin);
	my($len1) = $end1 - $begin1 + 1;
	my($len2) = $end2 - $begin2 + 1;
#print ">>$ovlen; $len1, $len2\n";
	if ($ovlen > min($len1, $len2) * 0.15) {
		return ($maxbegin, $minend);
	}
	return ();
}

sub conv_coord {
	my($coord1, $adsq) = @_;
	my($ali1_begin, $ali1_end, $ali2_begin, $ali2_end) = ($adsq->{repr_begin} - 1, $adsq->{repr_end}, $adsq->{tgt_begin} - 1, $adsq->{tgt_end});
	my($ret_coord2);
	## conv_tp_0_base
	$ret_coord2 = int($ali2_begin + ($coord1 - $ali1_begin) / ($ali1_end - $ali1_begin) * ($ali2_end - $ali2_begin) +.5);
	$ret_coord2;
}
