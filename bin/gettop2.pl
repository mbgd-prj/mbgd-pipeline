#!/usr/bin/perl -s

$MIN_OVLEN = 100;
$MIN_OVLP_RATIO = 0.2;

$EVAL_CUT = 0.001;

while (<>) {
	next if (/^#/);
	chomp;
#	my($qid, $sid, $ident, $alilen, $mismatch, $gap,
#		$qstart, $qend, $sstart, $send, $eval, $score, $title) = split;
	$hitseg = &readSeg($_);
	next if ($EVAL_CUT && $hitseg->{eval} > $EVAL_CUT);

	if ($skip_selfmatch) {
		$sid_orig = $hitseg->{sid};
		$sid_orig =~ s/\(\d+\)//;	## remove domain number
		next if ($hitseg->{qid} eq $sid_orig);
	} elsif ($skip_selfspmatch) {
		($qsp) = split(/:/, $hitseg->{qid});
		($ssp) = split(/:/, $hitseg->{sid});
		next if ($qsp eq $ssp);
	}

	if ($hitseg->{qid} ne $prev_qid) {
		&outputSegs(\@HitSegs);
		undef(@HitSegs);
		$prev_qid = $hitseg->{qid};
	}
	$hitseg->{from} = $hitseg->{qstart};
	$hitseg->{to} = $hitseg->{qend};

	if (! &checkOvlp($hitseg, \@HitSegs)) {
#print "ok:$hitseg->{qid}\n";
		push(@HitSegs, $hitseg);
	}
}
&outputSegs(\@HitSegs);

sub checkOvlp {
	my($hitseg, $HitSegs) = @_;
	my($flag) = 0;
	foreach $sg (@{$HitSegs}) {
		if (&ovlp($hitseg, $sg)) {
			$flag = 1; #last;
		}
	}
#if ($flag == 0 && @$HitSegs > 0) {print "******\n"};
	return $flag;
}

sub readSeg {
	my($str) = @_;
	my($qid, $sid, $ident, $alilen, $mismatch, $gap,
		$qstart, $qend, $sstart, $send, $eval, $score, $qlen, $slen, $title) = split(/\t/, $str);
	my($qdir) = '+1';
	my($sdir) = '+1';
	($qid) =~ s/ .*//;
	($sid) =~ s/ .*//;
	if ($qstart > $qend) {
		$tmp = $qstart; $qstart = $qend; $qend = $tmp;
		$qdir = -1;
	}
	if ($sstart > $send) {
		$tmp = $sstart; $sstart = $send; $send = $tmp;
		$sdir = -1;
	}
	{qid=>$qid, sid=>$sid, ident=>$ident, alilen=>$alilen,
		mismatch=>$mismatch, gap=>$gap, qstart=>$qstart, qend=>$qend, qdir=>$qdir,
		sstart=>$sstart, send=>$send, sdir=>$sdir, eval=>$eval, score=>$score, title=>$title};
}

sub ovlp {
	my($seg1, $seg2) = @_;
	my($ovlen) = &len( &max($seg1->{from}, $seg2->{from}), &min($seg1->{to}, $seg2->{to}) );
	my($len1) = &seglen($seg1);
	my($len2) = &seglen($seg2);
#print "$seg1->{from},$seg1->{to}; $seg2->{from},$seg2->{to}\n";
	return 1 if ($ovlen > $MIN_OVLEN || $ovlen > $MIN_OVLP_RATIO * &min($len1,$len2));
}

sub min {
	my($a,$b) = @_;
	$a < $b ? $a : $b;
}
sub max {
	my($a,$b) = @_;
	$a > $b ? $a : $b;
}
sub seglen {
	my($seg) = @_;
	&len($seg->{from}, $seg->{to});
}
sub len {
	my($from, $to) = @_;
	$to - $from + 1;
}
sub outputSegs {
	my($Segs) = @_;
	foreach $seg (@{$Segs}) {
		if ($outdir) {
			print join("\t", $seg->{qid}, $seg->{sid}, $seg->{ident}, $seg->{alilen},
					$seg->{mismatch}, $seg->{gap},
					$seg->{qstart}, $seg->{qend}, $seg->{qdir}, $seg->{sstart},
					$seg->{send}, $seg->{sdir}, $seg->{eval}, $seg->{score});
		} else {
			print join("\t", $seg->{qid}, $seg->{sid}, $seg->{ident}, $seg->{alilen},
					$seg->{mismatch}, $seg->{gap},
					$seg->{qstart}, $seg->{qend}, $seg->{sstart},
					$seg->{send}, $seg->{eval}, $seg->{score});
		}
		print "\t".$seg->{title} if ($seg->{title});
		print "\n";
	}
}
