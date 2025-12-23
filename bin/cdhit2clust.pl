#!/usr/bin/perl
while(<>){
	if (/>Cluster (\d+)/) {
		if (defined $clustid) {
			&output_seqs($clustid, $repr, \@memseqs, \%aliinfo);;
		}
		$clustid = $1;
		undef @memseqs;
		undef %aliinfo;
		$repr = '';
	} elsif (/^(\d+)\t(\d+)aa, >(\S+)\.\.\. (.*)$/) {
		$seqn = $1;
		$aalen = $2;
		$seqname = $3;
		$hitinfo = $4;
		if ($hitinfo eq '*') {
			$repr = $seqname;
		} else {
			if ($hitinfo =~ /at (\d+):(\d+):(\d+):(\d+)\/([\d\.]+)%/) {
				$aliinfo{$seqname} = {begin1=>$1, end1=>$2, begin2=>$3, end2=>$4, identity=>$5};
			}
			push(@memseqs, $seqname);
		}
	}
}
if (defined $clustid) {
	&output_seqs($clustid, $repr, \@memseqs, \%aliinfo);;
}

sub output_seqs {
	my($clustid, $repr, $memseqs, $aliinfo) = @_;
	print "* $repr\n";
	foreach $sqname (@memseqs) {
		print "  $sqname";
		if ($aliinfo->{$sqname}) {
			$info = $aliinfo->{$sqname};
			print "\t" . join("\t", $info->{begin2}, $info->{end2}, $info->{begin1}, $info->{end1}, $info->{identity});
		}
		print "\n";
	}
	print "\n";
}
