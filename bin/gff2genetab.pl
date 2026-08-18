#!/usr/bin/perl -s
#
use File::Basename;
$default_shape = 'linear' if (! $default_shape);

if ($seqfile) {
	foreach $f (split(/,/, $seqfile)) {
		my($name, $sp);
		open(F, $f) || warn "Can't open sequence file: $f\n";
		while(<F>) {
			if(/^>\s*(\S+)/) {
				$name = $1;
				if ($name =~ /:/) {
					($sp,$name) = split(/:/, $name);
				}
				$Length{$name} = 0;
			} else {
				chomp;
				$Length{$name} += length($_);
			}
		}
		close(F);
	}
}
while(<>){
	($seqname, $source, $feat, $from, $to, $score, $strand, $frame, $attrStr) = split(/\t/);
	if ($ARGV ne  $prevFile) {
		$spname= basename(dirname($ARGV));
		if ($lower_spname) {
			$spname =~ s/\..*$//;
			$spname = lc($spname);
		} else {
			$spname =~ s/\.[a-z]+$//;
		}
		print "##Genome\tsp:$spname\n";
		$prevFile = $ARGV;
		$prev_seqname = '';
		$fasta_flag = 0;
		$prevFile = $ARGV;
	}
	if (/^#/) {
		if (/^#!/) {
		} elsif (/sequence-region\s+(\S+)\s+(\d+)\s+(\d+)/) {
			$name = $1; $from = $2; $to = $3;
			$seqlen = $to - $from + 1;
			$seqLen{$name} = $seqlen;
		} elsif (/^##FASTA/) {
			$fasta_flag = 1;
		}
		next;
	} elsif ($fasta_flag) {
		next;
	} else {
		$attr = &processAttr($attrStr);
	}
	if ($seqname ne $prev_seqname) {
		if ($prev_seqname) {
			&print_endmark($shape);
		}
		print "##Chromosome";
		print "\tname:$seqname" if ($seqname);
		print "\tseq_length:$seqLen{$seqname}" if ($seqLen{$seqname});
		print "\n";
	}
	if ($feat eq 'region') {
		if ($attr->{Is_circular} eq 'true') {
			$shape = 'circular';
		} elsif ($attr->{Is_circular} eq 'false') {
			$shape = 'linear';
		} elsif ($attr->{Is_circular} eq '') {
			$shape = $default_shape;
		}
	} elsif ($feat eq 'CDS') {
		$name = &get_genename($attr, $tag);
		if ($attr->{'pseudo'} eq 'true') {
			next;
		}

		$dir = ($strand eq '+') ? 1 : -1;
		if ($name =~ /${spname}:/) {
			$name =~ s/${spname}://;
		} elsif ($name =~ /:/) {
			$name =~ s/:/_/;
		}

		if ($Length{$name}) {
			$aalen = $Length{$name};
		} elsif (! $attr->{partial}) {
			$aalen = int( (($to - 3) - $from + 1) / 3 ); ## containing a stop codon
		} else {
			$aalen = int( ($to - $from + 1) / 3 ); ## no stop codon
		}
		$pos = ($from + ($to - 3)) / 2;
		print join("\t", $spname, $name, $aalen, $pos, $dir), "\n";
	}
	$prev_seqname = $seqname;
}
if ($prev_seqname) {
	&print_endmark($shape);
}


sub processAttr {
	my($str) = @_;
	my(%Attr);
	foreach my $attr (split(/\s*;\s*/, $str)) {
		($varname, $val) = split(/=/, $attr);
		$Attr{$varname} = $val;
	}
	\%Attr;
}

sub get_genename {
	my($attr, $tag) = @_;
	my($name);
	if ($tag) {
		$name = $attr->{$tag};
	} elsif ($attr->{locus_tag}) {
		$name = $attr->{locus_tag};
	} elsif ($attr->{ID}) {
		$name = $attr->{ID};
	} elsif ($attr->{protein_id}) {
		$name = $attr->{protein_id};
	}
	if ($name =~ /:/) {
		my($spname, $genename) = split(/:/, $name);
		if ($lower_spname) {
			$name = join(':', lc($spname), $genename);
		}
		$name = join(':', $spname, $genename);
	}
	return($name);
}
sub print_endmark {
	my($shape) = @_;
	if ($shape eq '') {
		$shape = 'default_shape';
	}
	if ($shape eq 'linear') {
		print "#1\n";
	} elsif ($shape eq 'circular') {
		print "#2\n";
	}
}
