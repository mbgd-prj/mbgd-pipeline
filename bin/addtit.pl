#!/usr/bin/perl -s

$idcol1 = 0 if (! defined $idcol1);
$titcol1 = '$' if (! defined $titcol1);
$idcol2 = 1 if (! defined $idcol2);

$titfile = $ARGV[0];
$bloutfile = $ARGV[1];
die "Usage: $0 titfile bloutfile\n" if (! $titfile);
$MAXTITLEN = 120;

$bloutfile = "-" if (! $bloutfile);

if ($fasta) {
	open(F, $titfile) || die;
	while(<F>) {
		if (/^>/) {
			chomp;
			s/^>\s+//;
			
			$first_entry = substr($_, 1, $MAXTITLEN);
			($name, $title) = split(/\s+/, $first_entry, 2);
			$Title{$name} = $title;
		}
	}
	close(F);
} else {
	open(F, $titfile) || die;
	while(<F>) {
		chomp;
		@F = split(/\t/);

		$id = $F[$idcol1];
		if ($titcol1 eq '$') {
			$title = $F[$#F];
		} else {
			$title = $F[$titcol1];
		}
		$Title{$id} = $title;
	}
	close(F);
}

open(F, $bloutfile) || die "Can't open blastoutfile: $bloutfile\n";
while(<F>) {
	chomp;
	@F = split(/\t/);
	$id = $F[$idcol2];
	$id =~ s/\.\d+$//;
	print "$_\t$Title{$id}\n";
}
close(F);
