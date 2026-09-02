#!/usr/bin/perl -s

use LWP::Simple;
use HTTP::Tiny;
use JSON;
use File::Path qw( make_path remove_tree );
use Time::HiRes qw( sleep );

$MBGDAPI = "https://mbgdapi.nibb.ac.jp";

@species = @ARGV;

if (! @species) {
	print STDERR "Usage: $0 [ -outname='genomes' ] [ -keep_tmp ] sp1, sp2, ...\n";
	exit(0);
}

$tmpdir = "tmp_mbgd_genetab" if (! $tmpdir);

$outname = "genomes" if (! $outname);

make_path($tmpdir);

$hash_protseq = &get_proteinseq(@species);
&output_all(@species);

foreach $suffix ("genetab", "faa") {
	my($outfile) = "$outname.$suffix";
	open(O, ">$outfile") || die "Can't open $outfile for output\n";
	foreach $sp (@species) {
		my($infile) = "$tmpdir/$sp.$suffix";
		open(F, "$infile") || die "Can't open $infile\n";
		while (<F>) {
			print O $_;
		}
		close(F);
		if (! $keep_tmp) {
			unlink "$tmpdir/$sp.$suffix";
		}
	}
	close(O);
}
remove_tree($tmpdir) if (! $keep_tmp);

sub output_all {
	#output genetab and protseq
	my(@species) = @_;
	foreach $sp (@species) {
		my($genome) = &get_data($sp, 'genome');
		$genome = $genome->[0];

		my($outname) = "$tmpdir/$sp";
		open(GENETAB, ">$outname.genetab") || die "Can't open genetab for output\n";

		print GENETAB "##Genome\tsp:$genome->{sp}";
	       	print GENETAB "\tabbrev:$genome->{abbrev}" if ($genome->{abbrev});;
	       	print GENETAB "\torgname:$genome->{orgname}" if ($genome->{orgname});;
	       	print GENETAB "\tstrain:$genome->{strain}" if ($genome->{strain});;
		print GENETAB "\n";

		my($chromInfo) = &get_data($sp, 'chromosome');
		my($geneInfo) = &get_data($sp, 'gene');
		my(@chroms);

		foreach $c ( @{$chromInfo} ) {
			push(@chroms, $c);
		}
		@chroms = sort {$a->{seqno} <=> $b->{seqno}} @chroms;
		my(@genes);
		my($hash_genetab) = {};

		foreach $g ( @{$geneInfo} ) {
			if ($g->{type} eq 'CDS') {
				# output only when the gene is preset in protseq
				if ($hash_protseq->{$g->{name}}) {
					if (! $g->{aalen}) {
						$g->{aalen} = int( (($g->{to1} - 3) - $g->{from1} + 1) / 3 );
					}
					$g->{pos} = ($g->{from1} + ($g->{to1}-3)) / 2;
					push(@{$genes[$g->{seqno}-1]}, $g);
					$hash_genetab->{$g->{name}} = 1;
				}
			}
		}
		foreach $c (@chroms) {
			print GENETAB "##Chromosome\t";
			print GENETAB "name:$c->{name}" if ($c->{name});
			print GENETAB "\tseq_length:$c->{seq_length}" if ($c->{seq_length});
			print GENETAB "\n";
			foreach $g (@{ $genes[$c->{seqno} - 1] }) {
				print GENETAB join("\t", $g->{sp}, $g->{name}, $g->{aalen}, $g->{pos}, $g->{dir}), "\n"
			}
		}
		open(PROTSEQ, ">$outname.faa") || die "Can't open protseq for output\n";
		$pseq = &get_data($sp, 'proteinseq');
		foreach $g ( @{$pseq} ) {
			# output only when the gene is preset in genetab
			if ($hash_genetab->{$g->{name}}) {
				print PROTSEQ ">$g->{sp}:$g->{name} $g->{descr}\n";
				print PROTSEQ "$g->{seq}\n";
			}
		}
		close(GENETAB);
		close(PROTSEQ);
	}
}
sub get_proteinseq {
	my(@species) = @_;
	my(%hash_protseq);
	foreach $sp (@species) {
		$pseq = &get_data($sp, 'proteinseq');
		foreach $g ( @{$pseq} ) {
			$hash_protseq{$g->{name}} = $g;
		}
	}
	\%hash_protseq;
}
sub get_data {
	my($sp, $type) = @_;
	my $url = "$MBGDAPI/$type/$sp";

	sleep(0.3);

	my $http = HTTP::Tiny->new(timeout => 15);
	my $max_retries = 3;
	my $response;

	for my $attempt (1 .. $max_retries) {
		$response = $http->get($url);

		last if $response->{success};

		warn "Attempt $attempt failed for $url: $response->{status} $response->{reason}. Retrying...\n";
		sleep(1);
	}

	unless ($response->{success}) {
		die "HTTP Error fetching $url: Status $response->{status} - $response->{reason}\n";
	}

	my $data = eval { decode_json($response->{content}) };
	if ($@) {
		die "JSON Parse Error for $url: $@\nResponse content was:\n$response->{content}\n";
	}

	return $data;
}
