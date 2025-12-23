#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM -s SPECIES_FILE
";

my %OPT;
getopts('s:g:', \%OPT);

my %SPECIES = ();
if ($OPT{s}) {
    my $species_file = $OPT{s};
    open(SPECIES_FILE, $species_file) || die;
    while (<SPECIES_FILE>) {
	chomp;
	my ($sp) = split;
	$SPECIES{$sp} = 1;
    }
    close(SPECIES_FILE);
}

!@ARGV && -t and die $USAGE;
while (<STDIN>) {
    chomp;
    my @f = split("\t", $_, -1);
    if ($f[0] =~ /^\d+$/) {
	my $sp = $f[2];
	my $name = $f[3];
	my $from = $f[4];
	my $to = $f[5];
	my $motlib = $f[6];
	my $motid = $f[7];
	my $motname = $f[8];
	if ($SPECIES{$sp}) {
	    print join("\t", "$sp:$name", $from, $to, $motlib, $motid, $motname), "\n";
	}
    } else {
	next;
    }
}
