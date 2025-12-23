#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM FILE
";

my %OPT;
getopts('', \%OPT);

!@ARGV && -t and die $USAGE;
my @input = <>;
unless (@input) {
    @input = ("");
}

shuffle(\@input);
print @input;

sub shuffle {
    my ($a) = @_;
    my $i;
    for ($i=@$a; --$i; ) {
	my $j = int rand ($i+1);
	next if $i == $j;
	@$a[$i,$j] = @$a[$j,$i];
    }
}
