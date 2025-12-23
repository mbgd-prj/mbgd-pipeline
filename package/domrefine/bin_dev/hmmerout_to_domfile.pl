#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM
renumber domain IDs (increasing order of position)
";

use DomRefine::Read;
use DomRefine::Refine;

### Settings ###
my %OPT;
getopts('', \%OPT);

my $TMP_DCLST = define_tmp_file("$PROGRAM.dclst");
END {
    remove_tmp_file($TMP_DCLST);
}

### Main ###
-t and die $USAGE;
open(TMP_DCLST, ">$TMP_DCLST") || die;
while (<STDIN>) {
    if (/^#/) {
	next;
    }
    my @f = split;
    print TMP_DCLST "$f[4] $f[0] 0 $f[17] $f[18]\n";
}
close(TMP_DCLST);

# system "cat $TMP_DCLST";
print renumber_domain($TMP_DCLST);
