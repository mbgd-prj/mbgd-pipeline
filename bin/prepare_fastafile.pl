#!/usr/bin/perl -s

use File::Basename;

while(<>) {
	if (/^>\s*(\S+)\s*(.*)$/) {
		$name = $1;
		$rest = $2;
		if ($name =~ /:/) {
			($spname, $name) = split(/:/, $name);
		}
		if ($SPNAME) {
			$spname = $SPNAME;
		} else {
			$spname = basename(dirname($ARGV));
		}
		if ($lower_spname) {
			$spname = lc($spname);
		}
		print ">", join(" ", "$spname:$name", $rest), "\n";
	} else {
		print;
	}
}
