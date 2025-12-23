#!/usr/bin/perl -s

while(<>) {
    if (/^HomCluster/) {
        if ($output) {
            print "$output\n";
            $output = '';
        }
        $output = $_;
    }
    elsif (/^Cluster/) {
        if ($output !~ /^HomCluster\s+\d+[\r\n]*$/) {
            $output .= "\n";
        }
        $output .= $_;
    }
    elsif (/^$/) {
    }
    elsif (/[\*\+]\-\s+(.+)/) {
        $name = $1;
        if (! $Found{$name}) {
            ## eliminating duplicated entries
            $output .= "$name\n";
            $Found{$name} = 1;
        }
    }
}
if ($output) {
    print "$output\n";
}
