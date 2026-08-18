#!/usr/bin/perl 

my $CUTOFF = shift;
my $simf = shift;

if(!-e "$CUTOFF"){
   print STDERR "$0 <mmseqs2> <cutoff>\n";
   exit;
}
$simf = "-" if (! $simf);
  
my $coff=read_cutoff($CUTOFF);

open(FH, "$simf") || die("Can not open $simf");
while(<FH>) {
    next if /^#/;
    $_=~s/\r\n|\n|\r//ig;
    $_=~s/ +/\t/ig;
    my($qname, $tname, $ident, $alen, $nmis, $ngap, $from1, $to1, $from2, $to2,
           $eval, $score, $qlen, $slen, $desc) = split(/\t/,$_,15);

    if($coff->{"$tname"} && $score > $coff->{"$tname"}){
       print "$_\n";
    }
}
close(FH);

exit;


sub read_cutoff{
    my ($infile)=@_;
    my $info={};
    open my $fh, '<', $infile or die ;
    while(my $n=<$fh>){
        $n=~s/\r\n|\r|\n//;
        next if $n=~/^#/;
        my ($clustid,$thresh,$etc)=split(/\t/,$n,3);
        if($clustid){
           $info->{"$clustid"}=$thresh;
        }
    }
    close($fh);
    return $info;
}
