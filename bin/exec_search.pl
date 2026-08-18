#!/bin/bash
QSEQ=GCA_000802015.1.faa
PROFDB=prof/mbgd
MMOUT=`basename $QSEQ .faa`.mmout
OUTPUT=`basename $QSEQ .faa`.out
TITFILE=mbgdcluster.tit
CUTOFF=cutoff_min3
TMP=tmp
mmseqs easy-search $QSEQ $PROFDB $MMOUT $TMP
./exec_cutoff.pl $CUTOFF $MMOUT | ./gettop2.pl | ./addtit.pl $TITFILE > $OUTPUT
