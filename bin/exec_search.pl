#!/bin/bash

QSEQ=$1
PROFDB=$2
BINDIR=`dirname $0`
PATH=${BINDIR}:$PATH

if [ "$QSEQ" == "" ]; then
	echo "Usage: $0 qseq prof"
	exit
fi
if [ "$PROFDB" == "" ]; then
	PROFDB=db/mbgd_prof/mbgd
fi

MMOUT=`basename $QSEQ .faa`.mmout
OUTPUT=`basename $QSEQ .faa`.out
TITFILE=mbgdcluster.tit
CUTOFF=cutoff_min3
TMP=tmp
mmseqs easy-search $QSEQ $PROFDB $MMOUT $TMP
#cutoff_filter.pl $CUTOFF $MMOUT | gettop2.pl | addtit.pl $TITFILE > $OUTPUT
