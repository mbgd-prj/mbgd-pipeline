#!/bin/sh

dbname=$1
if [ "$dbname" == "" ]; then
	echo Usage: $0 dbname
	exit
fi
protfile=$dbname.fas
qseqfile=$dbname.qseq
qsub_script=qsub_diamond.sh
bindir=`dirname $0`
export PATH="${bindir}:$PATH"

domclust_dir=~/project/domclust/domclust/domclust
mergetree_dir=~/project/domclust/mergeTree/

SPECIES='H.pylori'
GENUS='Helicobacter'

function exec_prokka_all {
	seqdir=$1
	let idx=1
	for genomefile in $seqdir/*.fasta; do
		name=`printf MyData_%05d $idx`
		exec_prokka $genomefile $name
		let idx++
	done
}

function exec_prokka {
	infile=$1
	data_name=$2

	PROKKA_PATH=~/package/Bio/prokka/prokka-master
	PROKKA_CMD=$PROKKA_PATH/bin/prokka
	PROKKA_OPT="-species $SPECIES -genus $GENUS -compliant $ADDOPT "
	PROKKA_OUT="prokka"

	strain=`basename $infile .fasta`
	outdir=$PROKKA_OUT/$data_name
	PROKKA_OPT="$PROKKA_OPT --strain '$strain' --outdir $outdir --pref $data_name --locustag $data_name"
	qsub  <<EOF
	cd \$PBS_O_WORKDIR
	$PROKKA_CMD $PROKKA_OPT $infile
EOF
}
function exec_cdhit {
	infile=$1
	qseqfile=$2
	cdhit_opt="-c 0.98 -aS 0.98"
	qsub  -l mem=12gb,ncpus=12 -W block=true <<EOF
	cd \$PBS_O_WORKDIR
	cd-hit -i $infile -o $qseqfile -M 5000 -T \$NCPUS -d 0 -p 1 $cdhit_opt
EOF
}
function exec_diamond {
	qseqfile=$1
	dbname=$2
	outdir=output_${dbname}
	outfile=${dbname}.homfile

	diamond makedb --in $qseqfile --db $dbname
	$bindir/split_seq.pl -BLOCK_SIZE=1000000 -SEARCH_OPTIONS="--outfmt 6 --evalue 1e-10" -SEARCH_DB=$dbname -QSUB_SCRIPT_FILE=$qsub_script -RESULT_OUT_DIR=$outdir $qseqfile
	qsub -W block=true -V $qsub_script
	cat $outdir/homfile* > $outfile
	
}
function exec_domclust() {
	homfile=$1
	genefile=$2
	clustout=$3
	preclustfile=$4

	$domclust_dir/src/bin/domclust ${homfile} ${genefile} -S80 -C150 -V.8 -ai.9 -ao.6 -HO -n1 -o11 \
		-OgeneClustFile=${preclustfile}> ${clustout}.o11
	$mergetree_dir/bin/mergetree -o0 ${clustout}.o11 > ${clustout}.o1
	$bindir/clusttree2flat.pl ${clustout}.o1 >${clustout}.o0
}
function exec_corealigner() {
	genefile=$1
	clustfile=$2
	corefile=$3
	$bindir/corealign -domclustIn $clustfile $genefile > $corefile
}

##exec_prokka_all seq

#exec_cdhit $protfile $qseqfile
#$bindir/cdhit2clust.pl ${qseqfile}.clstr > ${dbname}.preclust

#exec_diamond $qseqfile $dbname

#$domclust_dir/src/build_input/fasta2genefile.pl $qseqfile

#exec_domclust $dbname.homfile $qseqfile.gene $dbname.cluster $dbname.preclust

exec_corealigner $dbname.cluster.o0 $qseqfile.gene $dbname.coaln

