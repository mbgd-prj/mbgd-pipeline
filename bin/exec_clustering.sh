#!/bin/sh

dbname=$1
if [ "$dbname" == "" ]; then
	dbname=query
fi
protfile=$dbname.fas
qseqfile=$dbname.qseq
genefile=$dbname.genetab

inseq_dir=seq
annot_dir=in_data
species_file=$inseq_dir/species.txt

# execute prokka to annotate genomes
exec_prokka=0
# reducing data size using cdhit
exec_cdhit=0
# refining the ortholog grouping using domrefine
exec_domrefine=0
# identifying genomic island in addition to core
exec_findIsland=0

lower_spname=0
protid_tag='locus_tag'
#protid_tag='protein_id'

# use PBS for job queuing
with_pbs=0

if [ $exec_cdhit == 1 ]; then
	preclust_file=$dbname.preclust
else
	qseqfile=$protfile
fi

#function realpath {
#        p=$1
#        [[ $p = /* ]] || p=`pwd`/$p
#        echo $p
#}

qsub_script=qsub_diamond.sh
bindir=`dirname $0`
bindir=`realpath $bindir`
export PATH="${bindir}:$PATH"

domclust_dir=${bindir}/../package/domclust/
domrefine_dir=${bindir}/../package/domrefine/

ncpus=8

domrefine_out=domrefine
aliprog=famsa

SPECIES='H.pylori'
GENUS='Helicobacter'

function get_species {
	local id=$1
	SPNAME=`grep $id $species_file | cut -f2`
	GENUS=`grep $id $species_file | cut -f3 | cut -d' ' -f1`
	SPECIES=`grep $id $species_file | cut -f3 | cut -d' ' -f2`
	STRAIN=`grep $id $species_file | cut -f4`
}

function exec_command_by_qsub {
	local comm=$1
	local qsub_opt=$2
	echo $comm
	echo $qsub_opt
	qsub $qsub_opt <<EOF
	cd \$PBS_O_WORKDIR
	$comm
EOF
}

function exec_prokka_all {
	local seqdir=$1
	local outdir=$2
	let idx=1
	if [ ! -d $seqdir ]; then
		echo "Directory not found: $seqdir"
	fi
	for genomefile in $seqdir/*; do
		name=`printf UG%05d $idx`
		exec_prokka $genomefile $name $outdir 
		let idx++
	done
	if [ "$with_pbs" == 1 ]; then
		exit
	fi
}

function exec_prokka {
	local infile=$1
	local data_name=$2
	local PROKKA_OUT=$3

	local outdir=$PROKKA_OUT/$data_name

	filebase=`basename $infile .fasta`
	if [ $species_file ]; then
		get_species $filebase
	fi

	echo "## Execute Prokka"
	PROKKA_PATH=/usr/local
#	PROKKA_CMD=$PROKKA_PATH/bin/prokka
	PROKKA_CMD=prokka
#	PROKKA_OPT="--species $SPECIES --genus $GENUS --compliant $ADDOPT "
	PROKKA_OPT="--species '$SPECIES' --genus '$GENUS' $ADDOPT "
	if [ "$STRAIN" != "" ]; then
		PROKKA_OPT="--strain '$STRAIN'"
	fi
	if [ "$PROKKA_OUT" == "" ]; then
		PROKKA_OUT="prokka"
	fi

	spname=$SPNAME
	if [ -f "$outdir/$spname.faa" ]; then
		return
	fi
	PROKKA_OPT="$PROKKA_OPT --strain '$strain' --outdir $outdir --pref $spname --locustag '$data_name'"
#	echo $PROKKA_CMD $PROKKA_OPT $infile
# 	for PBS
	if [ "$with_pbs" == 1 ]; then
		qsub <<EOF
		cd \$PBS_O_WORKDIR
		$PROKKA_CMD $PROKKA_OPT $infile
EOF
	else
		eval $PROKKA_CMD $PROKKA_OPT $infile
	fi
}
function prepare_from_gbk {
	local annot_dir=$1
	rm -f $protfile $genefile
	for genome_dir in $annot_dir/*; do
		faa=`echo $genome_dir/*faa|sed 's/  */,/g'`
#		echo "cat $genome_dir/*faa | sed 's/^\([^: \t]*\):/\L\1:/' >> $protfile"

		prepare_fastafile.pl -lower_spname=$lower_spname $genome_dir/*faa >> $protfile

#		if [ "$lower_spname" == "1" ]; then
#			cat $genome_dir/*faa | sed 's/^\([^: \t]*\):/\L\1:/' >> $protfile
#		else 
#			cat $genome_dir/*faa >> $protfile
#		fi 

#		gff2genetab.pl $genome_dir/*gff >> $genefile
		echo "gff2genetab.pl -seqfile=$faa $genome_dir/*gff >> $genefile"
		gff2genetab.pl -tag=$protid_tag -lower_spname=$lower_spname -seqfile=$faa $genome_dir/*gff >> $genefile
	done
}

function exec_cdhit {
	local infile=$1
	local qseqfile=$2
	local cdhit_opt="-c 0.98 -aS 0.98"

	local qsub_opt="-l mem=12gb,ncpus=$ncpus -W block=true"
	local command="cd-hit -i $infile -o $qseqfile -M 5000 -T $ncpus -d 0 -p 1 $cdhit_opt"

	echo "## Execute CD-Hit"

	if [ "$with_pbs" == 1 ]; then
		# for PBS
		local command="cd-hit -i $infile -o $qseqfile -M 5000 -T \$NCPUS -d 0 -p 1 $cdhit_opt"
		exec_command_by_qsub "$command" "$qsub_opt"
	else
		echo "$command"
		eval $command
	fi
}
function exec_diamond {
	local qseqfile=$1
	local dbname=$2
	local outdir=output_${dbname}
	local outfile=${dbname}.homfile
	local search_opt="--outfmt 6 --evalue 1e-10 --threads $ncpus --max-target-seqs 5000"

	echo "## Execute Diamond"

	if [ -s $outfile ]; then
		return
	fi

	echo "diamond makedb --in $qseqfile --db $dbname"
	diamond makedb --in $qseqfile --db $dbname

	if [ "$with_pbs" == 1 ]; then
		# for PBS with splitting sequences into chunks
		$bindir/split_seq.pl -BLOCK_SIZE=200000 -SEARCH_OPTIONS="$search_opt" -SEARCH_DB=$dbname -QSUB_SCRIPT_FILE=$qsub_script -RESULT_OUT_DIR=$outdir $qseqfile
		qsub -W block=true -V $qsub_script
		cat $outdir/homfile* > $outfile
	else 
		diamond blastp $search_opt --db $dbname --query $qseqfile >$outfile.tmp
		blast2homfile.pl -skip_sort $outfile.tmp >$outfile
		echo "diamond blastp $search_opt --db $dbname --query $qseqfile >$outfile.tmp"
		echo "blast2homfile.pl -skip_sort $outfile.tmp >$outfile"
		rm $outfile.tmp
	fi
}
function exec_domclust() {
	local homfile=$1
	local genefile=$2
	local clustout=$3
	local preclustfile=$4
	local outfmt=$5

	if [ "$outfmt" == "" ]; then
		outfmt="o0"
	fi

	if [ -s "${clustout}.$outfmt" ]; then
		return
	fi

	echo "## Execute DomClust"

	if [ -f "$preclustfile" ]; then
		preclustOpt="-OgeneClustFile=${preclustfile}"
	fi

	local qsub_opt="-l mem=96gb -W block=true"
	local command=`cat<<EOF
	$domclust_dir/bin/domclust ${homfile} ${genefile} -S80 -C150 -V.8 -ai.9 -ao.6 -HO -n1 -$outfmt $preclustOpt > ${clustout}.$outfmt
EOF`

	if [ "$with_pbs" == 1 ]; then
		# for PBS
		exec_command_by_qsub "$command" "$qsub_opt"
	else
		echo "$command"
		eval $command
	fi
}

function exec_domrefine() {
	local clustout=$1
	local seqfile=$2

	if [ -s ${domrefine_out}/cluster.domrefine.o0 ]; then
		return
	fi

	echo "## Execute DomRefine"

	mkdir -p $domrefine_out
	pushd $domrefine_out
	cwd=`pwd`
	echo "$domrefine_dir/domrefine -p $aliprog -c $cwd/cache -F -a -t -O -N $cwd/../$clustout.o11 $cwd/../$seqfile"
	$domrefine_dir/domrefine -p $aliprog -c $cwd/cache -F -a -t -O -N $cwd/../$clustout.o11 $cwd/../$seqfile
	popd
	if [ ! -f dbname.domrefine.o0 ]; then
		ln -s $domrefine_out/cluster.domrefine.o0 $dbname.domrefine.o0
	fi
}
function exec_corealigner() {
	local clustfile=$1
	local genefile=$2
	local corefile=$3
	local islfile=$4
#	COREALIGN_OPT="-ConsRatio=0.9 -NbrConsRatio=0.6 -NbrConsRatio2=0.6"

	echo "## Execute CoreAligner"

	echo "## Execute CoreAligner"

	echo "$bindir/corealign $COREALIGN_OPT -domclustIn $clustfile $genefile > $corefile"
	$bindir/corealign $COREALIGN_OPT -domclustIn $clustfile $genefile > $corefile
	if [ $exec_findIsland == 1 -a "$islfile" != "" ]; then
		echo "$bindir/findIsland $cluster_out $genefile $corefile > $islfile"
		$bindir/findIsland -domclustIn $corefile $cluster_out $genefile > $islfile
	fi
}

# to create annotated data by Prokka
if [ "$exec_prokka" == 1 ]; then
	exec_prokka_all $inseq_dir $annot_dir
else
	mkdir -p $annot_dir
fi
prepare_from_gbk $annot_dir


if [ "$exec_cdhit" == 1 ]; then
	if [ ! -s ${dbname}.preclust ]; then
		exec_cdhit $protfile $qseqfile
		$bindir/cdhit2clust.pl ${qseqfile}.clstr > ${dbname}.preclust
		$bindir/fasta2genefile.pl $qseqfile 
	fi
	genefile4domclust=$qseqfile.gene
else
	genefile4domclust=$genefile
fi

exec_diamond $qseqfile $dbname

if [ "$exec_domrefine" == 1 ]; then
	domclust_outfmt=o11
else
	domclust_outfmt=o0
fi

exec_domclust $dbname.homfile $genefile4domclust $dbname.cluster "$preclust_file" "$domclust_outfmt"

if [ "$exec_domrefine" == 1 ]; then
	exec_domrefine $dbname.cluster $qseqfile
	cluster_out=$dbname.domrefine.o0
else
	cluster_out=$dbname.cluster.o0
fi

if [ "$exec_cdhit" == 1 ]; then
	restore_seq.pl $cluster_out $preclust_file > ${cluster_out}.restore
	cluster_out=${cluster_out}.restore
fi

exec_corealigner $cluster_out $genefile $dbname.coaln $dbname.isl

