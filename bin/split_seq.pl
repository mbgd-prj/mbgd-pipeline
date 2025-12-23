#!/usr/bin/perl -s
use File::Basename;

$BLOCK_SIZE=100000 if (! $BLOCK_SIZE);

$QSUB_SCRIPT_FILE = "qsub_blast.sh" if (! $QSUB_SCRIPT_FILE);

$query_file = $ARGV[0];

die "Usage: $0 [-BLOCK_SIZE=($BLOCK_SIZE)] [-QUERY_OUT_DIR=(query_<query_file_name>)] [-RESULT_OUT_DIR=output_<query_file_name>)] [-QSUB_SCRIPT_FILE=($QSUB_SCRIPT_FILE)] query_file\n" if (! $query_file);

$QFILE_NAME =  basename($query_file);
$QFILE_NAME =~ s/\.\w+$//;  #remove suffix
$QUERY_OUT_DIR = "query_$QFILE_NAME" if (! $QUERY_OUT_DIR);
$OUTPUT_FILENAME = "$QUERY_OUT_DIR/$QFILE_NAME";

## to create qsub script
$RESULT_OUT_DIR = "output_${QFILE_NAME}" if (! $RESULT_OUT_DIR);
$SEARCH_OPTIONS = "-outfmt 6 -evalue 0.001" if (! $SEARCH_OPTIONS);
$SEARCH_DB = qq("##ENTER DB NAME##") if (! $SEARCH_DB);
$MAX_EXEC = 32;
$DEFAULT_NCPUS = 4;
$DEFAULT_MEMORY = 12;

$NCPUS = $DEFAULT_NCPUS if (! $NCPUS); 
$MEMORY = $DEFAULT_MEMORY if (! $MEMORY); 


mkdir $QUERY_OUT_DIR;

$fn = 1;
while (<>) {
	if (/^>/){
		if ($len >= $BLOCK_SIZE) {
			&output_query($out, $fn);
			$fn++;
			$out = '';
			$len = 0;
    		}
	} else {
		$seq = $_;
		chomp $seq;
		$seq =~ s/\s+//g;
		$len += length($seq);
	}
	$out .= $_;
}

if ($len > 0) {
	&output_query($out, $fn);
}

&output_qsub_script($fn);

sub output_query {
	my($out, $fn) = @_;
	open(O, ">$OUTPUT_FILENAME.$fn");
	print O $out;
	close(O);
}

## create qsub script
sub output_qsub_script {
	my($num_split) = @_;
	my($NCPU_SETTING);
	if ($NCPUS) {
		$NCPU_SETTING = "#PBS -l ncpus=$NCPUS";
	}
	if ($MEMORY) {
		$MEM_SETTING = "#PBS -l mem=${MEMORY}gb";
	}

	open(O, ">${QSUB_SCRIPT_FILE}");
	print O <<SCRIPT_END;
#!/bin/sh
#PBS -J 1-$num_split
#PBS -N blastjob
#PBS -S /bin/sh
#PBS -V
$NCPU_SETTING
$MEM_SETTING

cd \$PBS_O_WORKDIR
DB=${SEARCH_DB}
SEQDIR=${QUERY_OUT_DIR}
RESULT_OUT_DIR=${RESULT_OUT_DIR}
QUERY_OUT_DIR=${QUERY_OUT_DIR}

INFILE=\$QUERY_OUT_DIR/${QFILE_NAME}.\$PBS_ARRAY_INDEX
OUTFILE=${QFILE_NAME}.\$PBS_ARRAY_INDEX

if [ ! -d \$RESULT_OUT_DIR ]; then
	mkdir \$RESULT_OUT_DIR
fi
##blastp -db \$DB -query \$INFILE -out \$RESULT_OUT_DIR/blastp_\$OUTFILE -num_threads \$NCPUS $SEARCH_OPTIONS
diamond blastp --db \$DB --query \$INFILE --out \$RESULT_OUT_DIR/diamond_\$OUTFILE --threads \$NCPUS $SEARCH_OPTIONS

blast2homfile.pl -skip_sort \$RESULT_OUT_DIR/diamond_\$OUTFILE >\$RESULT_OUT_DIR/homfile_\$OUTFILE

SCRIPT_END
}
