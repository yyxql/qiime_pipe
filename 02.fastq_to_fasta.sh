#!/usr/bin/env bash 
# Convert FASTQ to Qiime-formatted FASTA files and subsample (optional)
#
# Usage: ./03.fastq_to_fasta.sh
#
# Outputs a single Qiime-compliant FASTA file (seqs_concat.fna.gz) using split_libraries_fastq.py from 
# Qiime. If SUBSAMPLE is specified in params.txt, the script randomly subsamples the requested 
# number of sequences from each sample. 
#
# Requires: Qiime 1.9+
source params.txt

function subsample_fasta {
# Subsamples sequences randomly from a FASTA file
# usage: subsample_fasta <FASTA> <sequences> 
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $1 |\
awk 'NR>1{ printf("%s",$0); n++; if(n%2==0) { printf("\n");} else { printf("\t");} }' |\
awk -v k=$2 'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x-1:int(rand()*x);if(s<k)R[s]=$0}END\
	{for(i in R)print R[i]}' |\
sed 's;\t;\n;' 
}

outdir="${SEQ_PATH}/slout"
if [ ! -d $outdir ]; then mkdir -p $outdir; fi

# Convert FASTQ files into QIIME formatted FASTA files
for i in `find $SEQ_PATH -maxdepth 1 -iname "*fastq.gz"`; do 
 sampleid=`basename $i | cut -d _ -f 1`
 split_libraries_fastq.py -i $i -o ${outdir}/${sampleid}\
	 -q 0\
	 --phred_offset 33\
	 --sample_ids $sampleid\
	 --barcode_type 'not-barcoded' 
 gzip -f ${outdir}/${sampleid}/seqs.fna
done

if [ -f $outdir/seqs_concat* ]; then rm $outdir/seqs_concat*; fi

# Subsample FASTA files
if [ "$SUBSAMPLE" == -9 ]; then 
 for i in `find $outdir -iname "seqs.fna.gz"`; do 
	gunzip -c $i >> $outdir/seqs_concat.fna
 done
else
 for i in `find $outdir -iname "seqs.fna.gz"`; do 
	gunzip -c $i | subsample_fasta /dev/stdin $SUBSAMPLE >> $outdir/seqs_concat.fna
 done
fi

# Output statistics
ids=`grep '>' ${outdir}/seqs_concat.fna  | cut -f1 -d '_' | uniq | wc -l`
seqs=`grep -v '>' ${outdir}/seqs_concat.fna | wc -l`
seqs_per_id=`echo $[seqs/ids]`

echo -e "$seqs sequences in $ids samples (mean: $seqs_per_id sequences per sample) in ${outdir}/seqs_concat.fna"

gzip $outdir/seqs_concat.fna
