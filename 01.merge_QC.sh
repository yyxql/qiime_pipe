#!/usr/bin/env bash 
# Merge paired end reads
#
# Usage: ./01.merge_QC.sh
#
# Outputs merged paired end reads with primers trimmed, discarding sequences with low-quality
# Illumina quality scores. 
#
# Requires: usearch, cutadapt
#
source params.txt

if [ ! -d $SEQ_PATH/ ]; then mkdir $SEQ_PATH;fi

for i in `find $RAW_SEQPATH/ -maxdepth 1 -iname "*${READF_ID}*fastq.gz"`; do
 R1=`echo $i | sed 's;.gz;;'` 
 R2=`echo $R1 | sed "s/${READF_ID}/${READR_ID}/"`

 if [ ! -f $R2.gz ]; then echo "$R2 is missing. Skipping to next read pair."; continue; fi
 sampleid=`basename $i | cut -d _ -f 1` # sampleid is the first field in the filename w/delimiter set to _

 gunzip ${R1}.gz ${R2}.gz
 
# Merge sequences and quality filter
 usearch -fastq_mergepairs $R1 -reverse $R2\
 -fastqout /tmp/usearch_merge.fastq\
 -fastq_allowmergestagger -fastq_truncqual $TRUNC -fastq_merge_maxee $MAXEE -fastq_minlen MINLENGTH\
 2> $SEQ_PATH/${sampleid}.usearch
 usearch -fastq_eestats /tmp/usearch_merge.fastq -output $SEQ_PATH/${sampleid}_merged.fastq.stats

# Trim primer
 cutadapt -g $PRIMER_F -a $PRIMER_R\
 -o ${SEQ_PATH}/${sampleid}_merged.fastq.gz\
 -e ${ERROR_RATE}\
 --trimmed-only\
 /tmp/usearch_merge.fastq\
 > ${SEQ_PATH}/${sampleid}.cutadapt
 rm /tmp/usearch_merge.fastq
 gzip $R1 $R2 
done

# Generate statistics
echo "SampleID,total_pairs,pairs_merged,p_combined,p_wprimer,final_reads"\
       	> merge_QC_stats.csv

for i in `find $SEQ_PATH/ -maxdepth 1 -iname "*.usearch"`; do
 sampleid=`basename $i | cut -d _ -f 1 | sed 's/.usearch//'` 
 total_pairs=`grep "Pairs" $i | sed -e 's; ;;g' -e 's;[A-Z];;g' -e 's;[a-z];;g'`
 pairs_merged=`grep "Converted" $i |sed 's/^[ \t]*//' |cut -d " " -f 1`
 pcombined=`grep "Converted" $i |cut -d "(" -f 2 |sed 's;%);;'`
 p_wprimer=`grep "Reads written" ${SEQ_PATH}/${sampleid}.cutadapt |cut -d '(' -f 3 |sed 's;%);;'`
 final_reads=`grep "Reads written" ${SEQ_PATH}/${sampleid}.cutadapt | sed 's/ //g' |\
 cut -d ':' -f 2 | cut -d '(' -f 1 |sed 's/,//'`
 echo "Proportion merged: $pcombined" >> $i
 echo "$sampleid,$total_pairs,$pairs_merged,$pcombined,$p_wprimer,$final_reads"\
	>> merge_QC_stats.csv
 
# Calculate sequence lengths
 sampleid=`basename $i | cut -d _ -f 1` 
 grep "Merged length" $i >> ${SEQ_PATH}/${sampleid}.seqlengths
done
