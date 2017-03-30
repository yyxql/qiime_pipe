#!/usr/bin/env bash 
# Chimera detection
#
# Usage: ./05.chimera_remove.sh
#
# Returns an OTU table in JSON BIOM format and a phylogeny with chimeric OTUs removed 
#
# Requires: ChimeraSlayer, Qiime 1.9+, biom, bc
source params.txt

chimsl=${QIIME_OUTDIR}/chimeraslayer
if [ ! -d "$chimsl" ]; then mkdir $chimsl;fi

# Check if core OTUs file is gzip'd
core_gz=`file $CORESET | grep gzip`
if [ ! -z "$core_gz" ]; then gunzip $CORESET; fi

# Identify chimeras using a reference set
ChimeraSlayer --query_NAST $QIIME_OUTDIR/pynast_aligned_seqs/rep_set_aligned.fasta --exec_dir $chimsl &> $chimsl/chimeraslayer.log
grep YES $chimsl/*.CPC | awk '{print $2'} > $chimsl/chimeric_otus.txt

# Filter PYNAST aligned sequences for chimeras
filter_fasta.py -v -f $QIIME_OUTDIR/pynast_aligned_seqs/rep_set_aligned.fasta -o $chimsl/non_chimeric_rep_set_aligned.fasta -s $chimsl/chimeric_otus.txt -n
filter_alignment.py -v -o $chimsl/pynast_aligned_seqs -i $chimsl/non_chimeric_rep_set_aligned.fasta

# Make the otu table command
make_otu_table.py -v -i $QIIME_OUTDIR/final_otu_map_mc2.txt -o $QIIME_OUTDIR/otu_table_mc2_nochimeras.biom -e $chimsl/chimeric_otus.txt

# Add taxa to OTU table command
biom add-metadata -i $QIIME_OUTDIR/otu_table_mc2_nochimeras.biom --observation-metadata-fp=$QIIME_OUTDIR/uclust_assigned_taxonomy/rep_set_tax_assignments.txt --sc-separated=taxonomy --observation-header=OTUID,taxonomy -o $QIIME_OUTDIR/otu_table_mc2_w_tax_nochimeras.biom

# Make phylogeny
make_phylogeny.py -i $chimsl/pynast_aligned_seqs/non_chimeric_rep_set_aligned_pfiltered.fasta -r midpoint -o $QIIME_OUTDIR/rep_set_rooted_nonchimeric.tre

# Convert to JSON biom format for compatability
biom convert -i ${QIIME_QIIME_OUTDIR}/otu_table_mc2_w_tax_nochimeras.biom -o ${QIIME_QIIME_OUTDIR}/otu_table_mc2_w_tax_nochimeras_json.biom --table-type="OTU table" --to-json

# Zip $CORESET
if [ ! -z "$core_gz" ]; then gzip $CORESET; fi

## Chimera stats
# Write biom table summary for all reads biom file and chimeras removed biom file
biom summarize-table -i ${QIIME_OUTDIR}/otu_table_mc2.biom > /tmp/all.sum
sed '1,15d' /tmp/all.sum | sort | cut -f1 -d ':' > /tmp/all.ids
sed '1,15d' /tmp/all.sum | sort | cut -f2 -d ' ' > /tmp/all.reads
sum_reads_all=`grep "Total count" /tmp/all.sum | cut -f3 -d ' '`

biom summarize-table -i ${QIIME_OUTDIR}/otu_table_mc2_nochimeras.biom > /tmp/nochim.sum
sed '1,15d' /tmp/nochim.sum | sort | cut -f1 -d ':' > /tmp/nochim.ids
sed '1,15d' /tmp/nochim.sum | sort | cut -f2 -d ' ' > /tmp/nochim.reads
sum_reads_nochim=`grep "Total count" /tmp/nochim.sum | cut -f3 -d ' '`

# Chimeric reads statistics output to chimera_stats.csv 
diff_id=`diff -q /tmp/all.ids /tmp/nochim.ids`
if [ -z "$diff_id" ]; then 
	# CSV header
	echo -e "id\tinput_reads\toutput_reads\tdifference\tp_chimeric" > ${QIIME_OUTDIR}/chimera_stats.csv
	# Iterate through parsed biom summary files and collate stats
	n_samples=`wc -l /tmp/nochim.reads | cut -f1 -d ' '`
	for i in `seq 1 1 $n_samples`; do
		id=`sed -n "${i}p" /tmp/all.ids`
		all_reads=`sed -n "${i}p" /tmp/all.reads`
		nochim_reads=`sed -n "${i}p" /tmp/nochim.reads`
		reads_diff=`echo "$all_reads - $nochim_reads" | bc`
		p_chim=`echo "1 - ${nochim_reads}/${all_reads}" | bc`
		echo -e "${id}\t${all_reads}\t${nochim_reads}\t${reads_diff}\t${p_chim}" >> ${QIIME_OUTDIR}/chimeraslayer/chimera_stats.csv
	done
	# Chimeric reads summed over all samples
	p_sum_chimeras=`echo "1 - $sum_reads_nochim / $sum_reads_all" | bc`
	diff_sum_chimeras=`echo "$sum_reads_all - $sum_reads_nochim" | bc`
	echo "${diff_sum_chimeras}/${sum_reads_all} chimeric reads (${p_sum_chimeras})"
else
	echo "Unable to calculate the number of chimeras present. ${QIIME_OUTDIR}/otu_table_mc2_nochimeras.biom and ${QIIME_OUTDIR}/otu_table_mc2.biom have different sample IDs."
	exit 0
fi
