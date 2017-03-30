#!/usr/bin/env bash 
# Pick OTUS
#
# Usage: ./04.pick_open_otus.sh
#
# Uses open reference OTU picking as implemented in Qiime. 
# Parameters to pass to Qiime are stored in $QIIME_PARAMS file.
#
# Requires: Qiime 1.9+
source params.txt

if [ ! -d "$QIIME_OUTDIR" ]; then mkdir -p $QIIME_OUTDIR; fi
if [ -f "${SEQ_PATH}/slout/seqs_concat.fna.gz" ]; then gunzip ${SEQ_PATH}/slout/seqs_concat.fna.gz; fi

pick_open_reference_otus.py -v -f -a -O $CPU -i ${SEQ_PATH}/slout/seqs_concat.fna -o $QIIME_OUTDIR -p $QIIME_PARAM --prefilter_percent_id=0

gzip ${SEQ_PATH}/slout/seqs_concat.fna 
