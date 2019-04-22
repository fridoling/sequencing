#!/bin/bash

shopt -s extglob

EXPFOLDER=~/Desktop/arbeit/IFOM/DataTest
DATE=20192104
SRC=/lustre/data/ANC/NGS/sequencing/src
BIN=/lustre/data/ANC/NGS/sequencing/bin
REF=/lustre/data/ANC/NGS/sequencing_nongit/data/ref_genomes/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa

FASTAFOLDER=$EXPFOLDER/rawdata/
TEMP=$EXPFOLDER/processed
BAMFOLDER=$EXPFOLDER/bamfiles

ls $FASTAFOLDER*R1_001.fastq.gz | sed 's/_R1_001.*//' |\
(while read FILENAME; do
   F1="$FILENAME"_R1_001.fastq.gz 
   F2="$FILENAME"_R2_001.fastq.gz 

   # get quality info for unprocessed sequences
   $BIN/FastQC/fastqc $F1 $F2 --outdir=$FASTAFOLDER

   # trim fasta sequences
   $SRC/prepare_fasta.sh $F1 $F2 -o $FASTAFOLDER
done)

ls $FASTAFOLDER*.trim_1P.fq.gz | sed 's/.trim_1P.fq.gz//' |\
(while read FILENAME; do
   F1="$FILENAME".trim_1P.fq.gz
   F2="$FILENAME".trim_2P.fq.gz
   $SRC/fasta2bam.sh $F1 $F2 -r $REF -o $BAMFOLDER
done)

ls $BAMFOLDER/*.bam |
(while read FILENAME; do
  $SRC/realign_bam.sh $FILENAME -r $REF
done)
