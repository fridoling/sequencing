#!/bin/bash

EXPFOLDER=/lustre/data/ANC/NGS/20190516_GA12/
DATE=20190516
SRC=/lustre/data/ANC/NGS/sequencing/src
BIN=/lustre/data/ANC/NGS/sequencing/bin
REF=/lustre/data/ANC/NGS/ref/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa

FASTAFOLDER=$EXPFOLDER/rawdata/
TEMP=$EXPFOLDER/processed
BAMFOLDER=$EXPFOLDER/bam_pipeline
LOGFOLDER=$EXPFOLDER/log

# unique name of logfile that includes date and time
LOGFILE=$LOGFOLDER/log`date +%Y%m%d_%H%M%S`.txt

# if log directory does not exist, create it
if [ ! -d $LOGFOLDER ]; then
   mkdir -m 775 $LOGFOLDER
fi

{
    echo -e "START EVERYTHING AT " $(date) "\n\n"

#    # rename fasta files using the date of experiment
#    echo -e "*** \t \t Step 1: rename original fasta files  \t \t ***"
#    $SRC/rename_fasta.sh $FASTAFOLDER $DATE
#    echo -e "*** RENAMING: done \n\n\n"

    ls $FASTAFOLDER*_1.fq.gz | sed 's/_1.fq.gz//' |\
	(while read FILENAME; do
	     F1="$FILENAME"_1.fq.gz 
	     F2="$FILENAME"_2.fq.gz

	     # check gz file integrity
	     if [ gzip -t $F1 ] && [ gzip -t $F2 ] ; then
		 echo -e "\n\n*** trim $FILENAME"
		 # get quality info for unprocessed sequences
		 $BIN/FastQC/fastqc $F1 $F2 --outdir=$LOGFOLDER
	     
		 # trim fasta sequences
		 $SRC/prepare_fasta.sh $F1 $F2 -o $FASTAFOLDER
	     else
		 echo -e "\n\nWARNING: $F1 or $F2 corrupted "
	     fi
	 done)
    echo -e "*** TRIMMING: done \n\n\n"

    
    ls $FASTAFOLDER*.trim_1P.fq.gz | sed 's/.trim_1P.fq.gz//' |\
	(while read FILENAME; do
	     F1="$FILENAME".trim_1P.fq.gz
	     F2="$FILENAME".trim_2P.fq.gz

	     echo -e "*** prepare bam from $FILENAME"

	     $SRC/fasta2bam.sh $F1 $F2 -r $REF -o $BAMFOLDER
	 done)
    echo -e "*** BAM CREATION: done \n\n\n"


    
    ls $BAMFOLDER/*.bam |
	(while read FILENAME; do
	     echo -e "*** realign $FILENAME"
	     $SRC/realign_bam.sh $FILENAME -r $REF

	     REALIGNEDFILE="${FILENAME%%.*}"_realigned.bam
	     echo -e "*** get information about $REALIGNEDFILE quality"
	     $BIN/qualimap_v2.2.1/qualimap bamqc -bam $REALIGNEDFILE -nw 400 -hm 3

	 done)
    echo -e "*** BAM REALIGNMENT and QUALITY CHECK: done \n\n\n"
    echo -e "EVERYTHING DONE AT " $(date) "\n\n"

} 2>&1 | tee $LOGFILE
