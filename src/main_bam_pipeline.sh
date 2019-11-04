#!/bin/bash

EXPFOLDER=/lustre/data/ANC/NGS/20190516_GA12/
DATE=20190516

EXPFOLDER=/lustre/data/ANC/NGS/20190929_GA5/
DATE=20190929

SRC=/lustre/data/ANC/NGS/sequencing/src
BIN=/lustre/data/ANC/NGS/sequencing/bin
REF=/lustre/data/ANC/NGS/ref/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa

SAMPLELIST=(A1 B1 C1 D1 E1 F1 G1 H1 A2 B2 C2 D2 E2 F2 G2 H2 A4 B4 C4 D4 E4 F4 G4  H4   A5   B5 C5  D5 E5 F5 G5 H5  A6   B6   C6   D6   E6  F6 G6 H6 A7   B7   C7   D7   E7   F7   G7   H7   A8   B8   C8   D8   E8   F8   G8   H8   A9   B9   C9   D9   E9   F9   G9   H9  A10  B10  C10  D10 E10  F10  G10  H10  A11  B11  C11  D11  E11  F11  G11  H11  A12  B12  C12  D12  E12  F12  G12  H12)

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
	(while read SAMPLENAME; do
	     # CHECK if SAMPLE is one of the SAMPLELIST
	     SAMPLE=$(basename "$SAMPLENAME")
	     SAMPLE="${SAMPLE/20190516_/}"
	     isthere=0
	     for i in "${SAMPLELIST[@]}"
	     do
		 if [[ "$SAMPLE" =~ "$i" ]] ; then
		     isthere=1
		 fi
	     done
	     # if the sample is one of the list, run
	     if [ $isthere == 1 ] ; then 
		 
		 echo -e "\n\n\n*** ANALYZE  $SAMPLENAME"
	     
		 # check gz file integrity
		 F1="$SAMPLENAME"_1.fq.gz 
		 F2="$SAMPLENAME"_2.fq.gz
		 if  gzip -t $F1  &&  gzip -t $F2  ; then
		     echo -e "*** trim $SAMPLENAME"
		     # get quality info for unprocessed sequences
		     $BIN/FastQC/fastqc $F1 $F2 --outdir=$LOGFOLDER
		     
		     # trim fasta sequences
		     $SRC/prepare_fasta.sh $F1 $F2 -o $FASTAFOLDER
		 else
		     echo -e "\n\nWARNING: $F1 or $F2 corrupted "
		     continue
		 fi
		 echo -e "*** done \n"

		 ## PREPARE BAM
		 F1="$SAMPLENAME".trim_1P.fq.gz
		 F2="$SAMPLENAME".trim_2P.fq.gz
		 echo -e "*** prepare bam from $SAMPLENAME"
		 $SRC/fasta2bam.sh $F1 $F2 -r $REF -o $BAMFOLDER
		 echo -e "*** done \n"

		 ## REALIGN BAM
		 SAMPLE=$(basename "$SAMPLENAME")
		 BAMFILE=$BAMFOLDER/"$SAMPLE"
		 echo -e "*** realign $BAMFILE".bam
		 $SRC/realign_bam.sh "$BAMFILE".bam -r $REF
		 echo -e "*** done \n"

		 REALIGNEDFILE="$BAMFILE"_realigned.bam
		 echo -e "*** get information about $REALIGNEDFILE quality"
		 $BIN/qualimap_v2.2.1/qualimap bamqc -bam $REALIGNEDFILE -nw 400 -hm 3
		 echo -e "*** done \n\n\n"
	     fi
	 done
	 )
    echo -e "EVERYTHING DONE AT " $(date) "\n\n"
} 2>&1 | tee $LOGFILE
