#!/bin/bash

# programs
TRIMMOPATH=/lustre/data/ANC/bin/Trimmomatic-0.36/trimmomatic-0.36.jar
ADAPTORPATH=/lustre/data/ANC/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa
# inputs
REF=/lustre/data/ANC/NGS/ref/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
MAINFOLDER=/lustre/data/ANC/NGS/20190225_tub2-150_facility
FASTAFOLDER=$MAINFOLDER'/rawdata/'
# outputs
BAMFOLDER=$MAINFOLDER'/bamfiles/'
BAMPREFIX=20190225__

# # # # # # # # # # # # # # # # # # # # # 
if [ ! -d $BAMFOLDER ]; then
  mkdir -p $BAMFOLDER;
fi

for samplefullfilename in $FASTAFOLDER/*_R1_*.fastq.gz
do
    # take only file name (remove path)
    samplefilename=${samplefullfilename##*/}
    # take sample name, to be used afterwards
    samplename=$( echo $samplefilename | cut -d"_" -f 1 )
    # create file sample name
    sample=$samplename
    sample=$sample"_"$( echo $samplefilename | cut -d"_" -f 2 )
    sample=$sample"_"$( echo $samplefilename | cut -d"_" -f 3 )

    echo
    echo
    echo
    echo PROCESSING $samplename FROM $sample

    # check alignment
    OUTPUTALIGN=$BAMFOLDER$BAMPREFIX$samplename.trim.align.bwa.bam
    fileprova=( $OUTPUTALIGN )
    if [ -f $fileprova ]
    then
	echo "$samplename" has been aligned
    else
	# trimming #
	INPUTTRIM=$FASTAFOLDER"$sample"_R1_001.fastq.gz\ $FASTAFOLDER"$sample"_R2_001.fastq.gz
	OUTPUTTRIM=$BAMFOLDER$samplename.trim_1P.fq.gz\ $BAMFOLDER$samplename.trim_1U.fq.gz\ $BAMFOLDER$samplename.trim_2P.fq.gz\ $BAMFOLDER$samplename.trim_2U.fq.gz
	fileprova=( $OUTPUTTRIM )
	if [ -f $fileprova ]
	then
	    echo "$samplename" has been trimmed
	else
	    # trim if you didn’t already
	    #java -jar $TRIMMOPATH PE -threads 2 -phred33 $INPUTTRIM $OUTPUTTRIM ILLUMINACLIP:$ADAPTORPATH:2:30:8:2:TRUE SLIDINGWINDOW:5:20 2> $BAMFOLDER"$samplename".trim_log.txt
	    java -jar $TRIMMOPATH PE -threads 2 -phred33 -trimlog $BAMFOLDER"$samplename".trim_log.txt $INPUTTRIM $OUTPUTTRIM ILLUMINACLIP:$ADAPTORPATH:2:30:8:2:TRUE SLIDINGWINDOW:5:20
	fi

	# do alignment
    	INPUTALIGN=$BAMFOLDER$samplename.trim_1P.fq.gz\ $BAMFOLDER$samplename.trim_2P.fq.gz
	bwa mem $REF $INPUTALIGN | samblaster -r | samtools view -buS - | samtools sort -m 5000000000 -o $OUTPUTALIGN
	samtools index $OUTPUTALIGN
	# remove trim files
	echo rm $OUTPUTTTRIM
	rm $OUTPUTTTRIM
    fi
done
