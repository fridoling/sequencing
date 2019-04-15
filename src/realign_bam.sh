#!/bin/bash

usage()
{
  echo "Usage: $0 [ -f REF ] [ -b BAM]"
  exit 2
}

set_variable()
{
  local varname=$1
  shift
  if [ -z "${!varname}" ]; then
    eval "$varname=\"$@\""
  else
    echo "Error: $varname already set"
    usage
  fi
}

unset REF BAM

while getopts 'f:b:?h' c
do
  case $c in
    f) set_variable REF $OPTARG ;;
    b) set_variable BAM $OPTARG ;;
    h|?) usage ;;
  esac
done

# add readgroups to bam file using picard
# TODO: any meaningful info to put there?
echo "adding readgroups..."
BAM_RG="${BAM%.*}"_RG.bam
java -jar /lustre/home/fgross/bin/GenomeAnalysisTK-3.8-0-ge9d806836/picard.jar AddOrReplaceReadGroups \
I=$BAM \
O=$BAM_RG RGLB=lib RGPL=illumina RGPU=NA RGSM=30
echo "...done."

# index resulting bam file
echo "indexing resulting file..."
samtools index $BAM_RG
echo "...done."

# create target intervals:
echo "creating  target intervals..."
IV="${BAM%.*}".intervals
java -jar /lustre/home/fgross/bin/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T RealignerTargetCreator \
-R $REF \
-I $BAM_RG \
-o $IV
echo "...done."

# realign around indels
echo "realigning..."
BAM_OUT="${BAM%.*}"_realigned.bam
java -jar /lustre/home/fgross/bin/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar -T IndelRealigner \
-R $REF \
-targetIntervals $IV \
-I $BAM_RG \
-o $BAM_OUT
echo "...done."
