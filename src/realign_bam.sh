#!/bin/bash

usage()
{
  echo "Usage: ${0##*/} bam -r REF"
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

BAM=$1

unset REF
OPTIND=2
while getopts 'r:?h' c
do
  case $c in
    r) set_variable REF $OPTARG ;;
    h|?) usage ;;
  esac
done

# check whether input file exists
if [ ! -f $BAM ]; then
   echo "Error: Please provide paths to an existing BAM file." 
   exit 2
fi

# check whether reference is provided
if [ -z $REF ]; then
   echo "Error: Please provide path to reference genome"
   usage
   exit 2
fi 

# path to sequencing repository
SRC="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PROJECT="$(dirname $SRC)"

PICARD=$PROJECT/bin/picard.jar
GATK=$PROJECT/bin/GenomeAnalysisTK.jar

# add readgroups to bam file using picard
# TODO: any meaningful info to put there?
echo "adding readgroups..."
BAM_RG="${BAM%.*}"_RG.bam
java -jar $PICARD AddOrReplaceReadGroups \
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
java -jar $GATK -T RealignerTargetCreator \
-R $REF \
-I $BAM_RG \
-o $IV
echo "...done."

# realign around indels
echo "realigning..."
BAM_OUT="${BAM%.*}"_realigned.bam
java -jar $GATK -T IndelRealigner \
-R $REF \
-targetIntervals $IV \
-I $BAM_RG \
-o $BAM_OUT
echo "...done."
