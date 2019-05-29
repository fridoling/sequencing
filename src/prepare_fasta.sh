#!/bin/bash

usage()
{
  echo "Usage: ${0##*/} fasta1 fasta2 [-o OUTDIR]"
  exit 2
}

if [ $# -lt 2 ]; then
    usage
fi

FASTA1=$1
FASTA2=$2

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

unset DATE OUTDIR
OPTIND=3
while getopts 'o:?h' c
do
  case $c in
    o) set_variable OUTDIR $OPTARG ;;
    h|?) usage ;;
  esac
done

# check whether input files exist
if [ ! -f $FASTA1 ] || [ ! -f $FASTA2 ]; then
   echo "Error: Please provide paths to existing fasta files." 
   exit 2
fi

# if output directory is not given, use current directory 
if [ -z $OUTDIR ]; then
   OUTDIR='.'
fi

# if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]; then
   mkdir -m 775 $OUTDIR
fi

# path to sequencing repository
SRC="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PROJECT="$(dirname $SRC)"

# path to used programs
TRIMMO=$PROJECT/bin/Trimmomatic-0.36/trimmomatic-0.36.jar

# path to adapter sequences
ADAPTERS=$PROJECT/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa
#ADAPTERS=$PROJECT/bin/Trimmomatic-0.36/adapters/Novogene.fa

# get file and samplenames for fasta files
RAWFASTAPATH=$(dirname "${FASTA1}")
F1NAME=${FASTA1##*/}
F2NAME=${FASTA2##*/}
BASENAME=${F1NAME%%.*}

# trim fasta files
java -jar $TRIMMO PE -threads 2 -phred33 -trimlog $OUTDIR"/$BASENAME".trim_log.txt \
   $FASTA1 $FASTA2 \
   -baseout $OUTDIR"/$BASENAME".trim.fq.gz  \
   ILLUMINACLIP:$ADAPTERS:2:30:8:2:TRUE SLIDINGWINDOW:5:20

