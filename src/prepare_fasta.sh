#!/bin/bash

usage()
{
  echo "Usage: $0 fasta1 fasta2 [-d DATE] [-o OUTDIR]"
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
while getopts 'd:o:?h' c
do
  case $c in
    d) set_variable DATE $OPTARG ;; 	
    o) set_variable OUTDIR $OPTARG ;;
    h|?) usage ;;
  esac
done

if [ -z $OUTDIR ]; then
   OUTDIR='.'
else
   mkdir -m 775 $OUTDIR
fi


# path to sequencing repository
SRC="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PROJECT="$(dirname $SRC)"

# path to used programs
FASTQC=$PROJECT/bin/FastQC/fastqc
TRIMMO=$PROJECT/bin/Trimmomatic-0.36/trimmomatic-0.36.jar

# path to adapter sequences
ADAPTERS=$PROJECT/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa

# get file and samplenames for fasta files
RAWFASTAPATH=$(dirname "${FASTA1}")
F1NAME=${FASTA1##*/}
F2NAME=${FASTA2##*/}
SAMPLE1=${F1NAME%%_*}
SAMPLE2=${F2NAME%%_*}

if [ -z "$DATE" ]; then
  S1NAME=${F1NAME%%.*}
  S2NAME=${F2NAME%%.*}
else
  S1NAME="$SAMPLE1"_"$DATE"
  S2NAME="$SAMPLE2"_"$DATE"
fi

# get quality info for unprocessed sequences
echo $RAWFASTAPATH
exit 1
$FASTQC $FASTA1 $FASTA2 --outdir=$RAWFASTAPATH

# trim fasta files
java -jar $TRIMMO PE -threads 2 -phred33 -trimlog $OUTDIR"/$S1NAME".trim_log.txt \
   $FASTA1 $FASTA2 \
   -baseout $OUTDIR"/$S1NAME".trim.fq.gz  \
   ILLUMINACLIP:$ADAPTERS:2:30:8:2:TRUE SLIDINGWINDOW:5:20

# get quality info for processed sequences
$FASTQC $OUTDIR"/$S1NAME".trim_1P.fq.gz $OUTDIR"/$S2NAME".trim_2P.fq.gz --outdir=$OUTDIR
