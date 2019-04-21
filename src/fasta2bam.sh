#!/bin/bash

usage()
{
  echo "Usage: $0 fasta1 fasta2 [-r REF] [-o OUTDIR]"
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

unset OUTDIR
OPTIND=3
while getopts 'r:o:?h' c
do
  case $c in
    r) set_variable REF $OPTARG ;;
    o) set_variable OUTDIR $OPTARG ;;
    h|?) usage ;;
  esac
done

# if output directory is not given, use current directory 
if [ -z $OUTDIR ]; then
   OUTDIR='.'
fi

# if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]; then
   mkdir -m 775 $OUTDIR
fi

F1NAME=${FASTA1##*/} # get file name
NAME=${F1NAME%%.*}   # remove file extension(s)
OUTPUTALIGN="$OUTDIR"/"$NAME".bam # set output path

# perform alignment, remove duplicates and sort
echo aligning "$NAME"...
{
bwa mem $REF $FASTA1 $FASTA2 |\
samblaster -r |\
samtools view -buS - |\
samtools sort -m 5000000000 -o $OUTPUTALIGN
} 2>"$OUTPUTALIGN".log.txt 

# index bam file
echo indexing "$NAME"...
samtools index $OUTPUTALIGN 2>>"$OUTPUTALIGN".log.txt 


