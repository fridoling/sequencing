#!/bin/bash

usage()
{
  echo "Usage: $0 r1.fastq.gz r2.fastq.gz"
  exit 2
}

if [ $# -lt 2 ]; then
    usage
fi

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

while getopts '?h' c
do
  case $c in
    h|?) usage ;;
  esac
done


# path to sequencing repository
PROJECT="$( cd ../"$(dirname "$0")" ; pwd -P )"

# path to used programs
FASTQC=$PROJECTPATH/bin/bin/FastQC/fastqc
TRIMMO=$PROJECTPATH/bin/Trimmomatic-0.36/trimmomatic-0.36.jar

# path to adapter sequences
ADAPTERS=$PROJECTPATH/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa

# get quality info for unprocessed sequences
$FASTQC $1 $2 



