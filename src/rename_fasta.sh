#!/bin/bash

usage()
{
  echo "Usage: ${0##*/} inputdir prefix "
  exit 2
}

if [ $# -lt 2 ]; then
    usage
fi

INPUTDIR=$1
PREFIX=$2

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

# check whether input directory exist
if [ ! -d $INPUTDIR ]; then
   echo "Error: Please provide paths to existing folder." 
   exit 2
fi



# rename files in INPUTDIR
for file in $INPUTDIR/*.gz
do
    filename=$(basename "$file")
    # if it is illumina name, take out SAMPLENAME
    if [[ $filename =~ _S[0-9]{1,2}_L[0-9]{3}_R[1-2]_001 ]]; then
	# sample name is what's before the first underscore
	samplename=${filename%%_*}
	newfilename="$PREFIX"_"$samplename".fastq.gz
    else if [[ $filename =~ _F ]]; then
	# if it is NOVOGENE 190516 name, take out SAMPLENAME
	# sample name is what's before the first _F
	samplename=${filename%%_F*}
	newfilename="$PREFIX"_"$samplename".fq.gz
	# otherwise simply add PREFIX to the name
	newfilename="$PREFIX"_$filename
    fi
    
    newfile=$INPUTDIR/$newfilename
    # check if the new filename is already taken
    if [ -f $newfile ]; then
	echo Skip $filename: $newfilename already exists
    else
	mv $file $newfile
    fi

done

