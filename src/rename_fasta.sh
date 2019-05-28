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
for file in $INPUTDIR/*[1-2].fq.gz
do
    filename=$(basename "$file")

#    # SAMPLE NAME
#    if [[ $filename == "$PREFIX"* ]] ; then
#	# if it already has PREFIX, do nothing
#	continue	
#    elif [[ $filename =~ _F ]]; then
#	# if it is NOVOGENE 190516 name, take out SAMPLENAME
#	# sample name is what's before the first _F(letter)
#	samplename=${filename%%_F[A-Z]*}
#	newfilename="$PREFIX"_"$samplename"
#   fi

    # SAMPLE NAME (if they already have a prefix)
    if [[ $filename =~ _F ]]; then
	# if it is NOVOGENE 190516 name, take out SAMPLENAME
	# sample name is what's before the first _F(letter)
	samplename=${filename%%_F[A-Z]*}
	newfilename="$samplename"
    fi
    
    # SUFFIX
    if [[ $filename == *"_1.fq.gz" ]] ; then
	SUFFIX="_1.fq.gz"	
    elif [[ $filename == *"_2.fq.gz" ]] ; then
	SUFFIX="_2.fq.gz"		
    fi    
    newfilename=$newfilename$SUFFIX
    
    newfile=$INPUTDIR/$newfilename
    # check if the new filename is already taken
    if [ -f $newfile ]; then
	echo Skip $filename: $newfilename already exists
    else
	mv $file $newfile
    fi
done

