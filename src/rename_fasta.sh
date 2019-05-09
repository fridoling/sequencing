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
    # check if file name already starts with PREFIX 
    if [[ $filename != "$PREFIX"* ]]; then
	newfilename="$PREFIX"_$filename
	newfile=$INPUTDIR/$newfilename
	# check if the new filename is already taken
	if [ -f $newfile ]; then
	    echo Skip $filename: $newfilename already exists
	else
	    mv $file $newfile
	fi
    fi
done

