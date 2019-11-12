#!/bin/bash

## Author(s): Camilla Ugolini
## Contact: camilla.ugolini@fmi.ch
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.

##


function usage {
    echo -e "usage : remove_regions_from_HiC.sh -i INPUT -b BIN  [-n NAME] [-h]"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo 
    echo "Bin data with exponential bins"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -i|--input INPUT:  *.CaTCH* file output of hicpro2catch.sh "
    echo "   -b|--binsize BINSIZE The binning size of matrix: "
    echo "   -r|--remove REMOVE: file containing the list of region to remove"
    echo "   [-n|--name NAME] : name of the sample, default inputfilename.filtered"
    echo "   [-h|--help]: help"
    exit;
}

# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
      "--input") set -- "$@" "-i" ;;
      "--binsize")   set -- "$@" "-b" ;;
      "--name") set -- "$@" "-n" ;;
      "--remove") set -- "$@" "-r" ;;
      "--help")   set -- "$@" "-h" ;;
       *)        set -- "$@" "$arg"
  esac
done

INPUT=""
BIN=""
NAME="filtered"
REMOVE=""

while getopts ":i:r:b:n:h" OPT
do
    case $OPT in
	i) INPUT=$OPTARG;;
        b) BIN=$OPTARG;;
        n) NAME=$OPTARG;;
	r) REMOVE=$OPTARG;;
        h) help ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            exit 1
            ;;
    esac
done


if [ $# -lt 6 ]
then
    usage
    exit
fi



if ! [ -f $INPUT ]; then

	echo "$INPUT file does not exist"
	exit

fi


if ! [ -f $REMOVE ]; then

        echo "$REMOVE file does not exist"
        exit

fi



echo "input file: $INPUT"
echo "binsize: $BIN"
echo "output name label: $NAME"
echo "file of region to remove: $REMOVE"



awk 'BEGIN{fn=0; r=0; b="'"$BIN"'"+0.}{

	if(FNR==1) fn++
	if(fn==1){
		chr[r]=$1
		s[r]=int($2/b)
		e[r]=int($3/b)+1
		r++
	}


	if(fn==2){
		bool=-1
		for(i=0;i<r;i++){
			if($1==chr[i]){
				if(($2>=s[i] && $2<=e[i]) || ($3>=s[i] && $3<=e[i])){
					bool=1
				}
				break
			}
		}
		
		if(bool==-1) print $0

		
	}

}' $REMOVE $INPUT > $INPUT.$NAME
