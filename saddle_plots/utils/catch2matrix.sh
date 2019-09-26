#!/bin/bash

## Author(s): Yinxiu Zhan
## Contact: yinxiu.zhan@fmi.ch
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.

##


function usage {
    echo -e "usage : catch2matrix.sh -i INPUT -b BIN -c CHROM [-n NAME] [-h]"
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
    echo "   -c|--chrom CHROM: chromosome name: "
    echo "   [-n|--name NAME] : name of the sample"
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
      "--chrom") set -- "$@" "-c" ;;
      "--help")   set -- "$@" "-h" ;;
       *)        set -- "$@" "$arg"
  esac
done

INPUT=""
BIN=""
NAME="OUTPUT"
CHROM=""

while getopts ":i:c:b:n:h" OPT
do
    case $OPT in
	i) INPUT=$OPTARG;;
        b) BIN=$OPTARG;;
        n) NAME=$OPTARG;;
	c) CHROM=$OPTARG;;
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



echo "input file: $INPUT"
echo "binsize: $BIN"
echo "output name label: $NAME"

if ! [ -d matrix ]; then
	mkdir matrix
fi

awk '
        BEGIN{
                coord1=0
                coord2=0
		row=0 
		column=0;      
                max=-99
		chr="'"$CHROM"'"
		bin="'"$BIN"'"+0.
        }{
                
                coord1=$2
                coord2=$3
                matr[coord1,coord2]=$4
		matr[coord2,coord1]=$4
                if(coord1<min){
                         min=coord1
                }
                if(coord2<min){
                         min=coord2
                }



                if(coord1>max){
                        max=coord1
                }
                if(coord2>max){
                        max=coord2
                }


	}END{

		print "#""'"$INPUT"'"

		for (i=0; i<=max;i++){ 
			printf("\t""HIC_bin"(i+1)"|dm6|"chr":"(i*bin)"-"((i+1)*bin-1))
			name[i]="HIC_bin"(i+1)"|dm6|"chr":"(i*bin)"-"((i+1)*bin-1)
		}
	
		printf "\n"
		for (row=0; row<=max; row++){
			printf "%s\t", name[row]
			for(column=0; column<=max; column++){
						printf("%f\t", matr[row,column]);
			}
      		     	printf("\n");
    		}

   				
				


	}' $INPUT > $INPUT.matrix.$NAME	

mv $INPUT.matrix.$NAME matrix/
