#!/bin/bash
# A bash script file for decomposing a read into a series of units

# The max mismatch rate between units and reads
MAX_DIS_RATIO=0.3
# Specify the path of uTR
uTR=./uTR
# The directory of the input fasta files and summary of repeat locations
inDir=debugAug2021
#inDir=TRsMay21_rep
stat_file=$inDir/hg38_0.8_repeat.txt
# The directory where fasta files for Asano-kun's EDDC algorithm are output
outDir=debugAug2021_result/MAX_DIS_RATIO_$MAX_DIS_RATIO
#outDir=TRsMay21_rep_result
unit_file=$outDir/units.txt

echo -n "" > $unit_file

for fastafile in $( ls $inDir/* | grep .fasta$); do
    # debugAug2021/2001_rep.fasta => 2001
    prefix=${fastafile%_*}  # remove remove the tail until ) (e.g., _rep.fasta)
    number=${prefix##*/}    # remove the head until / (e.g. ZZZZ/)
    
    # debugAug2021/2001.fasta => 2001
    #prefix=${fastafile%.*} # remove the tail until . (e.g. (.fasta)
    #number=${prefix##*/}    # remove the head until / (e.g. ZZZZ/)
    
    info=$(awk '$1=="'"$number"'"{print $1,$2,$3,$4,$8}' < "$stat_file")
    head=${info% *}
    string=${info##* }
    #echo $fastafile $prefix $number $info $head $string
    echo -n $head >> $unit_file # no line breaks
    #echo $uTR -f $fastafile -u $string  >> $unit_file
    $uTR -f $fastafile -u $string -o $outDir/${number}_EDDC.fasta -r $MAX_DIS_RATIO -t >> $unit_file 
done

exit 0
