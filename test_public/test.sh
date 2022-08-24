#!/bin/bash

# executable modules
uTR=../uTR
gen_units=gendata/gen
acc=testdata/acc

num_TRs=1000
error_ratio=0

list=("AC_AG" "ACC_GTT" "CCTGGT_CTTTT" "CCTGGT_CTTTT_CTTGT" "AAG_AGG" "AAAG_AGG" "AAG_AG" "AAAG_AG" "AAAG_AG_AAAG" "AAAG_AG_AGGG_AG_AAAG" "AAAAGAAAGAGAGGG_AGGGG" "AGGGG_AAAAGAAAGAGAGGG_AGGGG")

for units_name in ${list[@]}
do
    units=${units_name//\_/ }
    echo $units

    res="${units_name//[^_]}"
    numUnits=$(( ${#res}+1 ))

    TR_file=$units_name".fasta"
    TR_decomp_file=$units_name"_decomp.fasta"
    #echo $TR_file
    
    stat=$units_name"_unit_stat.txt"
    wallclock=$units_name"_time.txt"
    acc_table=$units_name"_accuracy.txt"
    rm $stat
    rm $wallclock
    rm $acc_table
    echo "Units       =" $units > $acc_table
    echo "Num. TRs    = " $num_TRs   >> $acc_table
    echo "Error ratio = " $error_ratio   >> $acc_table

    min_a=(2)
    max_a=(3 5 10 20 50 100 200 500)

    for min_unit_occ in ${min_a[@]}
    do
        for max_unit_occ in ${max_a[@]}
        do
            $gen_units -k $min_unit_occ -l $max_unit_occ -n $num_TRs -e $error_ratio -m $numUnits $units > $TR_file
            
            $uTR -f $TR_file -std -o $TR_decomp_file 1>> $stat 2>> $wallclock

            echo -n $min_unit_occ $max_unit_occ " " >> $acc_table #"Min Max num unit occ"
            $acc $TR_file $TR_decomp_file >> $acc_table
        done
    done
done

exit 0
