#!/bin/bash

# executable modules
check_uTR=check_uTR/acc
checkRM=check_RepeatMasker/checkRM
checkTRF=checkTRF/checkTRF
tmp=tmp
#tmp=tmp_backup

num_TRs=1000

listError=(0.0 0.01 0.03 0.05 0.1 0.15)
listAllowance=(0 0.01 0.02 0.03)
listPattern=("AC_AG" "ACC_GTT" "AAG_AG" "AAG_AGG" "AAAG_AG" "AAAG_AG_AAAG" "AAAG_AG_AGGG_AG_AAAG" "AGGGG_AAAAGAAAGAGAGGG_AGGGG")

for allowance in ${listAllowance[@]}
do
    acc_table_uTR=$tmp"/accuracy_uTR_allowance"$allowance".txt"
    acc_table_RM=$tmp"/accuracy_RM_allowance"$allowance".txt"
    acc_table_TRF=$tmp"/accuracy_TRF_allowance"$allowance".txt"
    rm $acc_table_uTR $acc_table_RM $acc_table_TRF
    
    for error_ratio in ${listError[@]}
    do
        for units_name in ${listPattern[@]}
        do
            units=${units_name//\_/ }
            echo $units
            
            min_a=(10) 
            max_a=(200)
            for min_unit_occ in ${min_a[@]}
            do
                for max_unit_occ in ${max_a[@]}
                do
                    run_name=$units_name"_"$min_unit_occ"_"$max_unit_occ"_"$error_ratio
                    echo $run_name

                    # uTR
                    TR_uTR_decomp_file=$tmp"/"$run_name"_decomp_uTR.fasta"
                    echo -n -e $run_name" " >> $acc_table_uTR #"Min Max num unit occ"
                    $check_uTR -i $TR_uTR_decomp_file -a $allowance >> $acc_table_uTR
                    
                    # RepeatMasker
                    parse_RM=$tmp"/"$run_name".fasta.result.txt"
                    match_RM=$tmp"/"$run_name".fasta.match.txt"
                    $checkRM -i $parse_RM -a $allowance -p > $match_RM
                    echo -n -e $run_name" " >> $acc_table_RM
                    tail -n1 $match_RM >> $acc_table_RM
                    
                    # TRF
                    result_TRF=$tmp"/"$run_name"_TRF.txt"
                    accuracy_TRF=$tmp"/"$run_name"_TRF_acc.txt"
                    $checkTRF -i $result_TRF -a $allowance -o $accuracy_TRF
                    echo -n -e $run_name" " >> $acc_table_TRF
                    tail -n1 $accuracy_TRF >> $acc_table_TRF
                done
            done
        done
    done
done

exit 0
