#!/bin/bash
# executable modules
uTR=../uTR
uTR_param="-stda"  

gen_units=gendata/gen
check_uTR=check_uTR/acc
checkTRF=checkTRF/checkTRF
tmp=tmp
mkdir $tmp
num_TRs=1000

listError=(0.0 0.01 0.03 0.05 0.1 0.15)

listPattern=("AC_AG" "ACC_GTT" "AAG_AG" "AAG_AGG" "AAAG_AG" "AAAG_AG_AAAG" "AAAG_AG_AGGG_AG_AAAG" "AGGGG_AAAAGAAAGAGAGGG_AGGGG")

wallclock_uTR="time_uTR.txt"
wallclock_RM="time_RM.txt"
wallclock_TRF="time_TRF.txt"

acc_table_uTR="accuracy_uTR.txt"
acc_table_RM="accuracy_RM.txt"
acc_table_TRF="accuracy_TRF.txt"

rm $wallclock_uTR $wallclock_RM $wallclock_TRF
rm $acc_table_uTR $acc_table_RM $acc_table_TRF

for error_ratio in ${listError[@]}
do
    for units_name in ${listPattern[@]}
    do
        units=${units_name//\_/ }
        min_a=(10)
        max_a=(200)
        for min_unit_occ in ${min_a[@]}
        do
            for max_unit_occ in ${max_a[@]}
            do
                run_name=$units_name"_"$min_unit_occ"_"$max_unit_occ"_"$error_ratio
                echo $run_name
                
                # Generate data of tandem repeats for the given pattern
                TR_file=$run_name".fasta"
                res="${units_name//[^_]}"
                numUnits=$(( ${#res}+1 ))
                $gen_units -k $min_unit_occ -l $max_unit_occ -n $num_TRs -e $error_ratio -m $numUnits $units > $TR_file
                
                # Run uTR
                TR_uTR_decomp_file=$run_name"_decomp_uTR.fasta"
                uTR_stat=$run_name"_unit_stat_uTR.txt"
                $uTR -f $TR_file $uTR_param -o $TR_uTR_decomp_file 1>> $uTR_stat 2>> $wallclock_uTR
                rm $uTR_stat
                # Compute the accuracy of uTR
                echo -n -e $run_name" " >> $acc_table_uTR
                $check_uTR -i $TR_uTR_decomp_file >> $acc_table_uTR
                
                # Run RepeatMasker
                parse_RM=$run_name".fasta.result.txt"
                match_RM=$run_name".fasta.match.txt"
                # Apply RepeatMasker to the string set
                SECONDS=0
                repeatmasker -e hmmer -noint -pa 4 -div 0 -xsmall $TR_file &> progress_report.txt
                run_time=$SECONDS
                echo -n -e $run_name" " >> $wallclock_RM
                echo $run_time >> $wallclock_RM
                #repeatmasker -e hmmer -noint -pa 4 -div 20 -xsmall $TR_file &> progress_report.txt
                # -e(ngine) [crossmatch|wublast|abblast|ncbi|rmblast|hmmer]
                # -noint Only masks low complex/simple repeats (no interspersed repeats)
                # -pa 4 The number of sequence batch jobs [50kb minimum] to run in parallel.
                # -div [number] Masks only those repeats < x percent diverged from consensus seq
                # -xsmall Returns repetitive regions in lowercase (rest capitals) rather than masked
                parse_RepeatMasker/parseRM -i $TR_file.out -o $parse_RM
                # Parse a repeatmasker .out file and output the result
                # Remove temporary files
                rm $TR_file.cat
                rm $TR_file.tbl
                rm $TR_file.out
                rm $TR_file.masked
                rm progress_report.txt
                # Parse the file to determine how RepeatMasker could identify mosaic tandem repeats
                check_RepeatMasker/checkRM -i $parse_RM -p > $match_RM
                # Print out the statistics
                echo -n -e $run_name" " >> $acc_table_RM
                tail -n1 $match_RM >> $acc_table_RM
                
                # Run TRF
                result_TRF=$run_name"_TRF.txt"
                accuracy_TRF=$run_name"_TRF_acc.txt"
                # Apply TRF to the string set
                SECONDS=0
                trf $TR_file 2 7 7 80 10 10 1000 -h -ngs > $result_TRF
                #trf $TR_file 2 7 7 80 10 50 500 -h -ngs > $result_TRF
                run_time=$SECONDS
                echo -n -e $run_name" " >> $wallclock_TRF
                echo $run_time >> $wallclock_TRF
                $checkTRF -i $result_TRF -o $accuracy_TRF
                # Print out the statistics
                echo -n -e $run_name" " >> $acc_table_TRF
                tail -n1 $accuracy_TRF >> $acc_table_TRF
                
                mv $TR_file       $tmp
                mv $TR_uTR_decomp_file $tmp
                mv $parse_RM      $tmp
                mv $match_RM      $tmp
                mv $result_TRF    $tmp
                mv $accuracy_TRF  $tmp
            done
        done
    done
done

mv $wallclock_uTR $tmp
mv $wallclock_RM  $tmp
mv $wallclock_TRF $tmp
mv $acc_table_uTR $tmp
mv $acc_table_RM $tmp
mv $acc_table_TRF $tmp

exit 0
