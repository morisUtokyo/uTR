#!/bin/bash

# executable modules
uTR=../uTR
gen_units=gendata/gen
check_uTR=check_uTR/acc

num_TRs=1000
error_ratio=0

list=("AC_AG" "ACC_GTT" "CCTGGT_CTTTT" "CCTGGT_CTTTT_CTTGT" "AAG_AGG" "AAAG_AGG" "AAG_AG" "AAAG_AG" "AAAG_AG_AAAG" "AAAG_AG_AGGG_AG_AAAG" "AAAAGAAAGAGAGGG_AGGGG" "AGGGG_AAAAGAAAGAGAGGG_AGGGG")

wallclock_uTR="time_uTR.txt"
wallclock_RM="time_RM.txt"
acc_table_uTR="accuracy_uTR.txt"
acc_table_RepeatMasker="accuracy_RepeatMasker.txt"
rm $wallclock $wallclock_RM $acc_table_uTR $acc_table_RepeatMasker

mkdir tmp

for units_name in ${list[@]}
do
    units=${units_name//\_/ }
    echo $units
    echo $units >> $wallclock_uTR
    echo $units >> $wallclock_RM
    echo $units >> $acc_table_uTR
    echo $units >> $acc_table_RepeatMasker

    min_a=(2)
    max_a=(20 50 100 200 500)
    #max_a=(3 5 10 20 50 100 200 500)

    for min_unit_occ in ${min_a[@]}
    do
        for max_unit_occ in ${max_a[@]}
        do
            run_name=$units_name"_"$min_unit_occ"_"$max_unit_occ
            echo $run_name
            param=$min_unit_occ" "$max_unit_occ"\t"
            
            # Generate data of tandem repeats for the given pattern
            TR_file=$run_name".fasta"
            res="${units_name//[^_]}"
            numUnits=$(( ${#res}+1 ))
            #echo $gen_units $min_unit_occ $max_unit_occ $num_TRs $error_ratio $numUnits $units $TR_file
            $gen_units -k $min_unit_occ -l $max_unit_occ -n $num_TRs -e $error_ratio -m $numUnits $units > $TR_file
            # Run uTR
            TR_uTR_decomp_file=$run_name"_decomp_uTR.fasta"
            uTR_stat=$run_name"_unit_stat_uTR.txt"
            
            echo -n -e $param >> $wallclock_uTR
            $uTR -f $TR_file -std -o $TR_uTR_decomp_file 1>> $uTR_stat 2>> $wallclock_uTR
            rm $uTR_stat
            # Compute the accuracy of uTR
            echo -n -e $param >> $acc_table_uTR #"Min Max num unit occ"
            $check_uTR $TR_file $TR_uTR_decomp_file >> $acc_table_uTR
            
            # Run RepeatMasker
            parse_RM=$run_name".fasta.result.txt"
            match_RM=$run_name".fasta.match.txt"
            
            # Apply RepeatMasker to the string set
            SECONDS=0
            repeatmasker -e hmmer -noint -pa 4 -div 0 -xsmall $TR_file &> progress_report.txt
            run_time=$SECONDS
            echo -n -e $param >> $wallclock_RM
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
            echo -n -e $param >> $acc_table_RepeatMasker
            tail -n1 $match_RM >> $acc_table_RepeatMasker
            
            mv $TR_file  tmp
            mv $TR_uTR_decomp_file tmp
            mv $parse_RM tmp
            mv $match_RM tmp
        done
    done
done

exit 0
