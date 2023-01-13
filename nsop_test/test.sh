#!/bin/bash

echo "----------------------------"
echo "Compile gendata"
cd gendata
make clean; make

echo "----------------------------"
echo "Compile uTR"
cd ../uTR
make clean; make
cd ..

echo "----------------------------"
echo "Calculate compression ratio"
gen_units=gendata/gen
uTR=uTR/uTR

num_TRs=1000
min_unit_occ=10
max_unit_occ=200

result="nsop_compression.csv"
rm $result


listError=(0.0 0.01 0.03 0.05 0.1 0.15)
listPattern=("AC_AG" "ACC_GTT" "AAG_AG" "AAG_AGG" "AAAG_AG" "AAAG_AG_AAAG" "AAAG_AG_AGGG_AG_AAAG" "AGGGG_AAAAGAAAGAGAGGG_AGGGG")

for error_ratio in ${listError[@]}
    do
    
    for units_name in ${listPattern[@]}
    do
        units=${units_name//\_/ }
        run_name=$units_name"_"$min_unit_occ"_"$max_unit_occ"_"$error_ratio
            
        # Generate data of tandem repeats for the given pattern
        TR_file=$run_name".fasta"
        res="${units_name//[^_]}"
        numUnits=$(( ${#res}+1 ))

        $gen_units -k $min_unit_occ -l $max_unit_occ -n $num_TRs -e $error_ratio -m $numUnits $units > $TR_file
        
        echo -n $run_name >> $result
        $uTR -f $TR_file >> $result
        rm $TR_file

    done
done

exit 0
