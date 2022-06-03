#!/bin/bash

# executable modules
uTR=../uTR
gen_units=gendata/gen
acc=testdata/acc

TR_file="TR_input.txt"
TR_decomp_file="TR_decomp.txt"

num_TRs=1000
error_ratio=0 #0.001

numUnits=2; units="AC AG"; keyUnit="AG";
#numUnits=2; units="ACC GTT"; keyUnit="ACC";
#numUnits=2; units="CCTGGT CTTTT"; keyUnit="CTTTT";
#numUnits=3; units="CCTGGT CTTTT CTTGT"; keyUnit="CTTTT";

#numUnits=2; units="AAG AGG"; keyUnit="AAG";
#numUnits=2; units="AAAG AGG"; keyUnit="AAAG";
#numUnits=3; units="AAAG AGG AAAG"; keyUnit="AAAG";

#numUnits=2; units="AAG AG"; keyUnit="AAG";
#numUnits=2; units="AAAG AG"; keyUnit="AAAG";
#numUnits=3; units="AAAG AG AAAG"; keyUnit="AAAG";
#numUnits=5; units="AAAG AG AGGG AG AAAG"; keyUnit="AAAG";

#numUnits=2; units="AAAGAGAGGGAAAAG AGGGG"; keyUnit="AAAGAGAGGGAAAAG";
#numUnits=3; units="AAAAAG AAAGAGAGGGAAAAG AGGGG"; keyUnit="AAAGAGAGGGAAAAG";

units_name=$units
units_name=${units_name// /\_}

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

        $uTR -f $TR_file -u $keyUnit -std -o $TR_decomp_file 1>> $stat 2>> $wallclock
        # -d print information for debugging

        echo -n $min_unit_occ $max_unit_occ " " >> $acc_table #"Min Max num unit occ"
        $acc $TR_file $TR_decomp_file >> $acc_table
        #diff $TR_file $TR_decomp_file | grep "> >" | wc
    done
done

exit 0
