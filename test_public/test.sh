#!/bin/bash

# generate benchmark data and decompose them into mosaic tandem repeats
bash test_gendata_decompose.sh

# Compute the accuracy of predicting mosaic tandem repeats, say U_i V_j W_k with allowing the values of i, j, and k can differ from the true values at most 1%, 2%, and 3%.
bash test_allowance.sh

exit 0
