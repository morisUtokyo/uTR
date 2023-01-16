#!/bin/bash

cd check_uTR; make clean; make;
cd ../gendata; make clean; make;
cd ../parse_RepeatMasker; make clean; make;
cd ../check_RepeatMasker; make clean; make;
cd ../checkTRF; make clean; make;
cd ..

exit 0
