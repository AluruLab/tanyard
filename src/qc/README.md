Input Data::
============

Files are organized in the folder structure 

    <EXP_DIR>/CEL/<ACCESS_NAME_DIR>/<CEL_FILES>

QC output is written to the folder

    <EXP_DIR>/QC/<ACCESS_NAME_DIR>/<QC_FILES>

All scripts take the <EXP_DIR> as input

TINGe QC Steps:
===============

Step 1 - Affymetrix QC

    ./affyqc <EXP_DIR>

Step 2 - RLE + NUSE

    ./plmstat <EXP_DIR>

Step 3 - Deleted residuals

    ./delres <EXP_DIR>

Step 4 - Summary

     ./summary <EXP_DIR>

Step 5 - Report from summary

    ./report <EXP_DIR>

TINGe Normalization Steps:
==========================

Step 6 - Create links
./links <EXP_DIR>

Step 7 - MAS-5
cd .FINAL
../mas5

Step 8 - Log2 mean adjustment
cd .FINAL
../submean2

Step 9 - Merge into array
cd .FINAL
../combine <OUT_CSV>

TINGe IQR Steps:
==========================

Step 10 - Quantile normalize
R_IN=<IN_CSV> R_OUT=<OUT_CSV> R < quantile.R

Step 11 - IQR filtering
R_IN=<IN_CSV> R_OUT=<OUT_CSV> R_IQR=<IQR> R < iqr-filter.R

Step 12 - Expression file generation
./mkexp <IN_CSV>
