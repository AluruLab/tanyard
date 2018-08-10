
Input Data::
============

Files are organized in the folders 

    <EXP_DIR>/CEL/<ACCESS_NAME_DIR>/<CEL_FILES>

QC files are in the folder

    <EXP_DIR>/QC/<ACCESS_NAME_DIR>/<QC_FILES>

Final list of files are in the folders(s)

    <EXP_DIR>/FINAL/<ACCESS_NAME_DIR>/<SOFT_LINK_TO_CEL_FILES>


TINGe IQR Steps:
==========================

Step 9 - Merge into array

    ./combine <EXP_DIR> <CSV_SUFFIX> <OUT_CSV>

Step 10 - Quantile normalize

    Rscript quantile.R <IN_CSV> <OUT_CSV>

Step 11 - IQR filtering

    Rscript iqr_filter.R <IN_CSV> <OUT_CSV> <IQR>

Step 12 - Expression file generation

    ./mkexp <IN_CSV>
