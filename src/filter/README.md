
Input Data::
============

`<EXP_DIR>` is the root directory, in which all datasets are placed.

Files are organized in the folders as below:

    <EXP_DIR>/CEL/<ACCESS_NAME_DIR>/<CEL_FILES>

where `<ACCESS_NAME_DIR>` is the accession id, under which the files were
submitted in ArrayExpress (AE) or GEO.

QC files are already generated from the QC step,
for each submission in AE/GEO, in the folder:

    <EXP_DIR>/QC/<ACCESS_NAME_DIR>/<QC_FILES>

Final list of QC-passed CEL files are softlinks and are expected to be
in the folders(s):

    <EXP_DIR>/FINAL/<ACCESS_NAME_DIR>/<SOFT_LINK_TO_CEL_FILES>



IQR Gene Filter Steps:
==========================

Step 9 - Merge into array

    ./combine <EXP_DIR> <CSV_SUFFIX> <OUT_CSV>

`<CSV_SUFFIX>` should be mas5 if normalization was done using MAS5, gcrma if normalization done with GCRMA.

Step 10 - Quantile normalize

    Rscript quantile.R <IN_CSV> <OUT_CSV>

Step 11 - IQR filtering

    Rscript iqr_filter.R <IN_CSV> <OUT_CSV> <IQR>

Step 12 - Build exp-clean-ht file

    make

Step 13 - Expression file generation

    ./mkexp <IN_CSV>
