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

All scripts take the `<EXP_DIR>` as input; 
The scripts assume that QC is completed and links are created in the 
FINAL directory.

Normalization Steps:
==========================

Step 7 - GCRMA Normalization

    ./gcrma <EXP_DIR>

Step 8 - Log2 mean adjustment

    ./submean2 <FINAL_DIR>

