Input Data::
============

Files are organized in the folders 

    <EXP_DIR>/CEL/<ACCESS_NAME_DIR>/<CEL_FILES>

QC files are in the folder

    <EXP_DIR>/QC/<ACCESS_NAME_DIR>/<QC_FILES>

Final list of files are in the folders(s)

    <EXP_DIR>/FINAL/<ACCESS_NAME_DIR>/<SOFT_LINK_TO_CEL_FILES>

All scripts take the <EXP_DIR> as input; 
Assumes that QC is completed and links are created.

TINGe Normalization Steps:
==========================

Step 2 - MAS-5

    ./gcrma <EXP_DIR>

Step 3 - Log2 mean adjustment

    ./submean2 <FINAL_DIR>

