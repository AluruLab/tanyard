Data Organization:
==================


Files are organized in the folders 

    <EXP_DIR>/CEL/<ACCESS_NAME_DIR>/<CEL_FILES>

QC files are in the folder

    <EXP_DIR>/QC/<ACCESS_NAME_DIR>/<QC_FILES>

Final list of files are in the folders(s)

    <EXP_DIR>/FINAL/<ACCESS_NAME_DIR>/<SOFT_LINK_TO_CEL_FILES>

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

Step 6 - Create links

    ./mkfinal <EXP_DIR>

