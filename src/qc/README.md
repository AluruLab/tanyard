Data Organization:
==================

`<EXP_DIR>` is the root directory, in which all datasets are placed.

Files are organized in the folders 

    <EXP_DIR>/CEL/<ACCESS_NAME_DIR>/<CEL_FILES>

where `<ACCESS_NAME_DIR>` is the accession id, under which the files were
submitted in ArrayExpress (AE) or GEO.

QC files are generated for each submission in AE/GEO in the folder:

    <EXP_DIR>/QC/<ACCESS_NAME_DIR>/<QC_FILES>

Final list of QC-passed CEL files are softlinks and are placed 
in the folders(s):

    <EXP_DIR>/FINAL/<ACCESS_NAME_DIR>/<SOFT_LINK_TO_CEL_FILES>

Micorarray QC Steps:
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

