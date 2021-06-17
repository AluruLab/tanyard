Workflow Scripts for Processing large-scale Microarray datasets.
=======================================

Scripts for quality control, MAS5 normalization, gene filtering and annotation are present in the
src/qc, src/mas5, src/filter and src/annotation directories. 
A workflow for constructing a gene expression matrix is as follows:

1. Download the datsaets from GEO/ArrayExpress(AE) and organize them as follows:
  `<EXP_DIR>` is the root directory, in which all datasets are placed.
  Files are organized in the folders: `<EXP_DIR>/CEL/<ACCESS_NAME_DIR>/<CEL_FILES>`,
  where `<ACCESS_NAME_DIR>` is the accession id, under which the files were
  submitted in ArrayExpress (AE) or GEO.
2. Run Quality control using the scripts in src/qc directory.
3. Run MAS5 normalization using the scripts in src/mas5 directory.
4. Run IQR Gene filter and quantile normalization using scripts  in src/filter directory.

Please see the README files in respective directories on how to run these scripts. 
To prepare annotation from probe elements, helper scripts are available in the src/annotation directory.


data/ directory contains the annotation, reference networks, classification tables of
A. thaliana datasets used in our paper
