
File Formats
------------

Annotation : TAB separated file (.tsv) mapping probe to AGI ids with header PROBE	ID; Multiple ids are delimited by ';'

Gold Standard Network : TAB separated file (tsv) having TF-TARGET relations having header TF,TARGET

Network file formats :
======================
eda format : has a header line; each record is of the format "source (pp) target = weight", where source and target are nodes in the network and weight is the corresponding edge weight.

adj format : adjacency list with a line for each node in the following format:
src tgt1 wt1 tgt2 wt2 ...
where tgt1, tgt2, .. are src's neighbors and wt1, wt2,,

tsv format : tab-seperated file with atleast three columns titled source, target, wt
