import sys

INPUT_REFSET = {
    "atrm"         : "atrm.tsv",
    "tf2cs"        : "tf2-chipseq.tsv",
    "tf2de"        : "tf2-de.tsv",
    "atreg"        : "AtRegNet-Single-Confirmed.tsv",
    "y1h"          : "Y1H-Reference.tsv",
    "ptregtfbs"    : "PlantRegMap-FnTFBS.tsv",
    "ptregcs"      : "PlantRegMap-Chip-Seq.tsv",
    "ptregmerged"  : "PlantRegMap-Merged.tsv",
    "ptregcslit"   : "PlantRegMap-Chip-Seq-Lit.tsv",
    "ptregmotfbs"  : "PlantRegMap-Motif-TFBS.tsv",
    "n2target"     : "n2target.tsv", 
    "n2targetgrn"  : "n2targetgrn-full.tsv",
    "n2targetgrn3" : "n2targetgrn-3P5.tsv",
    "n2targetgrn9" : "n2targetgrn-9P5.tsv",
    "biogridy1h"   : "BioGRID-Y1H.tsv"
}

INPUT_REFSET_COLS = {
    "atrm"         : (0, 1),
    "tf2cs"        : (0, 1),
    "tf2de"        : (0, 1),
    "atreg"        : (1, 3),
    "y1h"          : (0, 1),
    "ptregtfbs"    : (0, 1),
    "ptregcs"      : (0, 2),
    "ptregmerged"  : (0, 1),
    "ptregcslit"   : (0, 2),
    "ptregmotfbs"  : (0, 2),
    "n2target"     : (0, 1),
    "n2targetgrn"  : (0, 1),
    "n2targetgrn3" : (0, 1),
    "n2targetgrn9" : (0, 1),
    "biogridy1h"   : (0, 1)
}

def merge_refset(edge_set, ref_file, tf_col=0, tgt_col=1):
    with open(ref_file) as fptr:
        fptr.readline()
        for line in fptr:
            line = line.strip()
            elts = [x.strip() for x in line.split("\t")]
            if elts[tf_col] == elts[tgt_col]:
                continue
            if elts[tf_col] < elts[tgt_col]:
                uvx = (elts[tf_col], elts[tgt_col])
            else:
                uvx = (elts[tgt_col], elts[tf_col])
            if uvx not in edge_set:
                edge_set.add(uvx)
                #print("\t".join([elts[tf_col], elts[tgt_col]]))
                print("-".join([uvx[0], uvx[1]]))
    return edge_set

def main(ref_nets):
    #print(ref_nets)
    print("\t".join(["TF", "TARGET"]))
    if "all" in ref_nets:
        ref_nets = list(INPUT_REFSET.keys())
    edge_set = set([])
    for ref_key in ref_nets:
        if ref_key in INPUT_REFSET and ref_key in INPUT_REFSET_COLS:
            tf_col, tgt_col = INPUT_REFSET_COLS[ref_key]
            ref_file = INPUT_REFSET[ref_key]
            edge_set = merge_refset(edge_set, ref_file, tf_col, tgt_col)

if __name__ == "__main__":
    main(set(sys.argv[1:]))
