import pandas as pd
import argparse
from data_utils import load_full_annotation, load_gsnetwork

def match_columns(annot_df, atrmid):
    tdf = annot_df.loc[
            annot_df.TID.str.contains(atrmid, regex=False) |
            annot_df.ID.str.contains(atrmid, regex=False) |
            annot_df.AGI.str.contains(atrmid, regex=False) |
            annot_df.UNIGENE.str.contains(atrmid, regex=False) |
            annot_df.ENSEMBL.str.contains(atrmid, regex=False) , : ]
    if not tdf.empty:
        tdf['ATRMID'] = [ atrmid for _ in range(tdf.shape[0])]
        return tdf.loc[:, ['ATRMID', 'PROBE', 'TID', 'ID', 'AGI', 'ENSEMBL', 'UNIGENE']]
    else:
        return pd.DataFrame({
            'ATRMID'  : [atrmid],
            'PROBE'   : ['NO MATCH'],
            'TID'     : ['NO MATCH'],
            'ID'      : ['NO MATCH'],
            'AGI'     : ['NO MATCH'],
            'ENSEMBL' : ['NO MATCH'],
            'UNIGENE' : ['NO MATCH']
        })


def match_network_probes(annot_file, gs_file):
    #def match_df(name)
    annot_df = load_full_annotation(annot_file)
    gsnet_df = load_gsnetwork(gs_file)
    gs_nodes = list(set(gsnet_df.TF) | set(gsnet_df.TARGET))
    gs_nodes_df = pd.DataFrame({'ATRMID' : gs_nodes})
    gsnet_grouped = gs_nodes_df.groupby('ATRMID')
    return gsnet_grouped.apply(lambda x: match_columns(annot_df, x.name))
    #print(annot_df.columns)
    #print(pd.DataFrame(fdf))
    #rdf = grouped.apply()


if __name__ == "__main__":
    prog_desc = """
    Given a gold standard network with AGI ids, 
    Network intersection is computed as the intersection of edges of the input networks,
    after mapping to the annotation.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("annotation_file",
        help="annotation file (a tab seperated file) mapping probe to ids")
    parser.add_argument("gs_network_file",
        help="gold standard network (tab seperated file of TF-TARGET interactions)")
    parser.add_argument("-o", "--out_file", type=argparse.FileType('w'), required=True,
                        help="output file in tab-seperated format")
    args = parser.parse_args()
    match_df = match_network_probes(args.annotation_file,
                                    args.gs_network_file)
    #print(len(set(match_df.PROBE)))
    match_df.to_csv(args.out_file, sep='\t', index=False)
