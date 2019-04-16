import argparse
import pandas as pd
from data_utils import load_full_annotation, load_gsnetwork


def match_columns(annot_df: pd.DataFrame, atrmid: str) -> pd.DataFrame:
    tdf: pd.DataFrame = annot_df.loc[
        annot_df.TID.str.contains(atrmid, regex=False) |
        annot_df.ID.str.contains(atrmid, regex=False) |
        annot_df.AGI.str.contains(atrmid, regex=False) |
        annot_df.UNIGENE.str.contains(atrmid, regex=False) |
        annot_df.ENSEMBL.str.contains(atrmid, regex=False), : ]
    if not tdf.empty:
        tdf['ATRMID'] = [atrmid for _ in range(tdf.shape[0])]
        return tdf.loc[:, ['ATRMID', 'PROBE', 'TID', 'ID', 'AGI', 'ENSEMBL', 'UNIGENE']]
    return pd.DataFrame({
        'ATRMID'  : [atrmid],
        'PROBE'   : ['NO MATCH'],
        'TID'     : ['NO MATCH'],
        'ID'      : ['NO MATCH'],
        'AGI'     : ['NO MATCH'],
        'ENSEMBL' : ['NO MATCH'],
        'UNIGENE' : ['NO MATCH']
    })


def match_network_probes(annot_file: str, gs_file: str) -> pd.DataFrame:
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


def main(annotation_file: str, gs_network_file: str, out_file: str) -> None:
    match_df = match_network_probes(annotation_file,
                                    gs_network_file)
    #print(len(set(match_df.PROBE)))
    match_df.to_csv(out_file, sep='\t', index=False)


if __name__ == "__main__":
    PROG_DESC = """
    Given a gold standard network with AGI ids, find all the probes 
    Network intersection is computed as the intersection of edges of the input networks,
    after mapping to the annotation.
    """
    PARSER = argparse.ArgumentParser(description=PROG_DESC)
    PARSER.add_argument("annotation_file",
                        help="annotation file (a tab seperated file) mapping probe to ids")
    PARSER.add_argument("gs_network_file",
                        help="gold standard network (tab seperated file of TF-TARGET interactions)")
    PARSER.add_argument("-o", "--out_file", type=argparse.FileType('w'), required=True,
                        help="output file in tab-seperated format")
    ARGS = PARSER.parse_args()
    main(ARGS.annotation_file, ARGS.gs_network_file, ARGS.out_file)
