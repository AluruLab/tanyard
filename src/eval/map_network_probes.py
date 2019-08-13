import argparse
from typing import List
import data_utils as du

def main(annot_file: str, net_file: str, net_col_names: List[str], out_file: str):
    annot_df = du.load_annotation_alias(annot_file)
    rv_net = du.load_reveng_network(net_file)
    net_col_ids = [x+'_id' for x in net_col_names]
    net_col_alias = [x+'_alias' for x in net_col_names]
    if 'ALIAS' in annot_df.columns:
        rdf = du.map_probes_cols_idalias(rv_net, annot_df, net_col_names)
        rdf = rdf.loc[:, net_col_ids + net_col_alias + ['wt']]
    else:
        rdf = du.map_probes_cols(rv_net, annot_df, net_col_names)
        rdf = rdf.loc[:, net_col_ids + ['wt']]
    rdf.to_csv(out_file, sep="\t", index=False)


if __name__ == "__main__":
    PROG_DESC = """ Map the network columns to id and alias"""
    PARSER = argparse.ArgumentParser(description=PROG_DESC)
    PARSER.add_argument("annotation_file",
                        help="annotation file (a tab seperated file) mapping probe to ids")
    PARSER.add_argument("reveng_network_file",
                        help="""network build from a reverse engineering methods
                                (currenlty supported: eda, adj, tsv)""")
    PARSER.add_argument("-c", "--network_cols", type=str, default="source,target",
                        help="""comma seperated names of the output files;
                                should have as many names as the number of networks""")
    PARSER.add_argument("-o", "--out_file", type=argparse.FileType('w'), required=True,
                        help="output file in tab-seperated format")
    ARGS = PARSER.parse_args()
    MAP_COL_NAMES = ARGS.network_cols.split(",")
    main(ARGS.annotation_file, ARGS.reveng_network_file, MAP_COL_NAMES, ARGS.out_file)
