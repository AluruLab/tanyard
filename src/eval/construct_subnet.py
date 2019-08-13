
import argparse
from typing import Set
import pandas as pd
import networkx as nx
import data_utils as du

def main(annot_file: str, tf_list_file: str, net_file: str, out_file: str):
    annot_df = du.load_annotation_alias(annot_file)
    tflst_df = du.map_atid2probes(pd.read_csv(tf_list_file, sep=r'\s+'),
                                  annot_df)
    rv_net = du.load_reveng_network(net_file)
    #rv_net_nodes = set(rv_net.source) | set(rv_net.target)
    rv_net_graph: nx.Graph = nx.from_pandas_edgelist(rv_net, edge_attr='wt')
    #print(rv_net_graph.size())
    subnet_tfs = [x for x in tflst_df.PROBE if x in rv_net_graph]
    subnet_tf_nbrs: Set[str] = set([])
    for tfn in subnet_tfs:
        subnet_tf_nbrs.update(list(rv_net_graph.neighbors(tfn)))
    #print(len(subnet_tfs))
    print(len(subnet_tf_nbrs))
    rv_net_subgraph = rv_net_graph.subgraph(subnet_tf_nbrs|set(subnet_tfs))
    rdf: pd.DataFrame = nx.to_pandas_edgelist(rv_net_subgraph)
    print(annot_df.columns)
    if 'ALIAS' in annot_df.columns:
        rdf = du.map_probes_cols_idalias(rdf, annot_df, ['source', 'target'])
        rdf = rdf.loc[:, ['source_id', 'target_id',
                          'source_alias', 'target_alias', 'wt']]
    else:
        rdf = du.map_probes_cols(rdf, annot_df, ['source', 'target'])
        rdf = rdf.loc[:, ['source_id', 'target_id', 'wt']]
    rdf = rdf.loc[rdf.source_id != rdf.target_id, :]
    rdf['source_ind'] = rdf.source_id.isin(tflst_df.ID)
    rdf['target_ind'] = rdf.target_id.isin(tflst_df.ID)
    rdf.loc[rdf['source_ind'] | rdf['target_ind'],
            :].to_csv(out_file, sep="\t", index=False)


if __name__ == "__main__":
    PROG_DESC = """
    Given an input network and transcription factors, finds the 
    and common to all the networks.
    """
    PARSER = argparse.ArgumentParser(description=PROG_DESC)
    PARSER.add_argument("annotation_file",
                        help="""annotation file
                                (a tab seperated file mapping probe to ids)""")
    PARSER.add_argument("tf_list_file",
                        help="""list of transcription factors""")
    PARSER.add_argument("reveng_network_file",
                        help="""network build from a reverse engineering methods
                                (currenlty supported: eda, adj, tsv)""")
    PARSER.add_argument("-o", "--out_file",
                        type=str,
                        help="output file in png format")
    ARGS = PARSER.parse_args()
    main(ARGS.annotation_file, ARGS.tf_list_file,
         ARGS.reveng_network_file, ARGS.out_file)
