
import pandas as pd
import networkx as nx 
import argparse
from data_utils import load_annotation, map_atid2probes, load_reveng_network, map_probes_cols

def main(annot_file: str, tf_list_file: str, net_file: str, out_file:str):
    annot_df = load_annotation(annot_file)
    tflst_df = map_atid2probes(pd.read_csv(tf_list_file, sep= r'\s+'),
                               annot_df)
    rv_net = load_reveng_network(net_file)
    #rv_net_nodes = set(rv_net.source) | set(rv_net.target)
    rv_net_graph : nx.Graph = nx.from_pandas_edgelist(rv_net, edge_attr='wt')
    #print(rv_net_graph.size())
    subnet_tfs = [x for x in tflst_df.PROBE if x in rv_net_graph]
    subnet_tf_nbrs: set = set([])
    for tfn in subnet_tfs:
        subnet_tf_nbrs.update(list(rv_net_graph.neighbors(tfn)))
    #print(len(subnet_tfs))
    #print(len(subnet_tf_nbrs))
    rv_net_subgraph = rv_net_graph.subgraph(subnet_tf_nbrs|set(subnet_tfs))
    rdf = nx.to_pandas_edgelist(rv_net_subgraph);
    rdf = map_probes_cols(rdf, annot_df, ['source', 'target'])
    rdf.to_csv(out_file, sep="\t");


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser()
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


