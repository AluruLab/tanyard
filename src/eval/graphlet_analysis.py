from typing import List
import argparse
import networkx as nx
import pandas as pd
import data_utils as du


def count_common_tf_graphlets(tf_list: List[str], tf_subnet: nx.Graph):
    """
    Graphlets of type 1 :
       1.  Y --- TF1 ---  X
           A common TF for two non-TFs
    """
    ctf_graphlet_ctx = {}
    tf_list_set = set(tf_list)
    for tfx in tf_list:
        nbx = sum((1 if y in tf_list_set else 0 for y in tf_subnet.neighbors(tfx)))
        if nbx > 1:
            ctf_graphlet_ctx[tfx] = nbx * (nbx - 1)
        else:
            ctf_graphlet_ctx[tfx] = 0
    return ctf_graphlet_ctx

def count_commtgt_graphlets(tf_list: List[str], tf_subnet: nx.Graph):
    """
    Graphlets of type 2:
        2. TF1 --- X --- TF2
            A non-TF regulated by two TFs
    """
    ctg_comtgt_ctx = {}
    for tfx in tf_list:
        ctg_comtgt_ctx[tfx] = 0
    tf_list_set: set = set(tf_list)
    if len(tf_list) < 2:
        return ctg_comtgt_ctx
    tf_nbrs = {x : set(list(tf_subnet.neighbors(x))) for x in tf_list}
    for tfx in tf_list:
        xnbrs = tf_nbrs[tfx]
        for tfy in tf_list[1:]:
            ynbrs = tf_nbrs[tfy]
            nzlen = len(tf_list_set - (xnbrs & ynbrs))
            ctg_comtgt_ctx[tfx] += nzlen
            ctg_comtgt_ctx[tfy] += nzlen
    return ctg_comtgt_ctx


def count_passtf_graphlets(tf_list: List[str], tf_subnet: nx.Graph):
    """
    Graphlets of type 2:
        2. TF2 --- TF1 --- X
            A non-TF regulated by a TF, but not by another
    """
    ctg_passtf_ctx = {}
    for tfx in tf_list:
        ctg_passtf_ctx[tfx] = 0
    tf_list_set: set = set(tf_list)
    if len(tf_list) < 2:
        return ctg_passtf_ctx
    tf_nbrs = {x : set(list(tf_subnet.neighbors(x))) for x in tf_list}
    for tfx in tf_list:
        xnbrs = tf_nbrs[tfx]
        tfnbrs = tf_list_set & xnbrs
        ntfnbrs = xnbrs - tfnbrs
        ctg_passtf_ctx[tfx] += len(tfnbrs) * len(ntfnbrs)
    return ctg_passtf_ctx


def count_triangle_tf_graphlets(tf_list: List[str], tf_subnet: nx.Graph):
    """
    Graphlets of type 4 :
        4. TF1 --- TF2 --- TF3
            |               |
            |---------------|
            A triangle of TFs
    """
    ctf_triangle_ctx = {}
    for tfx in tf_list:
        ctf_triangle_ctx[tfx] = 0
    if len(tf_list) < 3:
        return ctf_triangle_ctx
    for tfx in tf_list:
        for tfy in tf_list[1:]:
            if not tf_subnet.has_edge(tfx, tfy):
                continue
            for tfz in tf_list[2:]:
                if not tf_subnet.has_edge(tfx, tfz):
                    continue
                if not tf_subnet.has_edge(tfy, tfz):
                    continue
                ctf_triangle_ctx[tfx] += 1
                ctf_triangle_ctx[tfy] += 1
                ctf_triangle_ctx[tfz] += 1
    return ctf_triangle_ctx

def count_network_graphlets(tflst_df: pd.DataFrame,
                            net_file: str):
    rv_net = du.load_reveng_network(net_file)
    #rv_net_nodes = set(rv_net.source) | set(rv_net.target)
    rv_net_graph: nx.Graph = nx.from_pandas_edgelist(rv_net, edge_attr='wt')
    #print(rv_net_graph.size())
    subnet_tfs = [x for x in tflst_df.PROBE if x in rv_net_graph]
    subnet_tf_nbrs: set = set([])
    for tfn in subnet_tfs:
        subnet_tf_nbrs.update(list(rv_net_graph.neighbors(tfn)))
    #print(len(subnet_tfs))
    #print(len(subnet_tf_nbrs))
    rv_net_subgraph = rv_net_graph.subgraph(subnet_tf_nbrs|set(subnet_tfs))
    ctf_cts = count_common_tf_graphlets(subnet_tfs, rv_net_subgraph)
    tgl_cts = count_triangle_tf_graphlets(subnet_tfs, rv_net_subgraph)
    ctg_cts = count_commtgt_graphlets(subnet_tfs, rv_net_subgraph)
    ptf_cts = count_passtf_graphlets(subnet_tfs, rv_net_subgraph)
    average_graphlets = [{'WT': 0.25 * (ctf_cts[x] + tgl_cts[x] + ctg_cts[x] + ptf_cts[x]),
                          'PROBE': x}
                         for x in subnet_tfs]
    #sorted(average_graphlets, reverse=True)
    rdf: pd.DataFrame = pd.DataFrame.from_records(average_graphlets)
    rdf = du.map_probes2atid(rdf, tflst_df)
    rdf.sort_values(by='WT', ascending=False, inplace=True)
    print(rdf.columns)
    return rdf


def main(annot_file: str, tf_list_file: str, net_files: str, out_file: str):
    annot_df = du.load_annotation(annot_file)
    tflst_df = du.map_atid2probes(pd.read_csv(tf_list_file, sep=r'\s+'),
                                  annot_df)
    for nx_file in net_files:
        rdf: pd.DataFrame = count_network_graphlets(tflst_df, nx_file)
        rdf.to_csv(out_file, sep="\t")
        break


# Algorithm:
# Four graphlet types
#   1.  Y --- TF1 ---  X
#       A common TF for two non-TFs
#   2. TF1 --- X --- TF2
#       A non-TF regulated by two TFs
#   3. TF2 --- TF1 --- X
#       A non-TF regulated by a TF, but not by another
#   4. TF1 --- TF2 --- TF3
#       |               |
#       |---------------|
#       A triangle of TFs
# For each TF :
#   For each of the four graphlet
#      - Count the number of graphlets found
#      - Compute REC for the graphlet
#   Compute RGD of the TF

if __name__ == "__main__":
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("annotation_file",
                        help="""annotation file
                                (a tab seperated file mapping probe to ids)""")
    PARSER.add_argument("tf_list_file",
                        help="""list of transcription factors""")
    PARSER.add_argument("reveng_network_files", nargs="+",
                        help="""network build from a reverse engineering methods
                                (currenlty supported: eda, adj, tsv)""")
    PARSER.add_argument("-n", "--dataset_names", type=str,
                        help="""comma seperated names of the network;
                                should have as many names as the number of networks""")
    PARSER.add_argument("-o", "--out_file",
                        type=str,
                        help="output file in png format")
    ARGS = PARSER.parse_args()
    main(ARGS.annotation_file, ARGS.tf_list_file,
         ARGS.reveng_network_files, ARGS.out_file)
