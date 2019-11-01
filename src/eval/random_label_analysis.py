import argparse
import numpy as np
import pandas as pd
import networkx as nx
import multiprocessing as mp
from typing import List
import data_utils as du
from mpi4py import MPI

def select_edges(net_df: pd.DataFrame, wt_attr_name: str = 'wt',
                 max_edges: int = None):
    if max_edges is None or max_edges >= net_df.shape[0]:
        return net_df
    cur_cols = net_df.columns
    maxwt_attr_name = wt_attr_name + '_max'
    net_df[maxwt_attr_name] = net_df[wt_attr_name].abs()
    net_df = net_df.nlargest(n=max_edges, columns=maxwt_attr_name)
    #print(net_df.iloc[0, :])
    return net_df.loc[:, cur_cols]

def rand_label_distribution(annot_file: str, gs_file: str, network_file: str,
                            wt_attr_name: str = 'wt',
                            max_edges :int = None,
                            nshuffles :int = 5) -> List[int]:
    annot_df = du.load_annotation(annot_file)
    gs_net = du.map_probes(du.load_gsnetwork(gs_file), annot_df)
    net_df: pd.DataFrame = select_edges(du.load_reveng_network(network_file), 
                                        wt_attr_name, max_edges)
    rev_net: nx.Graph = nx.from_pandas_edgelist(net_df, edge_attr=wt_attr_name)
    n_nodes = nx.number_of_nodes(rev_net)
    node_names_map = {y:x for x, y in zip(range(n_nodes), nx.nodes(rev_net))}
    #print(node_names_list[0], node_names_map['249052_at'])
    rev_net_nx = nx.convert_node_labels_to_integers(rev_net)
    gs_net_nx = [(node_names_map[x], node_names_map[y])
                  for x, y in zip(gs_net.TFPROBE, gs_net.TARGETPROBE)
                  if x in node_names_map and y in node_names_map]
    results = np.zeros(nshuffles)
    for idx in range(nshuffles):
        narray = np.arange(n_nodes)
        np.random.shuffle(narray)
        relabel_map = dict(zip(range(n_nodes),narray))
        rev_net_nx_rlbl = nx.relabel_nodes(rev_net_nx, mapping=relabel_map)
        results[idx] = sum(1 if rev_net_nx_rlbl.has_edge(x, y) else 0 for x,y in gs_net_nx)
    return results


def main(annot_file: str, gs_file: str, network_file: str, max_edges:int ) -> None:
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    prop_data = rand_label_distribution(annot_file, gs_file, network_file,
                                         'wt', max_edges)
    print(prop_data)
    #prop_df = pd.DataFrame(data=prop_data, index=cnames)
    #print(prop_df.transpose().to_csv(sep="\t", index=True))


if __name__ == "__main__":
    PROG_DESC = """Random relabelling anlysis"""
    ARGPARSER = argparse.ArgumentParser(description=PROG_DESC)
    ARGPARSER.add_argument("-x", "--max_edges", type=int,
                        help="""maximum number of edges""")
    ARGPARSER.add_argument("-o", "--out_file",
                           type=str,
                           help="output file in text format")
    ARGPARSER.add_argument("annotation_file",
                           help="""annotation file
                                   (a tab seperated file mapping probe to ids)""")
    ARGPARSER.add_argument("gs_network_file",
                           help="""gold standard network
                                   (tab seperated file of TF-TARGET interactions)""")
    ARGPARSER.add_argument("network_file")
    CMDARGS = ARGPARSER.parse_args()
    main(CMDARGS.annotation_file, CMDARGS.gs_network_file, CMDARGS.network_file, CMDARGS.max_edges)
