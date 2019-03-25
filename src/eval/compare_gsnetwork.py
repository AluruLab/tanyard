import pandas as pd
import networkx as nx
import argparse
from data_utils import load_annotation, load_reveng_network, load_gsnetwork, map_probes


def eval_network(annot_file, net_file, gs_file, max_dist):
    annot_df = load_annotation(annot_file)
    rv_net = load_reveng_network(net_file)
    gs_net = map_probes(load_gsnetwork(gs_file), annot_df)
    rv_net_graph = nx.from_pandas_edgelist(rv_net, edge_attr='wt')
    gs_nedges = gs_net.shape[0]
    gs_nodes = set(gs_net.TFPROBE) | set(gs_net.TARGETPROBE)
    gs_common_nodes = sum((1 if x in rv_net_graph else 0 for x in gs_nodes))
    gs_common_edges = sum((1 if (x in rv_net_graph and y in rv_net_graph) else 0 
                           for x,y in zip(gs_net.TFPROBE, gs_net.TARGETPROBE)))
    gs_spath = [(nx.shortest_path_length(rv_net_graph, x, y), x, y) 
                    if x in rv_net_graph and y in rv_net_graph else (float("inf"), x, y) 
                    for x,y in zip(gs_net.TFPROBE, gs_net.TARGETPROBE) ]
    #gs_spath.sort()
    dist_histogram = [0 for x in range(max_dist+1)]
    for x, _, _ in gs_spath:
        if x <= max_dist:
            dist_histogram[x] += 1
    dist_histogram_df = pd.DataFrame(data={'DIST'   : [x for x in range(max_dist+1)],
                                           'CNT.'   : dist_histogram,
                                           'DPCTX'  : [float(x)*100/gs_nedges for x in dist_histogram],
                                           'DPCTY'  : [float(x)*100/gs_common_edges for x in dist_histogram],
                                           'GSNDS'  : [len(gs_nodes) for _ in range(max_dist+1)],
                                           'GSEDG'  : [gs_nedges for _ in range(max_dist+1)],
                                           'GSRVN'  : [gs_common_nodes for _ in range(max_dist+1)],
                                           'GSRVE'  : [gs_common_edges for _ in range(max_dist+1)]})
    print(dist_histogram_df)
    return dist_histogram_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("annotation_file",
        help="annotation file (a tab seperated file) mapping probe to ids")
    parser.add_argument("gs_network_file",
        help="gold standard network (tab seperated file of TF-TARGET interactions)")
    parser.add_argument("reveng_network_file",
        help="network build from a reverse engineering methods (currenlty supported: eda)")
    parser.add_argument("-d", "--dist", type=int, default=3,
                        help="max. number of hops allowed")
    args = parser.parse_args()
    eval_network(args.annotation_file, args.reveng_network_file, args.gs_network_file, args.dist)
