import argparse
import pandas as pd
import networkx as nx
from data_utils import load_annotation, load_reveng_network, load_gsnetwork, map_probes

def shortest_path(net_graph, src, tgt):
    spath = None
    if src in net_graph and tgt in net_graph:
        try:
            spath = (nx.shortest_path(net_graph, src, tgt), src, tgt)
        except nx.NetworkXNoPath:
            spath = (None, src, tgt)
    else:
        if src in net_graph:
            spath = (None, src, None)
        elif tgt in net_graph:
            spath = (None, None, tgt)
        else:
            spath = (None, None, None)
    return spath

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

def eval_network(annot_file, net_file, gs_file, max_dist, max_edges):
    annot_df = load_annotation(annot_file)
    gs_net = map_probes(load_gsnetwork(gs_file), annot_df)
    gs_nedges = gs_net.shape[0]
    gs_tf_nodes = set(gs_net.TFPROBE)
    gs_tgt_nodes = set(gs_net.TARGETPROBE)
    gs_pos_edges = set([(x, y) for x in gs_tf_nodes for y in gs_tgt_nodes if x != y])
    gs_tru_edges = set([(x, y) for x, y in zip(gs_net.TFPROBE, gs_net.TARGETPROBE)])
    gs_fls_edges = list(gs_pos_edges - gs_tru_edges)
    gs_nodes = gs_tf_nodes | gs_tgt_nodes
    rv_net = select_edges(load_reveng_network(net_file), max_edges = max_edges)
    #rv_net_nodes = set(rv_net.source) | set(rv_net.target)
    rv_net_graph = nx.from_pandas_edgelist(rv_net, edge_attr='wt')
    gs_common_nodes = sum((1 if x in rv_net_graph else 0 for x in gs_nodes))
    gs_common_edges = sum((1 if (x in rv_net_graph and y in rv_net_graph) else 0
                           for x, y in zip(gs_net.TFPROBE, gs_net.TARGETPROBE)))
    gs_spath = [shortest_path(rv_net_graph, x, y)
                for x, y in zip(gs_net.TFPROBE, gs_net.TARGETPROBE)]
    #gs_spath.sort()
    dist_histogram = [0 for x in range(max_dist+1)]
    spath_graph_nodes = set(x for _, x, _ in gs_spath if x)
    spath_graph_nodes |= set(y for _, _, y in gs_spath if y)
    for spx, _, _ in gs_spath:
        if spx is not None and (len(spx)-1) <= max_dist:
            dist_histogram[len(spx)-1] += 1
            spath_graph_nodes.update(spx)
    spath_graph = nx.Graph()
    spath_graph.add_nodes_from(spath_graph_nodes)
    for spx, _, _ in gs_spath:
        if spx is not None and (len(spx)-1) <= max_dist:
            if len(spx) > 1:
                for s_node, t_node in zip(spx, spx[1:]):
                    spath_graph.add_edge(s_node, t_node)
    fp_cts = 0
    tp_cts = 0
    for x, y in gs_fls_edges:
        if x in rv_net_graph and y in rv_net_graph:
            try:
                if nx.shortest_path_length(rv_net_graph, x, y) == 1:
                    fp_cts += 1
            except nx.NetworkXNoPath:
                pass
    for x, y in gs_tru_edges:
        if x in rv_net_graph and y in rv_net_graph:
            try:
                if nx.shortest_path_length(rv_net_graph, x, y) == 1:
                    tp_cts += 1
            except nx.NetworkXNoPath:
                pass
    hist_data = {
        'DIST'  : [x for x in range(max_dist+1)],
        'EDGN'  : dist_histogram,
        'FPCTS' : [float(tp_cts)/(tp_cts + fp_cts) for _ in range(max_dist+1)],
        'PCTGS' : [float(x)*100/gs_nedges for x in dist_histogram],
        'PCTCM' : [float(x)*100/gs_common_edges for x in dist_histogram],
        'PCTSP' : [float(x)*100/spath_graph.number_of_edges()
                   for x in dist_histogram],
        'GRSP'  : [(spath_graph.number_of_nodes(), spath_graph.number_of_edges())
                   for _ in range(max_dist+1)],
        'GRGS'  : [(len(gs_nodes), gs_nedges) for _ in range(max_dist+1)],
        'GRCM'  : [(str(gs_common_nodes), str(gs_common_edges))
                   for _ in range(max_dist+1)]
    }
    dist_histogram_df = pd.DataFrame(data=hist_data)
    #print(net_file)
    #print(dist_histogram_df)
    return hist_data

def compare_eval_network(annot_file, net_files, gs_file, max_dist, max_edges):
    nhdat = [eval_network(annot_file, fx, gs_file, max_dist, max_edges)
             for fx in net_files]
    gs_cmp_data =   {str(net_files[x]) : 
                     [nhdat[x]['GRGS'][0][0], nhdat[x]['GRGS'][0][1]] + 
                     [nhdat[x]['GRCM'][0][0], nhdat[x]['GRCM'][0][1]] + 
                     [nhdat[x]['GRSP'][0][0], nhdat[x]['GRSP'][0][1]] + 
                     [nhdat[x]['EDGN'][y] for y in range(max_dist+1)] +
                     [nhdat[x]['FPCTS'][0]]
                     for x in range(len(net_files))}
    gs_cmp_data_df = pd.DataFrame(data=gs_cmp_data)
    print(gs_cmp_data_df.to_csv(sep='\t',index=False))

if __name__ == "__main__":
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("annotation_file",
                        help="""annotation file
                                (a tab seperated file mapping probe to ids)""")
    PARSER.add_argument("gs_network_file",
                        help="""gold standard network
                                (tab seperated file of TF-TARGET interactions)""")
    PARSER.add_argument("reveng_network_files", nargs="+",
                        help="""network build from a reverse engineering methods
                                (currenlty supported: eda, adj, tsv)""")
    PARSER.add_argument("-d", "--dist", type=int, default=3,
                        help="max. number of hops allowed")
    PARSER.add_argument("-x", "--max_edges", type=int,
                        help="""comma seperated names of the network;
                                should have as many names as the number of networks""")
    ARGS = PARSER.parse_args()
    compare_eval_network(ARGS.annotation_file, ARGS.reveng_network_files,
                         ARGS.gs_network_file, ARGS.dist, ARGS.max_edges)
