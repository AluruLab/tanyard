import argparse
import networkx as nx
import pandas as pd
import data_utils as du
import itertools

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


def find_connected_subnetwork(annotation_file, network_file,
                              tf_list_file, output_file):
    #
    tf_lst_df = pd.read_csv(tf_list_file)
    annot_df = du.load_annotation(annotation_file)
    rev_net = du.map_probes_cols(du.load_reveng_network(network_file, wt_attr_name='RANKAVG'),
                                 annot_df=annot_df,
                                 col_names=['source', 'target'], probe_suffix='_PROBE',
                                 id_suffix='_ID')
    rv_net_graph = nx.from_pandas_edgelist(rev_net,
                                           source='source_ID',
                                           target='target_ID')
    print(len(set(rev_net.source_ID) & set(tf_lst_df.TF)))
    print(len(set(rev_net.target_ID) & set(tf_lst_df.TF)))
    gs_spath = [shortest_path(rv_net_graph, x, y)
                for x, y in itertools.combinations(tf_lst_df.TF, 2)]
    spath_graph_nodes = set(x for _, x, _ in gs_spath if x)
    spath_graph_nodes |= set(y for _, _, y in gs_spath if y)
    print(rv_net_graph.size(), tf_lst_df.size, len(spath_graph_nodes))
    for spx, _, _ in gs_spath:
        if spx is not None:
            spath_graph_nodes.update(spx)
    spath_graph = nx.Graph()
    spath_graph.add_nodes_from(spath_graph_nodes)
    for spx, _, _ in gs_spath:
        if spx is not None and len(spx) > 1:
            for s_node, t_node in zip(spx[:-1], spx[1:]):
                spath_graph.add_edge(s_node, t_node)
    spath_df =  nx.to_pandas_edgelist(spath_graph)
    spath_df.to_csv(output_file, sep="\t", index=False)


def find_connected_subnetwork_probes(annotation_file, network_file,
                              tf_list_file, output_file):
    #
    annot_df = du.load_annotation(annotation_file)
    tf_lst_df = du.map_id2probe_cols(pd.read_csv(tf_list_file), 
                                 annot_df=annot_df,
                                 col_names=['TF'], probe_suffix='_PROBE',
                                 id_suffix='_ID')
    rev_net = du.load_reveng_network(network_file, wt_attr_name='RANKAVG')
    rv_net_graph = nx.from_pandas_edgelist(rev_net,
                                           source='source',
                                           target='target')
    print(len(set(rev_net.source) & set(tf_lst_df.TF_PROBE)))
    print(len(set(rev_net.target) & set(tf_lst_df.TF_PROBE)))
    gs_spath = [shortest_path(rv_net_graph, x, y)
                for x, y in itertools.combinations(tf_lst_df.TF_PROBE, 2)]
    spath_graph_nodes = set(x for _, x, _ in gs_spath if x)
    spath_graph_nodes |= set(y for _, _, y in gs_spath if y)
    print(rv_net_graph.size(), tf_lst_df.size, len(spath_graph_nodes))
    for spx, _, _ in gs_spath:
        if spx is not None:
            spath_graph_nodes.update(spx)
    spath_graph = nx.Graph()
    spath_graph.add_nodes_from(spath_graph_nodes)
    for spx, _, _ in gs_spath:
        if spx is not None and len(spx) > 1:
            for s_node, t_node in zip(spx[:-1], spx[1:]):
                spath_graph.add_edge(s_node, t_node)
    spath_df =  nx.to_pandas_edgelist(spath_graph)
    spath_df.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    PROG_DESC = """
    Add a genes to complete a subnetwork of a list of input genes
    """
    PARSER = argparse.ArgumentParser(description=PROG_DESC)
    PARSER.add_argument("annotation_file",
                        help="""annotation file
                                (a tab seperated file mapping probe to ids)""")
    PARSER.add_argument("tf_list_file",
                        help="""List of transcription factors in a file""")
    PARSER.add_argument("network_file", 
                        help="""network build from a reverse engineering methods
                                (currenlty supported: eda, adj, tsv)""")
    PARSER.add_argument("output_file", 
                        help="""sub nnetwork constructed from the complete networks""")
    PARSER.add_argument("-t", "--type", choices=['id', 'probe'], default='id',
                        help="""use id or probe for network node names""")
    ARGS = PARSER.parse_args()
    if ARGS.type == 'id':
        find_connected_subnetwork(ARGS.annotation_file, ARGS.network_file,
                                  ARGS.tf_list_file, ARGS.output_file)
    else:
        find_connected_subnetwork_probes(ARGS.annotation_file, ARGS.network_file,
                                         ARGS.tf_list_file, ARGS.output_file)
