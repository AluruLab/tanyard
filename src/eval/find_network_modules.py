import argparse
import pandas as pd
#import networkx as nx
import igraph as igx
import data_utils as du
#from data_utils import load_annotation, map_probes_cols, load_reveng_network, map_probes_cols
# from modularity_maximization import partition
# from modularity_maximization.utils import get_modularity

METHOD_FGREEDY = 'fast'
METHOD_LABELPR = 'label'
METHOD_LEIGENV = 'eigen'
METHOD_SPINGLS = 'spin'
METHOD_MULTILV = 'multi'
METHOD_INFOMAP = 'info'
METHOD_WALKTRP = 'walk'

METHOD_DESCRIPTION = {
        METHOD_FGREEDY : 'Fast Greedy',
        METHOD_LABELPR : 'Label Propagation',
        METHOD_LEIGENV : 'Largest Eigen Value',
        METHOD_SPINGLS : 'Spinglass',
        METHOD_MULTILV : 'Multi-level',
        METHOD_INFOMAP : 'Infomap',
        METHOD_WALKTRP : 'Random walk'
}

METHOD_KEYS = list(METHOD_DESCRIPTION.keys())

MODULE_RUNS = 0

def run_community_detection(method, in_graph, nspins=50, ntrials=100):
    global MODULE_RUNS
    node_clst = None
    if method == METHOD_FGREEDY:
        net_mods = in_graph.community_fastgreedy(weights='wt')
        node_clst = net_mods.as_clustering()
    elif method == METHOD_LABELPR:
        node_clst = in_graph.community_label_propagation(weights='wt')
    elif method == METHOD_LEIGENV:
        node_clst = in_graph.community_leading_eigenvector(weights='wt')
    elif method == METHOD_SPINGLS:
        node_clst = in_graph.community_spinglass(weights='wt', spins=nspins)
    elif method == METHOD_MULTILV:
        node_clst = in_graph.community_multilevel(weights='wt')
    elif method == METHOD_INFOMAP:
        node_clst = in_graph.community_infomap(edge_weights='wt', trials=ntrials)
    elif method == METHOD_WALKTRP:
        net_mods = in_graph.community_walktrap(weights='wt')
        node_clst = net_mods.as_clustering()
    else:
        print("method not supported")
        return None
    MODULE_RUNS += 1
    return node_clst

def find_graph_modules_recursive(method, in_graph, clust_id, nlevel, nlimit):
    clst_id_map = {}
    node_clst = run_community_detection(method, in_graph)
    for vert_list in node_clst:
        if len(vert_list) > nlimit and nlevel < 4:
            in_sub_graph = in_graph.subgraph(vert_list)
            rvclst_id_map, cid = find_graph_modules_recursive(method,  in_sub_graph, 
                                                              clust_id, nlevel + 1, nlimit)
            clust_id = cid
            for vid in rvclst_id_map:
               clst_id_map[vert_list[vid]] = rvclst_id_map[vid]
        else:
            for vid in vert_list:
                clst_id_map[vid] = clust_id
            clust_id += 1
    return (clst_id_map, clust_id)


def find_cluster_membership_recursive(method, rv_net_node_lst, rv_net_edge_map, nlimit=1000):
    rv_net_edges = [*rv_net_edge_map]
    rv_net_wts = [rv_net_edge_map[x] for x in rv_net_edges]
    rv_net_igraph = igx.Graph(edges=rv_net_edges, edge_attrs={'wt' : rv_net_wts})
    final_clst_ids, cid = find_graph_modules_recursive(method, rv_net_igraph, 1, 1, nlimit)
    rdf = pd.DataFrame({"GENE": [rv_net_node_lst[x] for x in final_clst_ids],
                        "MODULE_ID": [final_clst_ids[x] for x in final_clst_ids]})
    return rdf


def find_cluster_membership(method, rv_net_node_lst, rv_net_edge_map, nlimit=1000):
    rv_net_edges = [*rv_net_edge_map]
    rv_net_wts = [rv_net_edge_map[x] for x in rv_net_edges]
    rv_net_igraph = igx.Graph(edges=rv_net_edges, edge_attrs={'wt' : rv_net_wts})
    node_clst = run_community_detection(method, rv_net_igraph)
    final_clst_ids = {}
    clust_id = 1
    for vert_list in node_clst:
        if len(vert_list) > nlimit:
            sub_gx = rv_net_igraph.subgraph(vert_list)
            for cgx in node_clst_sgx:
                for vid in cgx:
                    final_clst_ids[vert_list[vid]] = clust_id
                clust_id += 1
        else:
            for vid in vert_list:
                final_clst_ids[vid] = clust_id
        clust_id += 1
    rdf = pd.DataFrame({"GENE": [rv_net_node_lst[x] for x in final_clst_ids],
                        "MODULE_ID": [final_clst_ids[x] for x in final_clst_ids]})
    return rdf


def select_edges(net_df: pd.DataFrame, wt_attr_name: str = 'wt',
                 max_edges: int = None,
                 reverse_order: bool = False):
    if max_edges is None or max_edges >= net_df.shape[0]:
        return net_df
    cur_cols = net_df.columns
    maxwt_attr_name = wt_attr_name + '_max'
    net_df[maxwt_attr_name] = net_df[wt_attr_name].abs()
    if reverse_order is True:
        net_df = net_df.nsmallest(n=max_edges, columns=maxwt_attr_name)
    else:
        net_df = net_df.nlargest(n=max_edges, columns=maxwt_attr_name)
    #print(net_df.iloc[0, :])
    return net_df.loc[:, cur_cols]


def main(annot_file: str, net_file: str,
        wt_attr: str, max_edges: int, reverse_order: bool,
        method: str, out_file: str):
    annot_df = du.load_annotation(annot_file)
    # tflst_df = du.map_atid2probes(pd.read_csv(tf_list_file, sep= r'\s+'),
    #                               annot_df)
    rv_net = select_edges(du.load_reveng_network(net_file, wt_attr_name=wt_attr),
                          wt_attr_name=wt_attr,
                          max_edges=max_edges,
                          reverse_order=reverse_order)
    print("Columns : ", rv_net.columns)
    rv_net_node_lst = list(set(rv_net.source) | set(rv_net.target))
    print("No. of Nodes : ", len(rv_net_node_lst))
    rv_net_node_map = {y:x for x, y in enumerate(rv_net_node_lst)}
    edge_rcds = rv_net.loc[:, ['source', 'target', wt_attr]].to_records(index=False)
    if reverse_order is True:
       rv_net_edge_map = {(rv_net_node_map[x], rv_net_node_map[y]):
                          (1/float(w)) for x, y, w in edge_rcds}
    else:
       rv_net_edge_map = {(rv_net_node_map[x], rv_net_node_map[y]):
                          float(w) for x, y, w in edge_rcds}
    print("No. of Edges : ", len(rv_net_edge_map))
#     rv_net_graph: nx.Graph = nx.from_pandas_edgelist(rv_net, edge_attr='wt')
#     comm_dict = partition(rv_net_graph)
#     for comm in set(comm_dict.values()):
#         print("Community %d"%comm)
    rdf = find_cluster_membership_recursive(method, rv_net_node_lst, rv_net_edge_map)
    out_df = du.map_probes_cols(net_df=rdf, annot_df=annot_df,
                                col_names=['GENE'], probe_suffix='_PROBE',
                                id_suffix='_ID')
    out_df.loc[:, ['GENE_ID', 'MODULE_ID']].to_csv(out_file, sep="\t", index=False)
    global MODULE_RUNS
    print("No. of Module Runs : ", MODULE_RUNS)
    #rdf.to_csv(out_file, sep="\t")


if __name__ == "__main__":
    PROG_DESC = """Construct modules  using four methods
    as given in the review : 
    
    https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0979-8
    """
    PARSER = argparse.ArgumentParser(description=PROG_DESC)
    PARSER.add_argument("annotation_file",
                        help="""annotation file
                                (a tab seperated file mapping probe to ids)""")
    PARSER.add_argument("reveng_network_file",
                        help="""network build from a reverse engineering methods
                                (currenlty supported: eda, adj, tsv)""")
    PARSER.add_argument("-t", "--wt_attr", type=str, default='wt',
                        help="name of weight attribute")
    PARSER.add_argument("-x", "--max_edges", type=int,
                        help="""Maximum number of edges""")
    PARSER.add_argument("-r", "--reverse_order", action='store_true',
                        help="""Order the edges ascending order""")
    PARSER.add_argument("-o", "--out_file",
                        type=str,
                        help="output file in png format")
    PARSER.add_argument("-m", "--method",
                        choices=METHOD_KEYS,
                        default=METHOD_FGREEDY,
                        help="Algorithm to use for community detection. Options available are : " +
                        ';'.join([" {} - {} ".format(k, v) for k, v in METHOD_DESCRIPTION.items()]))
    ARGS = PARSER.parse_args()
    print("""Inputs:
        Annotation : {} 
        Network    : {}
        Method     : {}
        Output     : {}""".format(
            ARGS.annotation_file, ARGS.reveng_network_file,
            METHOD_DESCRIPTION[ARGS.method], ARGS.out_file))
    main(ARGS.annotation_file, ARGS.reveng_network_file,
         ARGS.wt_attr, ARGS.max_edges, ARGS.reverse_order,
         ARGS.method, ARGS.out_file)
