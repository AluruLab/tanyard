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

def main(annot_file: str, net_file: str, method: str, out_file: str):
    annot_df = du.load_annotation(annot_file)
    # tflst_df = du.map_atid2probes(pd.read_csv(tf_list_file, sep= r'\s+'),
    #                               annot_df)
    rv_net = du.load_reveng_network(net_file)
    rv_net_node_lst = list(set(rv_net.source) | set(rv_net.target))
    print("No. of Nodes : ", len(rv_net_node_lst))
    rv_net_node_map = {y:x for x, y in enumerate(rv_net_node_lst)}
    edge_rcds = rv_net.loc[:, ['source', 'target', 'wt']].to_records(index=False)
    rv_net_edge_map = {(rv_net_node_map[x], rv_net_node_map[y]):(1/float(w)) for x, y, w in edge_rcds}
    print("No. of Edges : ", len(rv_net_edge_map))
#     rv_net_graph: nx.Graph = nx.from_pandas_edgelist(rv_net, edge_attr='wt')
#     comm_dict = partition(rv_net_graph)
#     for comm in set(comm_dict.values()):
#         print("Community %d"%comm)
    rv_net_edges = [*rv_net_edge_map]
    rv_net_wts = [rv_net_edge_map[x] for x in rv_net_edges]
    rv_net_igraph = igx.Graph(rv_net_edges)
    node_clst = None
    if method == METHOD_FGREEDY:
        net_mods = rv_net_igraph.community_fastgreedy(weights=rv_net_wts)
        node_clst = net_mods.as_clustering()
    elif method == METHOD_LABELPR:
        node_clst = rv_net_igraph.community_label_propagation(weights=rv_net_wts)
    elif method == METHOD_LEIGENV:
        node_clst = rv_net_igraph.community_leading_eigenvector(weights=rv_net_wts)
    elif method == METHOD_SPINGLS:
        node_clst = rv_net_igraph.community_spinglass(weights=rv_net_wts, spins=100)
    elif method == METHOD_MULTILV:
        node_clst = rv_net_igraph.community_multilevel(weights=rv_net_wts)
    elif method == METHOD_INFOMAP:
        node_clst = rv_net_igraph.community_infomap(edge_weights=rv_net_wts, trials=100)
    elif method == METHOD_WALKTRP:
        net_mods = rv_net_igraph.community_walktrap(weights=rv_net_wts)
        node_clst = net_mods.as_clustering()
    else:
        print("method not supported")
    rdf = pd.DataFrame({"GENE": rv_net_node_lst,
                        "CLUST_ID": node_clst.membership})
    out_df = du.map_probes_cols(net_df=rdf, annot_df=annot_df,
                                col_names=['GENE'], probe_suffix='_PROBE',
                                id_suffix='_ID')
    out_df.loc[:, ['GENE_ID', 'CLUST_ID']].to_csv(out_file, sep="\t")


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
         ARGS.method, ARGS.out_file)
