import argparse
import pandas as pd
import networkx as nx
import data_utils as du

def all_pair_edges(mod_genes):
    return pd.DataFrame(
            [(x, y) for x in mod_genes for y in mod_genes if x < y],
            columns=['SRC', 'TGT'])

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

def shorest_path_graph(mod_genes, annot_map, rv_net_graph):
    gs_gene_ids = {x: y for x, y in zip(mod_genes.PID_PROBE, mod_genes.PID_ID)  if x in rv_net_graph}
    gs_probes = set(list(mod_genes.PID_PROBE))
    rev_subgraph = rv_net_graph.subgraph(gs_probes)
    print('SUBNET ', len(rev_subgraph), nx.number_connected_components(rev_subgraph))
    spath_graph = nx.Graph()
    for probe in rev_subgraph.nodes():
        spath_graph.add_node(probe, ID=gs_gene_ids[probe],
                             COLOR='GREEN')
    for s_node, t_node in rev_subgraph.edges():
        spath_graph.add_edge(s_node, t_node, COLOR='GREEN')
    if nx.number_connected_components(rev_subgraph) == 1:
        return spath_graph
    gs_net = all_pair_edges(gs_gene_ids.keys())
    gs_spath = [shortest_path(rv_net_graph, x, y)
                for x, y in zip(gs_net.loc[:, 'SRC'], gs_net.loc[:, 'TGT'])]
    sorted(gs_spath, reverse=False, key=lambda x: len(x[0]) if x[0] else 0)
    spath_graph_nodes = set(x for _, x, _ in gs_spath if x)
    spath_graph_nodes |= set(y for _, _, y in gs_spath if y)
    ncsx = nx.number_connected_components(spath_graph)
    for spx, src, tgt in gs_spath:
        #spx = gs_spath[(x, y)][0]
        try:
            nlen = nx.shortest_path_length(spath_graph, src, tgt)
            if nlen >= 0:
                continue
        except nx.NetworkXNoPath:
            pass
        if spx is not None:
            #spath_graph_nodes.update(spx)
            if len(spx) > 1:
                for s_node, t_node in zip(spx, spx[1:]):
                    if s_node not in spath_graph:
                        spath_graph.add_node(s_node, ID=annot_map[s_node], COLOR='RED')
                    if t_node not in spath_graph:
                        spath_graph.add_node(t_node, ID=annot_map[t_node], COLOR='RED')
                    spath_graph.add_edge(s_node, t_node, 
                                         COLOR='YELLOW' if s_node not in gs_probes and s_node not in gs_probes else 'RED')
                    ncsx = nx.number_connected_components(spath_graph)
                    if ncsx == 1:
                        break
        if ncsx == 1:
            break
    return spath_graph

def main(annot_file, pathway_nodes_file, net_file,
         out_file, wt_attr, max_edges, reverse_order):
    annot_df = du.load_annotation(annot_file)
    pway_df = pd.read_csv(pathway_nodes_file, delimiter="\t", names=['PID'])
    print(pway_df.columns, pway_df.shape)
    pway_df = du.map_probes_cols(pway_df, annot_df, ['PID'],
                                right_onc='ID',
                                probe_suffix='_PROBE', id_suffix='_ID')
    print(pway_df.columns, pway_df.shape)
    annot_map = {x:y for x, y in zip(annot_df.PROBE, annot_df.ID)}
    rv_net = du.select_edges(du.load_reveng_network(net_file, wt_attr_name=wt_attr),
                             wt_attr_name=wt_attr,
                             max_edges=max_edges,
                             reverse_order=reverse_order)
    rv_net_graph = nx.from_pandas_edgelist(rv_net, edge_attr=wt_attr)
    spath_graph = shorest_path_graph(pway_df, annot_map, rv_net_graph)
    print("SPATH SIZE", len(spath_graph), nx.number_connected_components(spath_graph))
    nx.write_gml(spath_graph, out_file)

if __name__ == "__main__":
    PROG_DESC = """ Add edges to form a connected component """
    PARSER = argparse.ArgumentParser(description=PROG_DESC)
    PARSER.add_argument("annotation_file",
                        help="""annotation file
                                (a tab seperated file mapping probe to ids)""")
    PARSER.add_argument("pathway_nodes_file",
                        help="""Nodes of interest""")
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
                        type=str, required=True,
                        help="output file in tab-seperated format")
    ARGS = PARSER.parse_args()
    main(ARGS.annotation_file, ARGS.pathway_nodes_file, ARGS.reveng_network_file,
         ARGS.out_file, ARGS.wt_attr, ARGS.max_edges, ARGS.reverse_order)
