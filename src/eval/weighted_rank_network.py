from typing import List
import argparse
import pandas as pd
import numpy as np
import networkx as nx
import data_utils as du

def compute_tps(annot_file, gs_file, net_file,
                wt_attr, max_edges, reverse_order):
    annot_df = du.load_annotation(annot_file)
    gs_net = du.load_gsnetwork(gs_file)
    gs_net = gs_net.loc[:, ['TF', 'TARGET']]
    # print("Network size (Before Probe Mapping)", gs_net.shape)
    gs_net = du.map_probes(gs_net, annot_df)
    gs_tru_edges = set((x, y) for x, y in zip(gs_net.loc[:, "TFPROBE"],
                                              gs_net.loc[:, "TARGETPROBE"]))
    #
    #
    rv_net = select_edges(du.load_reveng_network(net_file, wt_attr_name=wt_attr),
                          wt_attr_name=wt_attr,
                          max_edges=max_edges,
                          reverse_order=reverse_order)
    rv_net_graph = nx.from_pandas_edgelist(rv_net, edge_attr=wt_attr)
    tp_edges = set([])
    for udx, vdx in gs_tru_edges:
        if udx in rv_net_graph and vdx in rv_net_graph:
            try:
                if nx.shortest_path_length(rv_net_graph, udx, vdx) == 1:
                    tp_edges.add((udx, vdx) if udx < vdx else (vdx, udx))
            except nx.NetworkXNoPath:
                pass
    return tp_edges

def find_tp_weights(annot_file, gs_file, network_files,
                    wt_attr, max_edges, reverse_order):
    all_tps = set([])
    all_net_tps = [compute_tps(annot_file, gs_file, net_file, wt_attr,
                    max_edges, reverse_order) for net_file in network_files]
    for net_tps in all_net_tps:
        all_tps = all_tps | net_tps
    print("""
        Network TPS : """,  [len(net_tps) for net_tps in all_net_tps],
                            len(all_tps))
    return [float(len(net_tps))/float(len(all_tps)) for net_tps in all_net_tps]


def select_smallest_edges(net_df: pd.DataFrame, wt_attr_name: str = 'wt',
                          max_edges: int = None):
    if max_edges is None or max_edges >= net_df.shape[0]:
        return net_df
    cur_cols = net_df.columns
    net_df = net_df.nsmallest(n=max_edges, columns=wt_attr_name, keep='all')
    #print(net_df.iloc[0, :])
    return net_df.loc[:, cur_cols]

def select_edges(net_df: pd.DataFrame, wt_attr_name: str = 'wt',
                 max_edges: int = None, reverse_order: bool = False):
    if max_edges is None or max_edges >= net_df.shape[0]:
        return net_df
    cur_cols = net_df.columns
    maxwt_attr_name = wt_attr_name + '_max'
    net_df[maxwt_attr_name] = net_df[wt_attr_name].abs()
    if reverse_order is True:
        net_df = net_df.nlargest(n=max_edges, columns=maxwt_attr_name)
    else:
        net_df = net_df.nlargest(n=max_edges, columns=maxwt_attr_name)
    #print(net_df.iloc[0, :])
    return net_df.loc[:, cur_cols]


def count_nonnan(row_x):
    return np.count_nonzero(~np.isnan(row_x))


def weighted_rank_network(network_files, network_names=None, network_weights=None,
                          wt_attr: str = 'wt', max_edges: int = None,
                          max_out_edges: int = None, max_avg: int = None,
                          reverse_order: bool = False):
    if network_weights is None:
        network_weights = [1 for _ in range(len(network_files))]
    if network_names is None:
        network_names = ['net_' + str(ix) for ix in range(len(network_files))]
    cmb_network = pd.DataFrame(columns=['source', 'target'])
    for nx_name, nx_file in zip(network_names, network_files):
        ndf = select_edges(du.load_reveng_network(nx_file, wt_attr),
                           wt_attr, max_edges, reverse_order)
        ndf = ndf.rename(columns={wt_attr: nx_name})
        if reverse_order is True:
            ndf[nx_name] = ndf[nx_name].rank(ascending=True)
        else:
            ndf[nx_name] = ndf[nx_name].rank(ascending=False)
        cmb_network = cmb_network.merge(ndf, how='outer', on=['source', 'target'])
        #print(str(nx_name), nx_file, ndf.shape, cmb_network.columns, cmb_network.shape)
    for nx_wt, nx_name in zip(network_weights, network_names):
        mx_rank = np.nanmax(cmb_network[nx_name])
        cmb_network.loc[cmb_network[nx_name].isna(), nx_name] = mx_rank+1
        cmb_network[nx_name] = nx_wt * cmb_network[nx_name]
    cmb_network['RANKAVG'] = cmb_network[network_names].sum(axis=1)
    cmb_network.loc[:, 'RANKAVG'] = cmb_network['RANKAVG'].div(float(sum(network_weights)))
    if max_avg:
        cmb_network = cmb_network.loc[cmb_network.RANKAVG <= max_avg, : ]
    else:
        cmb_network = select_smallest_edges(cmb_network, 'RANKAVG', max_out_edges)
    cmb_network.sort_values('RANKAVG', inplace=True)
    return cmb_network

def main(network_files: List[str], network_names: str,
         network_weights: str, wt_attr: str, max_edges: int,
         max_out_edges: int, max_avg: int, out_file: str,
         reverse_order: bool, annotation_file: str, gs_network_file: str) -> bool:
    if network_names:
        network_names = network_names.split(",")
        if len(network_names) != len(network_files):
            print("Length of network names should be equal to length of network files")
            return False
    if network_weights is not None:
        network_weights = [float(x) for x in network_weights.split(",")]
        if len(network_weights) != len(network_files):
            print("Length of network weights should be equal to length of network files")
            return False
    if (network_weights is None) and (gs_network_file is not None):
        if annotation_file is None:
            print("Given GS Network, need annotation.")
            return False
        network_weights = find_tp_weights(annotation_file, gs_network_file, network_files,
                                            wt_attr, max_edges, reverse_order)
    print("""
       PARSED : network_names   : %s
       PARSED/COMPUTED : network_weights : %s """ % (str(network_names),
                                                     str(network_weights)))
    combine_df = weighted_rank_network(network_files, network_names,
                                       network_weights, wt_attr, max_edges,
                                       max_out_edges, max_avg, reverse_order)
    combine_df.to_csv(out_file, sep='\t', index=False)
    return True


if __name__ == "__main__":
    PROG_DESC = """
    Compute a weighted rank integrated of the input networks.
    Network integration is computed as the weighted average of the ranks
    of edges of the input networks.
    Outputs a tab-seperated values with weights corresponding to each network
    in a separate column, and maximum and average weight.
    """
    PARSER = argparse.ArgumentParser(description=PROG_DESC)
    PARSER.add_argument("network_files", nargs="+",
                        help="""network build from a reverse engineering methods
                                (currenlty supported: eda, adj, tsv)""")
    PARSER.add_argument("-n", "--network_names", type=str,
                        help="""comma seperated names of the network;
                                should have as many names as the number of networks""")
    PARSER.add_argument("-w", "--network_weights", type=str,
                        help="""weights for  input networks;
                                should have as many weights as the number of networks""",
                        default=None)
    PARSER.add_argument("-g", "--max_avg", type=float, # default=1500000.0,
                        help="""maximum average rank allowed """)
    PARSER.add_argument("-x", "--max_edges", type=int,
                        help="""maximum edges allowed in the input network""")
    PARSER.add_argument("-y", "--max_out_edges", type=int,
                        help="""maximum edges allowed in the output network""")
    PARSER.add_argument("-t", "--wt_attr", type=str, default='wt',
                        help="name of weight attribute")
    PARSER.add_argument("-r", "--reverse_order", action='store_true',
                        help="""Order the edges ascending order""")
    PARSER.add_argument("-a", "--annotation_file", type=str, default=None,
                        help="Annotation File")
    PARSER.add_argument("-s", "--gs_network_file", type=str, default=None,
                        help="Gold standard network file")
    PARSER.add_argument("-o", "--out_file",
                        type=argparse.FileType(mode='w'), required=True,
                        help="output file in tab-seperated format")
    ARGS = PARSER.parse_args()
    print("""
       ARG : network_files : %s
       ARG : network_names : %s
       ARG : network_weights : %s
       ARG : max_edges : %s
       ARG : max_out_edges : %s
       ARG : max_avg : %s
       ARG : wt_attr : %s
       ARG : reverse_order : %s
       ARG : annotation_file : %s
       ARG : gs_network_file : %s
       ARG : out_file : %s """ %
          (str(ARGS.network_files), str(ARGS.network_names), str(ARGS.network_weights),
           str(ARGS.max_edges), str(ARGS.max_out_edges), str(ARGS.max_avg),
           str(ARGS.wt_attr), str(ARGS.reverse_order),
           str(ARGS.annotation_file), str(ARGS.gs_network_file), 
           str(ARGS.out_file)))
    if not main(ARGS.network_files, ARGS.network_names, ARGS.network_weights,
                ARGS.wt_attr, ARGS.max_edges, ARGS.max_out_edges,
                ARGS.max_avg, ARGS.out_file, ARGS.reverse_order,
                ARGS.annotation_file, ARGS.gs_network_file):
        PARSER.print_usage()
