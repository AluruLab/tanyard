import pandas as pd
import numpy as np
from typing import List
from data_utils import load_reveng_network
import argparse


def select_smallest_edges(net_df: pd.DataFrame, wt_attr_name: str = 'wt',
        max_edges: int = None):
    if max_edges is None or max_edges >= net_df.shape[0]:
        return net_df
    cur_cols = net_df.columns
    net_df = net_df.nsmallest(n=max_edges, columns=wt_attr_name, keep='all')
    #print(net_df.iloc[0, :])
    return net_df.loc[:, cur_cols]

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


def count_nonnan(row_x):
    return np.count_nonzero(~np.isnan(row_x))


def weighted_rank_network(network_files, network_names=None, network_weights=None,
                          wt_attr:str = 'wt', max_edges: int = None,
                          max_out_edges: int = None, max_avg: int = None):
    if network_weights is None:
        network_weights = [1 for _ in range(len(network_files))]
    if network_names is None:
        network_names = ['net_' + str(ix) for ix in range(len(network_files))]
    cmb_network = pd.DataFrame(columns=['source', 'target'])
    for nx_name, nx_file in zip(network_names, network_files):
        ndf = select_edges(load_reveng_network(nx_file, wt_attr),
                           wt_attr, max_edges)
        ndf = ndf.rename(columns={wt_attr: nx_name})
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
         max_out_edges: int, max_avg: int, out_file: str) -> bool:
    if network_names:
        network_names = network_names.split(",")
        if len(network_names) != len(network_files):
            print("Length of network names should be equal to length of network files")
            return False
    if network_weights:
        network_weights = [float(x) for x in network_weights.split(",")]
        if len(network_weights) != len(network_files):
            print("Length of network weights should be equal to length of network files")
            return False
    print("""
       PARSED : network_names   : %s
       PARSED : network_weights : %s """ % (str(network_names), str(network_weights)))
    combine_df = weighted_rank_network(network_files, network_names, network_weights,
                                       wt_attr, max_edges, max_out_edges, max_avg)
    combine_df.to_csv(out_file, sep='\t', index=False)
    return True


if __name__ == "__main__":
    PROG_DESC = """
    Finds a union of input networks.
    Network union is computed as the union of edges of the input networks.
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
                        help="""comma seperated names of the network;
                                should have as many names as the number of networks""")
    PARSER.add_argument("-g", "--max_avg", type=float, # default=1500000.0,
                        help="""maximum average rank allowed """)
    PARSER.add_argument("-x", "--max_edges", type=int,
                        help="""maximum edges allowed in the input network""")
    PARSER.add_argument("-y", "--max_out_edges", type=int,
                        help="""maximum edges allowed in the output network""")
    PARSER.add_argument("-t", "--wt_attr", type=str, default='wt',
                        help="name of weight attribute")
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
       ARG : out_file : %s """ %
          (str(ARGS.network_files), str(ARGS.network_names), str(ARGS.network_weights),
           str(ARGS.max_edges), str(ARGS.max_out_edges), str(ARGS.max_avg),
           str(ARGS.wt_attr), str(ARGS.out_file)))
    if not main(ARGS.network_files, ARGS.network_names, ARGS.network_weights,
                ARGS.wt_attr, ARGS.max_edges, ARGS.max_out_edges,
                ARGS.max_avg, ARGS.out_file):
        PARSER.print_usage()
