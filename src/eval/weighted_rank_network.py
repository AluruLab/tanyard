import pandas as pd
import numpy as np
from typing import List
from data_utils import load_reveng_network
import argparse

def count_nonnan(row_x):
    return np.count_nonzero(~np.isnan(row_x))


def consensus_network(network_files, network_names=None, network_weights=None, max_avg: int = None):
    if network_weights is None:
        network_weights = [1 for _ in range(len(network_files))]
    if network_names is None:
        network_names = ['net_' + str(ix) for ix in range(len(network_files))]
    cmb_network = pd.DataFrame(columns=['source', 'target'])
    for nx_name, nx_file in zip(network_names, network_files):
        ndf = load_reveng_network(nx_file, nx_name)
        ndf[nx_name] = ndf[nx_name].rank(ascending=False)
        cmb_network = cmb_network.merge(ndf, how='outer', on=['source', 'target'])
        #print(str(nx_name), nx_file, ndf.shape, cmb_network.columns, cmb_network.shape)
    for nx_wt, nx_name in zip(network_weights, network_names):
        mx_rank = np.nanmax(cmb_network[nx_name])
        cmb_network.loc[cmb_network[nx_name].isna(), nx_name] = mx_rank+1
        cmb_network[nx_name] = nx_wt * cmb_network[nx_name]
    cmb_network['RANKAVG'] = cmb_network[network_names].sum(axis=1)
    cmb_network.loc[:, 'RANKAVG'] = cmb_network['RANKAVG'].div(float(sum(network_weights)))
    cmb_network = cmb_network.loc[cmb_network.RANKAVG <= max_avg, : ]
    cmb_network.sort_values('RANKAVG', inplace=True)
    return cmb_network

def main(network_files: List[str], network_names: str,
        network_weights: str, max_avg: int, out_file: str) -> bool:
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
    combine_df = consensus_network(network_files, network_names, network_weights, max_avg)
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
    PARSER.add_argument("-x", "--max_avg", type=float, default=1500000.0,
                        help="""maximum average rank allowed """)
    PARSER.add_argument("-o", "--out_file",
                        type=argparse.FileType(mode='w'), required=True,
                        help="output file in tab-seperated format")
    ARGS = PARSER.parse_args()
    print("""
       ARG : network_files : %s
       ARG : network_names : %s
       ARG : network_weights : %s
       ARG : max_avg : %s
       ARG : out_file : %s """ %
          (str(ARGS.network_files), str(ARGS.network_names), str(ARGS.network_weights),
           str(ARGS.max_avg), str(ARGS.out_file)))
    if not main(ARGS.network_files, ARGS.network_names, ARGS.network_weights,
                ARGS.max_avg, ARGS.out_file):
        PARSER.print_usage()
