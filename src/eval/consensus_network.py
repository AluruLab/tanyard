import pandas as pd
import numpy as np
from typing import List
from data_utils import load_reveng_network
import argparse

def count_nonnan(row_x):
    return np.count_nonzero(~np.isnan(row_x))


def consensus_network(network_files, network_names=None, max_networks: int = None):
    if network_names is None:
        network_names = ['net_' + str(ix) for ix in range(len(network_files))]
    cmb_network = pd.DataFrame(columns=['source', 'target'])
    for nx_name, nx_file in zip(network_names, network_files):
        ndf = load_reveng_network(nx_file, nx_name)
        cmb_network = cmb_network.merge(ndf, how='outer', on=['source', 'target'])
        #print(str(nx_name), nx_file, ndf.shape, cmb_network.columns, cmb_network.shape)
    cmb_network['count'] = cmb_network[network_names].apply(count_nonnan, axis=1)
    return cmb_network.loc[cmb_network.count > max_networks, : ]

def main(network_files: List[str], network_names: List[str],
         max_networks: int, out_file: str) -> bool:
    if network_names:
        network_names = network_names.split(",")
        if len(network_names) == len(network_files):
            combine_df = consensus_network(network_files, network_names, max_networks)
        else:
            print("Length of network names should be equal to length of network files")
            return False
    else:
        combine_df = consensus_network(network_files, None, max_networks)
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
    PARSER.add_argument("-x", "--max_networks", type=int,
                        help="""comma seperated names of the network;
                                how many networks """)
    PARSER.add_argument("-o", "--out_file",
                        type=argparse.FileType(mode='w'), required=True,
                        help="output file in tab-seperated format")
    ARGS = PARSER.parse_args()
    print("""
       ARG : network_files : %s
       ARG : network_names : %s
       ARG : max_edges : %s
       ARG : out_file : %s """ %
          (str(ARGS.network_files), str(ARGS.network_names),
           str(ARGS.max_edges), str(ARGS.out_file)))
    if not main(ARGS.network_names, ARGS.network_files,
                ARGS.out_file, ARGS.max_edges):
        PARSER.print_usage()
