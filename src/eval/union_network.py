import argparse
import numpy as np
import pandas as pd
from data_utils import load_reveng_network

def abs_max(row_x):
    row_max = np.max(row_x)
    row_min = np.min(row_x)
    if row_max > abs(row_min):
        return row_max
    return row_min

def combine_network(network_files, network_names=None):
    if network_names is None:
        network_names = ['wt_'+str(ix) for ix in range(len(network_files))]
    cmb_network = pd.DataFrame(columns=['source', 'target'])
    for nx_name, nx_file in zip(network_names, network_files):
        ndf = load_reveng_network(nx_file, nx_name)
        cmb_network = cmb_network.merge(ndf, how='outer', on=['source', 'target'])
        #print(str(nx_name), nx_file, ndf.shape, cmb_network.columns, cmb_network.shape)
    cmb_network['wt'] = cmb_network[network_names].apply(abs_max, axis=1)
    cmb_network['avgwt'] = cmb_network[network_names].mean(axis=1)
    return cmb_network


def main(network_names, network_files, out_file):
    if network_names:
        network_names = network_names.split(",")
        if len(network_names) == len(network_files):
            combine_df = combine_network(network_files, network_names)
        else:
            print("Length of network names should be equal to length of network files")
            return False
    else:
        combine_df = combine_network(network_files)
    combine_df.to_csv(out_file, sep='\t', index=False)
    return True


if __name__ == "__main__":
    PROG_DESC = """
    Finds a union of input networks.
    Network union is computed as the union of edges of the input networks.
    Outputs a tab-seperated values with weights corresponding to each network
    in a sperate column, and maximum and average weight.
    """
    PARSER = argparse.ArgumentParser(description=PROG_DESC)
    PARSER.add_argument("network_files", nargs="+",
                        help="""network build from a reverse engineering methods
                                (currenlty supported: eda, adj, tsv)""")
    PARSER.add_argument("-n", "--network_names", type=str,
                        help="""comma seperated names of the network;
                                should have as many names as the number of networks""")
    PARSER.add_argument("-o", "--out_file",
                        type=argparse.FileType(mode='w'), required=True,
                        help="output file in tab-seperated format")
    ARGS = PARSER.parse_args()
    if not main(ARGS.network_names, ARGS.network_files, ARGS.out_file):
        PARSER.print_usage()
