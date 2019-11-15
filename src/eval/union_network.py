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


def select_edges(net_df: pd.DataFrame, wt_attr_name: str = 'wt',
                 max_edges: int = None):
    if max_edges is None or max_edges >= net_df.shape[0]:
        return net_df
    cur_cols = net_df.columns
    maxwt_attr_name = wt_attr_name + '_max'
    net_df[maxwt_attr_name] = net_df[wt_attr_name].abs()
    net_df = net_df.nlargest(n=max_edges, columns=maxwt_attr_name)
    return net_df.loc[:, cur_cols]




def combine_network(network_files, network_names=None, max_edges: int = None, avg_wt=True, max_wt=True):
    if network_names is None:
        network_names = ['wt_'+str(ix) for ix in range(len(network_files))]
    cmb_network = pd.DataFrame(columns=['source', 'target'])
    for nx_name, nx_file in zip(network_names, network_files):
        ndf = select_edges(load_reveng_network(nx_file, nx_name), nx_name, max_edges)
        cmb_network = cmb_network.merge(ndf, how='outer', on=['source', 'target'])
        #print(str(nx_name), nx_file, ndf.shape, cmb_network.columns, cmb_network.shape)
    if max_wt is True:
        cmb_network['wt'] = cmb_network[network_names].apply(abs_max, axis=1)
    if avg_wt is True:
        cmb_network['avgwt'] = cmb_network[network_names].mean(axis=1)
    return cmb_network


def main(network_names, network_files, out_file, max_edges, avg_wt, max_wt):
    if network_names:
        network_names = network_names.split(",")
        if len(network_names) == len(network_files):
            combine_df = combine_network(network_files, network_names, max_edges, avg_wt, max_wt)
        else:
            print("Length of network names should be equal to length of network files")
            return False
    else:
        combine_df = combine_network(network_files, None, max_edges)
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
    PARSER.add_argument("-x", "--max_edges", type=int,
                        help="""max number of edges to output""")
    PARSER.add_argument("-g", "--no_wt_avg", action='store_false',
                        help="""compute the average wt. (default: True)""")
    PARSER.add_argument("-m", "--no_wt_max", action='store_false',
                        help="""compute the average max wt. (default: True)""")
    PARSER.add_argument("-o", "--out_file",
                        type=argparse.FileType(mode='w'), required=True,
                        help="output file in tab-seperated format")
    ARGS = PARSER.parse_args()
    print("""
       ARG : network_files : %s
       ARG : network_names : %s
       ARG : max_edges     : %s
       ARG : wt_avg        : %s
       ARG : wt_max        : %s
       ARG : out_file      : %s """ %
          (str(ARGS.network_files), str(ARGS.network_names),
           str(ARGS.max_edges), str(ARGS.no_wt_avg), str(ARGS.no_wt_max),
           str(ARGS.out_file)))
    if not main(ARGS.network_names, ARGS.network_files,
                ARGS.out_file, ARGS.max_edges, ARGS.no_wt_avg, ARGS.no_wt_max):
        PARSER.print_usage()
