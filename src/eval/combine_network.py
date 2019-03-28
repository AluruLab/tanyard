import pandas as pd 
import argparse
from data_utils import load_annotation, load_reveng_network

def combine_network(network_files, network_names=None):
    if network_names is None:
        network_names = ['wt_'+str(ix) for ix in range(len(network_files))]
    cmb_network = pd.DataFrame(columns=['source','target'])
    for nx_name, nx_file in zip(network_names, network_files):
        ndf = load_reveng_network(nx_file, nx_name)
        cmb_network = cmb_network.merge(ndf, how='outer', on=['source', 'target'])
        #print(str(nx_name), nx_file, ndf.shape, cmb_network.columns, cmb_network.shape)
    cmb_network['maxwt'] = cmb_network[network_names].max(axis=1)
    cmb_network['wt'] = cmb_network[network_names].mean(axis=1)
    return cmb_network


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("network_files", nargs="+",
        help="network build from a reverse engineering methods (currenlty supported: eda)")
    parser.add_argument("-n", "--network_names", type=str,
                        help="names of the network")
    parser.add_argument("-o", "--out_file", type=str,
                        help="output file in tab-seperated format")
    args = parser.parse_args()
    if args.network_names:
       network_files = args.network_files
       network_names = args.network_names.split(",")
       if len(network_names) == len(network_files):
           combine_df = combine_network(network_files, network_names)
       else:
           print("Length of network names should be equal to length of network files")
    else:
       combine_df = combine_network(args.network_files)
    if args.out_file:
        combine_df.to_csv(args.out_file, sep='\t', index=False)
