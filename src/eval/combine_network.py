import pandas as pd 
import argparse
from data_utils import load_annotation, load_reveng_network

def combine_network(network_files):
    for nx_file in network_files:
        ndf = load_reveng_network(nx_file)
        print(nx_file, ndf.shape)
    return ndf


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("network_files", nargs="+",
        help="network build from a reverse engineering methods (currenlty supported: eda)")
    parser.add_argument("-o", "--out_file", type=str,
                        help="output file in tab-seperated format")
    args = parser.parse_args()
    combine_df = combine_network(args.gs_network_file)
    if args.out_file:
        combine_df.to_csv(args.out_file, sep='\t', index=False)
