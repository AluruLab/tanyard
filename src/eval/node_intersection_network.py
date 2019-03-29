import argparse
import pandas as pd
import networkx as nx
import argparse
from data_utils import load_annotation, load_reveng_network, load_gsnetwork, map_probes


def common_network(annot_file, gs_file, network_files):
    annot_df = load_annotation(annot_file)
    gsnet_df = load_gsnetwork(gs_file)
    gs_net = map_probes(gsnet_df, annot_df)
    common_nodes = set(gs_net.TFPROBE) | set(gs_net.TARGETPROBE)
    for net_file in network_files:
        rv_net = load_reveng_network(net_file)
        rv_net_nodes = set(rv_net.source) | set(rv_net.target)
        common_nodes = common_nodes & rv_net_nodes
    common_df = gs_net.loc[(gs_net.TFPROBE.isin(common_nodes)) & (gs_net.TARGETPROBE.isin(common_nodes)) , :]
    return common_df.loc[:, ['TF', 'TARGET']]


if __name__ == "__main__":
    prog_desc = """
    Finds a common gold standard network as the intersection of input network nodes.
    Network intersection is computed as the intersection of nodes of the input networks,
    in common with the gold standard networks, after mapping to the annotation.
    """
    parser = argparse.ArgumentParser(description=prog_desc)
    parser.add_argument("annotation_file",
        help="annotation file (a tab seperated file) mapping probe to ids")
    parser.add_argument("gs_network_file",
        help="gold standard network (tab seperated file of TF-TARGET interactions)")
    parser.add_argument("network_files", nargs="+",
        help="network build from a reverse engineering methods (currenlty supported: eda, tsv, adj)")
    parser.add_argument("-o", "--out_file", type=argparse.FileType(mode='w'), required=True,
                        help="output file in tab-seperated format")
    args = parser.parse_args()
    common_df = common_network(args.annotation_file, args.gs_network_file, args.network_files)
    common_df.to_csv(args.out_file, sep='\t', index=False)
