import pandas as pd
import networkx as nx
import argparse
from data_utils import load_annotation, load_reveng_network, load_gsnetwork, map_probes


def common_network(annot_file, gs_file, network_files):
    annot_df = load_annotation(annot_file)
    gsnet_df = load_gsnetwork(gs_file)
    gs_net = map_probes(gsnet_df, annot_df)
    common_df = None
    for net_file in network_files:
        rv_net = load_reveng_network(net_file)
        rv_net_nodes = set(rv_net.source) | set(rv_net.source)
        if common_df is None or common_df.shape[0] == 0:
            common_df = gs_net.loc[(gs_net.TFPROBE in rv_net_nodes) & (gs_net.TARGETPROBE in rv_net_nodes) , :]
        else:
            common_df = common_df.loc[(common_df.TFPROBE in rv_net_nodes) & (common_df.TARGETPROBE in rv_net_nodes) , :]
    print(common_df.shape)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("annotation_file",
        help="annotation file (a tab seperated file) mapping probe to ids")
    parser.add_argument("gs_network_file",
        help="gold standard network (tab seperated file of TF-TARGET interactions)")
    parser.add_argument("network_files", nargs="+",
        help="network build from a reverse engineering methods (currenlty supported: eda)")
    args = parser.parse_args()
    common_network(args.annotation_file, args.gs_network_file, args.network_files)
