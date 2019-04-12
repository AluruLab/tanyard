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
    common_df = gs_net.loc[(gs_net.TFPROBE.isin(common_nodes)) &
                           (gs_net.TARGETPROBE.isin(common_nodes)), :]
    return common_df.loc[:, ['TF', 'TARGET']]


def main(annotation_file, gs_network_file, network_files, out_file):
    common_df = common_network(annotation_file, gs_network_file, network_files)
    common_df.to_csv(out_file, sep='\t', index=False)

if __name__ == "__main__":
    PROG_DESC = """
    Finds a common gold standard network as the intersection of input network nodes.
    Network intersection is computed as the intersection of nodes of the input networks,
    in common with the gold standard networks, after mapping to the annotation.
    """
    PARSER = argparse.ArgumentParser(description=PROG_DESC)
    PARSER.add_argument("annotation_file",
                        help="annotation file (a tab seperated file) mapping probe to ids")
    PARSER.add_argument("gs_network_file",
                        help="gold standard network (tab seperated file of TF-TARGET interactions)")
    PARSER.add_argument("network_files", nargs="+",
                        help="network build from a reverse engineering methods"
                        "(currenlty supported: eda, tsv, adj)")
    PARSER.add_argument("-o", "--out_file", type=argparse.FileType(mode='w'), required=True,
                        help="output file in tab-seperated format")
    ARGS = PARSER.parse_args()
    main(ARGS.annotation_file, ARGS.gs_network_file,
         ARGS.network_files, ARGS.out_file)
