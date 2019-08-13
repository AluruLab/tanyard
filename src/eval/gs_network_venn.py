
import argparse
from typing import List
import networkx as nx
import matplotlib
import venn
import data_utils as du
matplotlib.use('Agg')

def get_pair(rcd):
    prx = (rcd['source'], rcd['target'])
    if prx[0] < prx[1]:
        return prx
    return (prx[1], prx[0])

def get_edge_lists(ndf):
    return [get_pair(rcd) for rcd in ndf[:, ['source', 'target']].to_dict('records')]

def get_network_names(network_names, nf_len):
    if network_names:
        network_names = network_names.split(",")
        if len(network_names) > nf_len:
            network_names = network_names[0:nf_len]
    if (not network_names) or len(network_names) < nf_len:
        network_names = ['net_' + str(ix) for ix in range(nf_len)]
    return network_names

def shortest_path(net_graph, src, tgt):
    spath = None
    if src > tgt:
        tmx = tgt
        tgt = src
        src = tmx
    if src in net_graph and tgt in net_graph:
        try:
            spath = (nx.shortest_path(net_graph, src, tgt), src, tgt)
        except nx.NetworkXNoPath:
            spath = (None, src, tgt)
    else:
        if src in net_graph:
            spath = (None, src, None)
        elif tgt in net_graph:
            spath = (None, None, tgt)
        else:
            spath = (None, None, None)
    return spath

def gs_shortest_paths(rv_net, gs_net, gs_nodes):
    rv_net_graph = nx.from_pandas_edgelist(rv_net, edge_attr='wt')
    # gs_common_nodes = sum((1 if x in rv_net_graph else 0 for x in gs_nodes))
    # gs_common_edges = sum((1 if (x in rv_net_graph and y in rv_net_graph) else 0
    #                        for x, y in zip(gs_net.TFPROBE, gs_net.TARGETPROBE)))
    gs_spath = [shortest_path(rv_net_graph, x, y)
                for x, y in zip(gs_net.TFPROBE, gs_net.TARGETPROBE)]
    return gs_spath

def filter_spath_list(gs_spath, max_dist):
    ft_lst = []
    for spx, udx, vdx in gs_spath:
        if spx is not None and (len(spx)-1) <= max_dist:
            ft_lst.append((udx, vdx))
    return ft_lst

def main(annot_file: str, gs_file: str, network_files: List[str],
         network_names: str, max_dist: int, out_file: str) -> None:
    if len(network_files) < 2:
        return False
    if len(network_files) >= 4:
        network_files = network_files[0:4]
    network_names = get_network_names(network_names, len(network_files))
    annot_df = du.load_annotation(annot_file)
    gs_net = du.map_probes(du.load_gsnetwork(gs_file), annot_df)
    gs_nedges = gs_net.shape[0]
    gs_nodes = set(gs_net.TFPROBE) | set(gs_net.TARGETPROBE)
    #
    network_dfs = [du.load_reveng_network(nx) for nx in network_files]
    st_lists = [gs_shortest_paths(df, gs_net, gs_nodes) for df in network_dfs]
    #
    ft_lists = [filter_spath_list(stx, max_dist) for stx in st_lists]
    venn_labels = venn.get_labels(ft_lists)
    print(venn_labels)
    if out_file:
        if len(network_files) == 4:
            fig, _ = venn.venn4(venn_labels, names=network_names)
            fig.savefig(out_file)
        if len(network_files) == 3:
            fig, _ = venn.venn3(venn_labels, names=network_names)
            fig.savefig(out_file)
        if len(network_files) == 2:
            fig, _ = venn.venn2(venn_labels, names=network_names)
            fig.savefig(out_file)
    return True


if __name__ == "__main__":
    PROG_DESC = """Network venn diagram of genes found in three networks"""
    ARGPARSER = argparse.ArgumentParser(description=PROG_DESC)
    ARGPARSER.add_argument("-n", "--network_names", type=str,
                           help="""comma seperated names of the network;
                                   should have as many names as the number of networks""")
    ARGPARSER.add_argument("-o", "--out_file",
                           type=str,
                           help="output file in png format")
    ARGPARSER.add_argument("-d", "--dist", type=int, default=1,
                           help="max. number of hops allowed")
    ARGPARSER.add_argument("annotation_file",
                           help="""annotation file
                                   (a tab seperated file mapping probe to ids)""")
    ARGPARSER.add_argument("gs_network_file",
                           help="""gold standard network
                                   (tab seperated file of TF-TARGET interactions)""")
    ARGPARSER.add_argument("network_files", nargs="+")
    CMDARGS = ARGPARSER.parse_args()
    print("""
       ARG : annotation_file : %s
       ARG : gs_network_file : %s
       ARG : network_files   : %s
       ARG : network_names   : %s
       ARG : max dist        : %s
       ARG : out_file        : %s """ %
          (str(CMDARGS.annotation_file), str(CMDARGS.gs_network_file),
           str(CMDARGS.network_files), str(CMDARGS.network_names),
           str(CMDARGS.dist), str(CMDARGS.out_file)))
    if not main(CMDARGS.annotation_file, CMDARGS.gs_network_file,
                CMDARGS.network_files, CMDARGS.network_names,
                CMDARGS.dist, CMDARGS.out_file):
        ARGPARSER.print_usage()
