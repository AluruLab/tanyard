
import argparse
import matplotlib
matplotlib.use('Agg')
import venn
from data_utils import load_reveng_network

def get_pair(rcd):
    px = (rcd['source'], rcd['target'])
    if px[0] < px[1]:
        return px
    return (px[1], px[0])

def get_edge_lists(ndf):
    return [get_pair(rcd) for rcd in ndf[:, ['source', 'target']].to_dict('records')]

def make_network_names(network_names):
    if network_names:
        network_names = network_names.split(",")
        if len(network_names) > 4:
            network_names = network_names[0:4]
    if network_names is None:
        network_names = ['net_' + str(ix) for ix in range(len(network_files))]
    return network_names

def main(network_files: str, network_names: str, out_file: str) -> None:
    if len(network_files) > 4:
        network_files = network_files[0:4]
    else:
        return False
    network_names = get_network_names(network_names)
    network_dfs = [load_reveng_network(nx) for nx in network_files]
    st_lists = [get_edge_lists for df in network_dfs]
    venn_labels = venn.get_labels(st_lists)
    print(venn_labels)
    fig, ax = venn.venn4(venn_labels, names=network_names)
    fig.savefig(out_file)
    return True


if __name__ == "__main__":
    PROG_DESC = """Network venn diagram"""
    ARGPARSER = argparse.ArgumentParser(description=PROG_DESC)
    ARGPARSER.add_argument("-n", "--network_names", type=str,
                        help="""comma seperated names of the network;
                                should have as many names as the number of networks""")
    ARGPARSER.add_argument("-o", "--out_file",
                        type=str, required=True,
                        help="output file in png format")
    ARGPARSER.add_argument("network_files", nargs="+")
    CMDARGS = ARGPARSER.parse_args()
    print("""
       ARG : network_files : %s
       ARG : network_names : %s
       ARG : out_file : %s """ %
          (str(CMDARGS.network_files), str(CMDARGS.network_names),
           str(CMDARGS.out_file)))
    if not main(CMDARGS.network_files, CMDARGS.network_names, CMDARGS.out_file):
        ARGPARSER.print_usage()
