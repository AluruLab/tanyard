import argparse
import matplotlib
matplotlib.use('Agg')
import venn
import data_utils as du
import pandas as pd

def select_edges(net_df: pd.DataFrame, wt_attr_name: str = 'wt',
                 max_edges: int = None,
                 reverse_order: bool = False):
    if max_edges is None or max_edges >= net_df.shape[0]:
        return net_df
    cur_cols = net_df.columns
    maxwt_attr_name = wt_attr_name + '_max'
    net_df[maxwt_attr_name] = net_df[wt_attr_name].abs()
    if reverse_order is True:
       net_df = net_df.nsmallest(n=max_edges, columns=maxwt_attr_name)
    else:
       net_df = net_df.nlargest(n=max_edges, columns=maxwt_attr_name)
    #print(net_df.iloc[0, :])
    return net_df.loc[:, cur_cols]


def load_network(net_file: str, wt_attr: str = 'wt',
                 max_edges: int = None,
                 reverse_order: bool = False):
    return select_edges(du.load_reveng_network(net_file, wt_attr_name=wt_attr),
                        wt_attr_name=wt_attr,
                        max_edges=max_edges,
                        reverse_order=reverse_order)

def get_pair(rcd):
    px = (rcd['source'], rcd['target'])
    if px[0] < px[1]:
        return px
    return (px[1], px[0])

def get_edge_lists(ndf):
    return [get_pair(rcd) for rcd in ndf.loc[:, ['source', 'target']].to_dict('records')]

def get_network_names(network_names):
    if network_names:
        network_names = network_names.split(",")
        if len(network_names) > 4:
            network_names = network_names[0:4]
    if network_names is None:
        network_names = ['net_' + str(ix) for ix in range(len(network_files))]
    return network_names


def main(network_files: str, network_names: str, out_file: str,
         wt_attr: str='wt', max_edges: int=None, reverse_order: bool=False) -> None:
    if len(network_files) < 2:
        return False
    if len(network_files) > 4:
        network_files = network_files[0:4]
    network_names = get_network_names(network_names)
    network_dfs = [load_network(nx, wt_attr, max_edges, reverse_order) for nx in network_files]
    st_lists = [get_edge_lists(df) for df in network_dfs]
    venn_labels = venn.get_labels(st_lists)
    print(venn_labels)
    if len(network_files) == 4:
        fig, ax = venn.venn4(venn_labels, names=network_names)
        fig.savefig(out_file)
    if len(network_files) == 3:
        fig, _ = venn.venn3(venn_labels, names=network_names)
        fig.savefig(out_file)
    if len(network_files) == 2:
        fig, _ = venn.venn2(venn_labels, names=network_names)
        fig.savefig(out_file)
    return True


if __name__ == "__main__":
    PROG_DESC = """Network venn diagram"""
    ARGPARSER = argparse.ArgumentParser(description=PROG_DESC)
    ARGPARSER.add_argument("-n", "--network_names", type=str,
                        help="""comma seperated names of the network;
                                should have as many names as the number of networks""")
    ARGPARSER.add_argument("-t", "--wt_attr", type=str, default='wt',
                           help="name of weight attribute")
    ARGPARSER.add_argument("-x", "--max_edges", type=int,
                           help="""Maximum number of edges""")
    ARGPARSER.add_argument("-r", "--reverse_order", action='store_true',
                           help="""Order the edges ascending order""")
    ARGPARSER.add_argument("-o", "--out_file",
                        type=str, required=True,
                        help="output file in png format")
    ARGPARSER.add_argument("network_files", nargs="+")
    CMDARGS = ARGPARSER.parse_args()
    print("""
       ARG : network_files : %s
       ARG : network_names : %s
       ARG : wt_attr : %s
       ARG : max_edges : %s
       ARG : reverse_order : %s
       ARG : out_file : %s """ %
          (str(CMDARGS.network_files), str(CMDARGS.network_names),
           str(CMDARGS.wt_attr), str(CMDARGS.max_edges), str(CMDARGS.reverse_order),
           str(CMDARGS.out_file)))
    if not main(CMDARGS.network_files, CMDARGS.network_names, CMDARGS.out_file,
                CMDARGS.wt_attr, CMDARGS.max_edges, CMDARGS.reverse_order):
        ARGPARSER.print_usage()
