import argparse
import matplotlib
import pandas as pd
import venn
import data_utils as du
matplotlib.use('Agg')

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
    stpair = (rcd['source'], rcd['target'])
    if stpair[0] < stpair[1]:
        return stpair
    return (stpair[1], stpair[0])

def get_edge_lists(ndf):
    return [get_pair(rcd) for rcd in ndf.loc[:, ['source', 'target']].to_dict('records')]

def get_network_names(network_names, network_files):
    if network_names:
        network_names = network_names.split(",")
        if len(network_names) > 6:
            network_names = network_names[0:6]
    if network_names is None:
        network_names = ['net_' + str(ix) for ix in range(len(network_files))]
    return network_names


def main(network_files: str, network_names: str, out_file: str,
         wt_attr: str = 'wt', max_edges: int = None,
         reverse_order: bool = False, percent: bool = True) -> None:
    if len(network_files) < 2:
        return False
    if len(network_files) > 6:
        network_files = network_files[0:6]
    network_names = get_network_names(network_names, network_files)
    network_dfs = [load_network(nx, wt_attr, max_edges, reverse_order) for nx in network_files]
    st_lists = [get_edge_lists(df) for df in network_dfs]
    if percent is True:
        fillv = ["number", "percent"]
    else:
        fillv = ["number"]
    venn_labels = venn.get_labels(st_lists, fill=fillv)
    print(venn_labels)
    if len(network_files) == 6:
        fig, _ = venn.venn6(venn_labels, names=network_names)
        fig.savefig(out_file)
    if len(network_files) == 5:
        fig, _ = venn.venn5(venn_labels, names=network_names)
        fig.savefig(out_file)
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
    ARGPARSER.add_argument("-p", "--percent", action='store_true',
                           help="""Added percentage in the venn diagram""")
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
       ARG : percent         : %s
       ARG : out_file : %s """ %
          (str(CMDARGS.network_files), str(CMDARGS.network_names),
           str(CMDARGS.wt_attr), str(CMDARGS.max_edges), str(CMDARGS.reverse_order),
           str(CMDARGS.percent), str(CMDARGS.out_file)))
    if not main(CMDARGS.network_files, CMDARGS.network_names,
                CMDARGS.out_file, CMDARGS.wt_attr, CMDARGS.max_edges,
                CMDARGS.reverse_order, CMDARGS.percent):
        ARGPARSER.print_usage()
