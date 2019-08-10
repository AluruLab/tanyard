from typing import Iterable, Set, Any, List
import argparse
import pandas as pd
from data_utils import load_annotation, load_reveng_network, load_gsnetwork, map_probes_cols


def common_network(network_files: Iterable[str]):
    common_df: pd.DataFrame = pd.DataFrame()
    for net_file in network_files:
        rv_net: pd.DataFrame = pd.read_csv(net_file, sep=r'\s+')
        if common_df.empty or common_df.shape[0] == 0:
            mcols = set(rv_net.columns) - set(['wt'])
            print(mcols)
            common_df = rv_net.loc[:, list(mcols)]
        else:
            mcols = set(rv_net.columns) - set(['wt'])
            rv_net = rv_net.loc[:, list(mcols)]
            common_df = common_df.merge(rv_net, how='inner')
    distinct_dfs = []
    for net_file in network_files:
        rv_net: pd.DataFrame = pd.read_csv(net_file, sep=r'\s+')
        rv_net = rv_net.loc[:, list(common_df.columns)]
        set_diff_df = pd.concat([common_df, rv_net]).drop_duplicates(keep=False)
        distinct_dfs.append(set_diff_df)
    return (common_df, distinct_dfs)


def main(network_files: Iterable[str], out_file: str,
         dist_files: List[str]) -> None:
    (common_df, distinct_dfs) = common_network(network_files)
    common_df.to_csv(out_file, sep='\t', index=False)
    for dist_df, dist_file in zip(distinct_dfs, dist_files):
        dist_df.to_csv(dist_file, sep="\t", index=False)


if __name__ == "__main__":
    PROG_DESC = """
    Finds a common  network as the intersection of input networks.
    Network intersection is computed as the intersection of edges of the input networks.
    """
    PARSER = argparse.ArgumentParser(description=PROG_DESC)
    PARSER.add_argument("network_files", nargs="+",
                        help="""network build from a reverse engineering methods
                                (currenlty supported: eda, tsv, adj)""")
    PARSER.add_argument("-o", "--out_file",
                        type=argparse.FileType(mode='w'), required=True,
                        help="output file in tab-seperated format")
    PARSER.add_argument("-n", "--network_distinct_files", type=str,
                        help="""comma seperated names of the output files;
                                should have as many names as the number of networks""")
    ARGS = PARSER.parse_args()
    distinct_file_names = ARGS.network_distinct_files
    if distinct_file_names:
        distinct_file_names = distinct_file_names.split(",")
        if len(distinct_file_names) != len(ARGS.network_files):
            distinct_file_names = None
    else:
        distinct_file_names = None
    if not distinct_file_names:
        distinct_file_names = ['net_' + str(ix) for ix in range(len(ARGS.network_files))]
    main(ARGS.network_files, ARGS.out_file, distinct_file_names)
