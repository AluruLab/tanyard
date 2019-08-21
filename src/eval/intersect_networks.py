from typing import Iterable, List
import argparse
import pandas as pd


def common_network(network_files: Iterable[str]):
    common_df: pd.DataFrame = pd.DataFrame()
    common_cols = set([])
    for net_file in network_files:
        rv_net: pd.DataFrame = pd.read_csv(net_file, sep=r'\s+')
        if common_df.empty or common_df.shape[0] == 0:
            mcols = set(rv_net.columns) - set(['wt'])
            common_df = rv_net.loc[:, list(mcols)]
            print(mcols)
            common_cols = set(mcols)
        else:
            mcols = set(rv_net.columns) - set(['wt'])
            rv_net = rv_net.loc[:, list(mcols)]
            common_df = common_df.merge(rv_net, how='inner')
            print(mcols)
            common_cols = common_cols & set(mcols)
    distinct_dfs = []
    common_cols = common_cols - set(['wt'])
    for net_file in network_files:
        rv_net: pd.DataFrame = pd.read_csv(net_file, sep=r'\s+')
        set_diff_df = pd.concat([common_df.loc[:, common_cols], 
            rv_net.loc[:, common_cols]]).drop_duplicates(keep=False)
        distinct_dfs.append(set_diff_df)
    return (common_df, distinct_dfs)


def main(network_files: Iterable[str], out_file: str,
         dist_files: List[str]) -> None:
    (common_df, distinct_dfs) = common_network(network_files)
    common_df.to_csv(out_file, sep='\t', index=False)
    if DISTINCT_FILE_NAMES is not None:
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
    DISTINCT_FILE_NAMES = ARGS.network_distinct_files
    if DISTINCT_FILE_NAMES:
        DISTINCT_FILE_NAMES = DISTINCT_FILE_NAMES.split(",")
        if len(DISTINCT_FILE_NAMES) != len(ARGS.network_files):
            DISTINCT_FILE_NAMES = None
    else:
        DISTINCT_FILE_NAMES = None
    #if not DISTINCT_FILE_NAMES:
    #    DISTINCT_FILE_NAMES = ['net_' + str(ix) for ix in range(len(ARGS.network_files))]
    main(ARGS.network_files, ARGS.out_file, DISTINCT_FILE_NAMES)
