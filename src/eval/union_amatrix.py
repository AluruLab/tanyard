from typing import Iterable, List
import argparse
import pandas as pd
import numpy as np

def read_index_names(net_file: str) -> List[str]:
    with open(net_file) as net_fptr:
        for header_line in net_fptr:
            return header_line.strip().replace('"', '').split()


def union_index_names(network_files: Iterable[str]) -> List[str]:
    union_net_columns = set([])
    for net_fx in network_files:
        net_columns = read_index_names(net_fx)
        union_net_columns.update(net_columns)
    return sorted(list(union_net_columns))


def main(network_files: Iterable[str], out_file: str) -> None:
    union_row_names = union_index_names(network_files)
    union_col_names = union_row_names
    nrows = len(union_row_names)
    ncols = len(union_col_names)
    print("NROWS %d NCOLS %d" % (nrows, ncols))
    union_df = pd.DataFrame(np.zeros((nrows, ncols)),
                            columns=union_col_names,
                            index=union_row_names)
    for net_fx in network_files:
        net_df = pd.read_csv(net_fx, sep="\t")
        net_col_names = net_df.columns
        net_row_names = net_df.index
        union_sub_df = union_df.loc[net_row_names, net_col_names]
        net_df_max = net_df[abs(net_df) > abs(union_sub_df)]
        union_sub_df_max = union_sub_df[abs(union_sub_df) >= abs(net_df)]
        union_df.loc[net_row_names, net_col_names] = union_sub_df_max.add(
            net_df_max, fill_value=0)
    union_df.to_csv(out_file, sep="\t")


if __name__ == "__main__":
    PROG_DESC = """
    Finds a union of input networks.
    Network union is computed as the union of edges of the input networks.
    Outputs a tab-seperated values with weights corresponding to each network
    in a sperate column, and maximum and average weight.
    """
    PARSER = argparse.ArgumentParser(description=PROG_DESC)
    PARSER.add_argument("network_files", nargs="+",
                        help="""network build from a reverse engineering methods
                                (currenlty supported: eda, adj, tsv)""")
    PARSER.add_argument("-o", "--out_file",
                        type=argparse.FileType(mode='w'), required=True,
                        help="output file in matrix format")
    ARGS = PARSER.parse_args()
    if main(ARGS.network_files, ARGS.out_file):
        PARSER.print_usage()
