from typing import Iterable
import argparse
import pandas as pd

def main(network_files: Iterable[str], out_file: str) -> None:
    union_df = pd.DataFrame()
    grow_names = union_df.index
    gcol_names = union_df.columns
    for net_fx in network_files:
        if union_df.empty:
            union_df = pd.read_csv(net_fx, sep="\t")
            gcol_names = union_df.columns
            grow_names = union_df.index
        else:
            net_df = pd.read_csv(net_fx, sep="\t")
            net_df = net_df.loc[grow_names, gcol_names]
            net_df_max = net_df[abs(net_df) > abs(union_df)]
            union_df_max = union_df[abs(union_df) >= abs(net_df)]
            union_df = union_df_max.add(net_df_max, fill_value=0)
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
