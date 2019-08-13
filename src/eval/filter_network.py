import argparse
import pandas as pd
import numpy as np

def select_edges(net_df: pd.DataFrame, wt_attr_name: str = 'wt',
                 min_weight: float = None):
    if min_weight is None or min_weight < np.min(net_df[wt_attr_name]):
        return net_df
    return net_df.loc[net_df[wt_attr_name] >= min_weight]

def main(network_file, output_file, min_weight):
    net_df = pd.read_csv(network_file, sep="\t")
    print(net_df.columns)
    print(net_df.dtypes)
    print(min_weight)
    out_df = select_edges(net_df, min_weight=min_weight)
    out_df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    PROG_DESC = """Compute network with edge weights above given Threshold"""
    ARGPARSER = argparse.ArgumentParser(description=PROG_DESC)
    ARGPARSER.add_argument("network_file")
    ARGPARSER.add_argument("-t", "--min_weight", type=float, default=0.0,
                           help="""minimum weight threshold for the network""")
    ARGPARSER.add_argument("output_file")
    CMDARGS = ARGPARSER.parse_args()
    main(CMDARGS.network_file, CMDARGS.output_file, CMDARGS.min_weight)
