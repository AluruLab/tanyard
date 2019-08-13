import argparse
import pandas as pd
import numpy as np
from data_utils import load_reveng_network


def select_edges_ct(net_df: pd.DataFrame, wt_attr_name: str = 'wt',
                    max_edges: int = None):
    if max_edges is None or max_edges >= net_df.shape[0]:
        return net_df
    cur_cols = net_df.columns
    maxwt_attr_name = wt_attr_name + '_max'
    net_df[maxwt_attr_name] = net_df[wt_attr_name].abs()
    net_df = net_df.nlargest(n=max_edges, columns=maxwt_attr_name)
    return net_df.loc[:, cur_cols]


def select_edges_wt(net_df: pd.DataFrame, wt_attr_name: str = 'wt',
                    min_weight: float = None):
    if min_weight is None or min_weight < np.min(net_df[wt_attr_name]):
        return net_df
    return net_df.loc[net_df[wt_attr_name] >= min_weight]

def main(network_file, output_file, min_weight, max_edges):
    net_df = load_reveng_network(network_file)
    print(net_df.columns)
    print(net_df.dtypes)
    print(min_weight)
    print(max_edges)
    if min_weight:
        tmp_df = select_edges_wt(net_df, min_weight=min_weight)
    else:
        tmp_df = net_df
    if max_edges:
        out_df = select_edges_ct(tmp_df, max_edges=max_edges)
    else:
        out_df = tmp_df
    out_df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    PROG_DESC = """Compute network with edge weights above given Threshold"""
    ARGPARSER = argparse.ArgumentParser(description=PROG_DESC)
    ARGPARSER.add_argument("network_file")
    ARGPARSER.add_argument("-t", "--min_weight", type=float, default=0.0,
                        help="""minimum weight threshold for the network""")
    ARGPARSER.add_argument("-x", "--max_edges", type=int,
                        help="""maximum eges in the network""")
    ARGPARSER.add_argument("output_file")
    CMDARGS = ARGPARSER.parse_args()
    main(CMDARGS.network_file, CMDARGS.output_file, CMDARGS.min_weight, CMDARGS.max_edges)
