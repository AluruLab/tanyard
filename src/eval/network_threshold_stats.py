import argparse
import numpy as np
import pandas as pd
import networkx as nx


def abs_max(row_x):
    row_max = np.max(row_x)
    row_min = np.min(row_x)
    if row_max > abs(row_min):
        return row_max
    return row_min


def select_edges(net_df: pd.DataFrame, wt_attr_name: str = 'wt',
                 max_edges: int = None):
    if max_edges is None or max_edges >= net_df.shape[0]:
        return net_df
    cur_cols = net_df.columns
    maxwt_attr_name = wt_attr_name + '_max'
    net_df[maxwt_attr_name] = net_df[wt_attr_name].abs()
    net_df = net_df.nlargest(n=max_edges, columns=maxwt_attr_name)
    return net_df.loc[:, cur_cols]


def main(network_file, output_file):
    mn_edges = 1000000
    mx_edges = 3000000
    ndiv = 100
    net_df = pd.read_csv(network_file, sep="\t" )
    print(net_df.columns)
    net_df = select_edges(net_df, max_edges=mx_edges)
    df_edges = mx_edges - mn_edges
    th_lst = [float(df_edges*x)/ndiv for x in range(0, ndiv+1)]
    edge_ct = [0 for _ in range(len(th_lst))]
    node_ct = [0 for _ in range(len(th_lst))]
    density_ct = [0.0 for _ in range(len(th_lst))]
    minwt_ct = [0.0 for _ in range(len(th_lst))]
    for i,j in enumerate(th_lst):
        x = mn_edges + int(j)
        sub_net_df = select_edges(net_df, max_edges=x)
        rev_net: nx.Graph = nx.from_pandas_edgelist(sub_net_df, edge_attr='wt')
        node_ct[i] = nx.number_of_nodes(rev_net)
        edge_ct[i] = nx.number_of_edges(rev_net)
        density_ct[i] = nx.density(rev_net)
        minwt_ct[i] = np.min(sub_net_df['wt'])
    out_df = pd.DataFrame({"MAX_EDGES": th_lst,
                           "EDGE_COUNT": edge_ct,
                           "NODE_COUNT": node_ct,
                           "DENSITY": density_ct,
                           "MWT": minwt_ct})
    out_df.to_csv(output_file)


if __name__ == "__main__":
    PROG_DESC = """Compute network density vs Threshold"""
    ARGPARSER = argparse.ArgumentParser(description=PROG_DESC)
    ARGPARSER.add_argument("network_file")
    ARGPARSER.add_argument("output_file")
    CMDARGS = ARGPARSER.parse_args()
    main(CMDARGS.network_file, CMDARGS.output_file)
