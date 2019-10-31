import argparse
import pandas as pd
import networkx as nx
from typing import List
from data_utils import load_reveng_network


def select_edges(net_df: pd.DataFrame, wt_attr_name: str = 'wt',
                 max_edges: int = None):
    if max_edges is None or max_edges >= net_df.shape[0]:
        return net_df
    cur_cols = net_df.columns
    maxwt_attr_name = wt_attr_name + '_max'
    net_df[maxwt_attr_name] = net_df[wt_attr_name].abs()
    net_df = net_df.nlargest(n=max_edges, columns=maxwt_attr_name)
    #print(net_df.iloc[0, :])
    return net_df.loc[:, cur_cols]

def net_properties(network_file: str, wt_attr_name: str = 'wt',
                   max_edges:int = None) -> List[str]:
    net_df: pd.DataFrame = select_edges(load_reveng_network(network_file), 
                                        wt_attr_name, max_edges)
    rev_net: nx.Graph = nx.from_pandas_edgelist(net_df, edge_attr=wt_attr_name)
    no_nodes = nx.number_of_nodes(rev_net)
    no_edges = nx.number_of_edges(rev_net)
    net_props =  [no_nodes,
                  no_edges,
                  nx.density(rev_net),
                  nx.number_connected_components(rev_net),
                  2.0*no_edges/no_nodes]
    return [str(x) for x in net_props]

def main(network_files: str, max_edges:int ) -> None:
    cnames = ["Nodes", "Edges", "Density", "CC", "Avg.Degree"]
    prop_data = {fx : net_properties(fx, 'wt', max_edges) for fx in network_files}
    prop_df = pd.DataFrame(data=prop_data, index=cnames)
    print(prop_df.transpose().to_csv(sep="\t", index=True))


if __name__ == "__main__":
    PROG_DESC = """Describe network statistics"""
    ARGPARSER = argparse.ArgumentParser(description=PROG_DESC)
    ARGPARSER.add_argument("network_file", nargs="+")
    ARGPARSER.add_argument("-x", "--max_edges", type=int,
                        help="""maximum number of edges""")
    CMDARGS = ARGPARSER.parse_args()
    main(CMDARGS.network_file, CMDARGS.max_edges)
