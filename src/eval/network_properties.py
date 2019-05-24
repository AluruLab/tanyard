import argparse
import pandas as pd
import networkx as nx
from typing import List
from data_utils import load_reveng_network


def net_properties(network_file: str) -> List[str]:
    net_df: pd.DataFrame = load_reveng_network(network_file)
    rev_net: nx.Graph = nx.from_pandas_edgelist(net_df, edge_attr='wt')
    net_props =  [nx.number_of_nodes(rev_net),
                  nx.number_of_edges(rev_net),
                  nx.density(rev_net)]
    return [str(x) for x in net_props]

def main(network_files: str) -> None:
    cnames = ["Nodes", "Edges", "Density"]
    prop_data = {fx : net_properties(fx) for fx in network_files}
    prop_df = pd.DataFrame(data=prop_data, index=cnames)
    print(prop_df.to_csv(sep="\t", index=False))


if __name__ == "__main__":
    PROG_DESC = """Describe network statistics"""
    ARGPARSER = argparse.ArgumentParser(description=PROG_DESC)
    ARGPARSER.add_argument("network_file", nargs="+")
    CMDARGS = ARGPARSER.parse_args()
    main(CMDARGS.network_file)
