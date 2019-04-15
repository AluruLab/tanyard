import argparse
import networkx as nx
from data_utils import load_reveng_network


def main(network_file: str):
    net_df = load_reveng_network(network_file)
    rev_net = nx.from_pandas_edgelist(net_df, edge_attr='wt')
    network_props = """ Nodes {} Edges {} Density {} """
    print(network_props % (nx.number_of_nodes(rev_net),
                           nx.number_of_edges(rev_net),
                           nx.density(rev_net)))


if __name__ == "__main__":
    PROG_DESC = """Describe network statistics"""
    ARGPARSER = argparse.ArgumentParser(description=PROG_DESC)
    ARGPARSER.add_argument("network_file")
    CMDARGS = ARGPARSER.parse_args()
    main(CMDARGS.network_file)
    