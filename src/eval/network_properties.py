import argparse
import pandas as pd
import networkx as nx
from typing import List
import data_utils as du


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
                   max_edges: int = None) -> List[str]:
    net_df: pd.DataFrame = select_edges(du.load_reveng_network(network_file),
                                        wt_attr_name, max_edges)
    rev_net: nx.Graph = nx.from_pandas_edgelist(net_df, edge_attr=wt_attr_name)
    no_nodes = nx.number_of_nodes(rev_net)
    no_edges = nx.number_of_edges(rev_net)
    net_props = [no_nodes,
                 no_edges,
                 nx.density(rev_net),
                 nx.number_connected_components(rev_net),
                 2.0*no_edges/no_nodes]
    return [str(x) for x in net_props]

def graph_properties(network_files: str, wt_attr_name: str,
                     max_edges: int, out_file) -> None:
    cnames = ["Nodes", "Edges", "Density", "CC", "Avg.Degree"]
    prop_data = {fx : net_properties(fx, wt_attr_name, max_edges) for fx in network_files}
    prop_df = pd.DataFrame(data=prop_data, index=cnames)
    print(prop_df.transpose().to_csv(sep="\t", index=True))
    if out_file:
        prop_df.transpose().to_csv(out_file, sep="\t", index=True)

def net_node_properties(network_file: str, probe_map, wt_attr_name: str = 'wt',
                        max_edges: int = None) -> List[str]:
    net_df: pd.DataFrame = select_edges(du.load_reveng_network(network_file),
                                        wt_attr_name, max_edges)
    rev_net: nx.Graph = nx.from_pandas_edgelist(net_df, edge_attr=wt_attr_name)
    net_props = {
        'PROBE'      : {x: x for x in rev_net},
        'DEGREE'     : {x: nx.degree(rev_net, x) for x in rev_net},
        'BETWEENESS' : nx.betweenness_centrality(rev_net)
    }
    if probe_map is not None:
        net_props['ID'] = {x: probe_map[x] for x in rev_net}
    return net_props


def graph_node_properties(network_files: str, wt_attr_name: str, max_edges: int,
                          probe_file: str, out_file) -> None:
    if probe_file is not None:
        probe_df = du.load_annotation(probe_file)
        probe_map = {x: y for x, y in zip(probe_df.PROBE, probe_df.ID)
                     if x != y and y != 'no_match'}
    else:
        probe_map = None
    prop_data = net_node_properties(network_files[0], probe_map, wt_attr_name,
                                    max_edges)
    prop_df = pd.DataFrame(data=prop_data)
    if out_file:
        prop_df.to_csv(out_file, sep="\t", index=True)


if __name__ == "__main__":
    PROG_DESC = """Describe network statistics"""
    ARGPARSER = argparse.ArgumentParser(description=PROG_DESC)
    ARGPARSER.add_argument("network_file", nargs="+")
    ARGPARSER.add_argument("-x", "--max_edges", type=int,
                           help="""maximum number of edges""")
    ARGPARSER.add_argument("-d", "--node_prop", action='store_true',
                           help="""Order the edges ascending order""")
    ARGPARSER.add_argument("-t", "--wt_attr", type=str, default='wt',
                           help="name of weight attribute")
    ARGPARSER.add_argument("-p", "--probe_map_file", type=str,
                           help="Probe Mapping File")
    ARGPARSER.add_argument("-o", "--out_file",
                           type=argparse.FileType(mode='w'), required=True,
                           help="output file in tab-seperated format")
    CMDARGS = ARGPARSER.parse_args()
    if CMDARGS.node_prop is True:
        graph_node_properties(CMDARGS.network_file, CMDARGS.max_edges, CMDARGS.wt_attr,
                              CMDARGS.probe_map_file, CMDARGS.output_file)
    else:
        graph_properties(CMDARGS.network_file, CMDARGS.max_edges, CMDARGS.wt_attr,
                         CMDARGS.output_file)
