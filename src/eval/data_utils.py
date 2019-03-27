import pandas as pd
import networkx as nx

def load_full_annotation(annot_file):
    """
    Read annotation file (a tab seperated file), and load the annotation
    data frame.
    There can be multiple IDs associated with same probe, which is proveided
    in the input file delimited by ';'
    """
    return pd.read_csv(annot_file, sep="\t", header=0)


def load_annotation(annot_file):
    """
    Read annotation file (a tab seperated file), and load the annotation
    data frame.
    Includes two columns ID and PROBE
    There can be multiple IDs associated with same probe, which is proveided
    in the input file delimited by ';'
    """
    input_df = pd.read_csv(annot_file, sep="\t", header=0)
    annot_df = pd.DataFrame(input_df.ID.str.split(';').tolist(),
                            index=input_df.PROBE).stack()
    annot_df = annot_df.reset_index()[[0, 'PROBE']]
    annot_df.columns = ['ID', 'PROBE']
    return annot_df


def load_gsnetwork(gs_file):
    """
    TAB seperated file w. header TF,TARGET
    """
    return pd.read_csv(gs_file, sep="\t")


def load_eda_network(eda_file):
    """
    Load network from eda file. Eda file lists the edges in the following format



    Parameters
    ----------

    """
    tmp_df =  pd.read_csv(eda_file, sep=' ', usecols=[0,2,4], skiprows=[0],
                          names=['source','target','wt'] )
    tmp_rcds = [(x, y, z) if x < y else (y, x, z) for x, y, z in tmp_df.to_records()]
    return pd.DataFrame(tmp_rcds, columns=['source','target','wt'])


def load_adj_network(adj_file, comments=["#", ">"], delimiter="\t"):
    """
    Load network with adjacency list;
    Adjacency file lists the target lists for each source node:
    source_node target_node1 weight1 ....

    Parameters
    ----------
    comments : list of comment indicators

    delimiter : seperators between the fields
    """
    edge_list = []
    with open(adj_file) as adjfp:
        for line in adjfp:
            if any((True if line.startswith(cx) else False for cx in comments)):
                continue
            vlist=line.strip().split(delimiter)
            if len(vlist) >= 3:
                tgt_length = (len(vlist) - 1) / 2
                edge_list.append([(vlist[0], vlist[2*ix+1], float(vlist[2*ix+2])) for ix in range(tgt_length)])
    edge_list = [(x, y, z) if x < y else (y, x, z) for x, y, z in edge_list]
    return pd.DataFrame(edge_list, columns=['source','target','wt'])


def load_reveng_network(net_file):
    if net_file.endswith(".eda"):
        return load_eda_network(net_file)
    elif net_file.endswith(".adj"):
        return load_adj_network(net_file)
    else:
        return None

def map_probes(gs_net, annot_df):
    annot_filter_df = annot_df.loc[annot_df.ID != 'no_match', :]
    gs_net_mapped = gs_net.merge(annot_filter_df, left_on='TF', right_on='ID')
    gs_net_mapped.columns = ['TF', 'TARGET', 'TFID', 'TFPROBE']
    gs_net_mapped = gs_net_mapped.merge(annot_filter_df, left_on='TARGET',
                                        right_on='ID')
    gs_net_mapped.columns = ['TF', 'TARGET', 'TFID', 'TFPROBE', 'TARGETID', 'TARGETPROBE']
    return gs_net_mapped
