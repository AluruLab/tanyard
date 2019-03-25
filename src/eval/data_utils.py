import pandas as pd
import networkx as nx

def load_annotation(annot_file):
    """
    Read annotation file (a tab seperated file), and load the annotation
    data frame.
    There can be multiple IDs associated with same probe, which is proveided
    in the input file delimited by ';'
    """
    input_df = pd.read_csv(annot_file, sep="\t", header=None,
                           names=["PROBE", "ID"])
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
    Load network from eda file
    """
    return pd.read_csv(eda_file, sep = ' ', usecols = [0,2,4], skiprows = [0],
                       names = ['source','target','wt'] )

def load_reveng_network(net_file):
    if net_file.endswith(".eda"):
        return load_eda_network(net_file)
    return None

def map_probes(gs_net, annot_df):
    annot_filter_df = annot_df.loc[annot_df.ID != 'no_match', :]
    gs_net_mapped = gs_net.merge(annot_filter_df, left_on='TF', right_on='ID')
    gs_net_mapped.columns = ['TF', 'TARGET', 'TFID', 'TFPROBE']
    gs_net_mapped = gs_net_mapped.merge(annot_filter_df, left_on='TARGET',
                                        right_on='ID')
    gs_net_mapped.columns = ['TF', 'TARGET', 'TFID', 'TFPROBE', 'TARGETID', 'TARGETPROBE']
    return gs_net_mapped
