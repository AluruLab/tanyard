from typing import List, Tuple
import numpy as np
import pandas as pd


def load_full_annotation(annot_file: str) -> pd.DataFrame:
    """
    Read annotation file (a tab seperated file), and load the annotation
    data frame.
    There can be multiple IDs associated with same probe, which is proveided
    in the input file delimited by ';'

    Parameters
    ----------
    annot_file : Annotation file (tab seperated)

    Returns
    -------
    pandas DataFrame with annotation

    """
    return pd.read_csv(annot_file, sep="\t", header=0)


def load_annotation(annot_file: str) -> pd.DataFrame:
    """
    Read annotation file (a tab seperated file), and load the annotation
    data frame.
    Includes two columns ID and PROBE
    There can be multiple IDs associated with same probe, which is proveided
    in the input file delimited by ';'

    Parameters
    ----------
    annot_file : Annotation file (tab seperated)

    Returns
    -------
    pandas DataFrame with expanded annotation
    """
    input_df = pd.read_csv(annot_file, sep="\t", header=0)
    annot_df = pd.DataFrame(input_df.ID.str.split(';').tolist(),
                            index=input_df.PROBE).stack()
    annot_df = annot_df.reset_index()[[0, 'PROBE']]
    annot_df.columns = ['ID', 'PROBE']
    return annot_df


def load_gsnetwork(gs_file: str) -> pd.DataFrame:
    """
    TAB seperated file w. header TF,TARGET

    Parameters
    ----------
    annot_file : Annotation file (tab seperated w. header 'TF' and 'TARGET)

    Returns
    -------
    pandas DataFrame with columns, TF, TARGET and weight (if provided)
    """
    tdf = pd.read_csv(gs_file, sep="\t")
    return tdf.loc[tdf.TF != tdf.TARGET, :]

def remove_dupe_rows(in_df: pd.DataFrame, wt_attr_name: str) -> pd.DataFrame:
    """
    Given a data frame of three columns : "source", "target" and wt_attr_name
    such that only the row with max weight is returned

    Parameters
    ----------
    in_df : Input data frame with columns 'source', 'target', wt_attr_name

    wt_attr_name : weight column name in the data frame returned

    Returns
    -------
    pandas DataFrame with three columns: 'source', 'target', wt_attr_name such that
    only the row with max weight is returned
    """
    in_df = in_df.sort_values(by=[wt_attr_name], ascending=False)
    return in_df.drop_duplicates(subset=['source', 'target'], keep='first')


def order_network_rows(in_df: pd.DataFrame, wt_attr_name: str) -> pd.DataFrame:
    """
    Given a data frame of three columns : "source", "target" and wt_attr_name
    returns rows such that source < target

    Parameters
    ----------
    in_df : Input data frame with columns 'source', 'target', wt_attr_name

    wt_attr_name : weight column name in the data frame returned

    Returns
    -------
    pandas DataFrame with three columns: 'source', 'target', wt_attr_name such that
    source entry < target entry
    """
    if (in_df.source < in_df.target).all():
        return remove_dupe_rows(in_df, wt_attr_name)
    tmp_rcds = [(x, y, z) if x < y else (y, x, z) for x, y, z in in_df.to_dict('split')['data']]
    tmp_df = pd.DataFrame(tmp_rcds, columns=['source', 'target', wt_attr_name])
    return remove_dupe_rows(tmp_df, wt_attr_name)

def load_eda_network(eda_file: str, wt_attr_name: str = 'wt',
                     delimiter: str = r'\s+') -> pd.DataFrame:
    """
    Load network from edge attribute file file (.eda).
    Eda file lists the edges in the following format
    The first line has the string "Weight". The edges are listed with the following
    format:
    source (itype) target = weight

    where source and target are node ids, itype is the interaction type and weight is
    the edge weight.

    Example:
    Weight
    244901_at (pp) 244902_at = 0.192777
    244901_at (pp) 244926_s_at = 0.0817807


    Parameters
    ----------
    eda_file : Path to the input .eda file

    wt_attr_name : weight column name in the data frame returned

    delimiter : seperators between the fields

    Returns
    -------
    pandas DataFrame with three columns: 'source', 'target', wt_attr_name

    """
    tmp_df = pd.read_csv(eda_file, sep=delimiter, usecols=[0, 2, 4], skiprows=[0],
                         names=['source', 'target', wt_attr_name])
    return order_network_rows(tmp_df, wt_attr_name)


def load_tsv_network(tsv_file: str, wt_attr_name: str = 'wt') -> pd.DataFrame:
    """
    Load network with a tsv file.
    The tsv file is expected to have three columns 'source', 'target', wt_attr_name
    If the columns doesn't exist in the data frame, then the first three
    columns are returned.

    Parameters
    ----------
    tsv_file : Path to the input .tsv file

    wt_attr_name : weight column name in the data frame returned

    delimiter : seperators between the fields

    Returns
    -------
    pandas DataFrame with three columns: 'source', 'target', wt_attr_name
    """
    tmp_df = pd.read_csv(tsv_file, sep='\t', header=0)
    if 'source' in tmp_df.columns and 'target' in tmp_df.columns and wt_attr_name in tmp_df.columns:
        return order_network_rows(tmp_df.loc[:, ['source', 'target', wt_attr_name]],
                                  wt_attr_name)
    if 'source' in tmp_df.columns and 'target' in tmp_df.columns and 'wt' in tmp_df.columns:
        tmp_df = tmp_df.loc[:, ['source', 'target', 'wt']]
        tmp_df = tmp_df.rename(columns={'wt': wt_attr_name})
        return order_network_rows(tmp_df, wt_attr_name)
    if len(tmp_df.columns) >= 3:
        cnames = ['source', 'target', wt_attr_name]
        tmp_df = tmp_df.iloc[:, [0, 1, 2]]
        rn_cols = {x : y for x, y in zip(tmp_df.columns, cnames)}
        tmp_df = tmp_df.rename(columns=rn_cols)
        return order_network_rows(tmp_df, wt_attr_name)
    return pd.DataFrame({'source': [], 'target': [], wt_attr_name: []})


def load_adj_network(adj_file: str, wt_attr_name: str = 'wt',
                     delimiter: str = None,
                     comments: List[str] = None) -> pd.DataFrame:
    """
    Load network with adjacency list;
    Adjacency file lists the target lists for each source node:
    source_node target_node1 weight1 ....

    Parameters
    ----------
    adj_file : Path to the input .adj file

    wt_attr_name : weight column name in the data frame returned

    comments : list of comment indicators

    delimiter : seperators between the fields (default: one or more white spaces)

    Returns
    -------
    pandas DataFrame with three columns: 'source', 'target', wt_attr_name
    """
    if comments is None:
        comments = ["#", ">"]
    edge_list: List[Tuple[str, str, float]] = []
    with open(adj_file) as adjfp:
        for line in adjfp:
            if any((line.startswith(cx) for cx in comments)):
                continue
            vlist = line.strip().split(delimiter)
            if len(vlist) >= 3:
                tgt_length = int((len(vlist) - 1) / 2)
                edge_list.extend([(vlist[0], vlist[2*ix+1], float(vlist[2*ix+2]))
                                  for ix in range(tgt_length)])
    edge_list = [(x, y, z) if x < y else (y, x, z) for x, y, z in edge_list]
    return pd.DataFrame(edge_list, columns=['source', 'target', wt_attr_name])



def load_mat_network(mat_file: str, wt_attr_name: str = 'wt',
                     delimiter: str = r'\s+') -> pd.DataFrame:
    """
    Load network with adjacency matrix;
    Adjacency matrix lists the target wts for each source node:
    source_node1 source_node2 source_node_3 ...
    target_node1 weight11 weight12 weight13 ...
    target_node2 weight21 weight22 weight23 ...
    target_node3 weight31 weight32 weight33 ...
    ...
    ...

    Parameters
    ----------
    adj_file : Path to the input .adj file

    wt_attr_name : weight column name in the data frame returned

    delimiter : seperators between the fields (default: one or more white spaces)

    Returns
    -------
    pandas DataFrame with three columns: 'source', 'target', wt_attr_name
    """
    mat_df = pd.read_csv(mat_file, sep=delimiter, index_col=0)
    mat_cnames = mat_df.columns
    mat_size = mat_df.shape[0]
    net_df = pd.DataFrame({
        'source': np.repeat(mat_cnames, mat_size),
        'target': np.tile(mat_cnames, mat_size),
        wt_attr_name: mat_df.values.flatten()})
    return net_df[net_df.source < net_df.target]


def load_reveng_network(net_file: str, wt_attr_name: str = 'wt') -> pd.DataFrame:
    """
    Load network with one of a edge attribute file (.eda),
    tab-seperated file (.tsv), adjacency matrix (.mat),
    adjacency list (.adj)

    Parameters
    ----------
    net_file : Path to the input .eda/.tsv/.mat/.adj file

    wt_attr_name : weight column name in the data frame returned

    delimiter : seperators between the fields

    Returns
    -------
    pandas DataFrame with three columns: 'source', 'target', wt_attr_name
    """
    if net_file.endswith(".eda"):
        return load_eda_network(net_file, wt_attr_name)
    if net_file.endswith(".adj"):
        return load_adj_network(net_file, wt_attr_name)
    if net_file.endswith(".tsv"):
        return load_tsv_network(net_file, wt_attr_name)
    if net_file.endswith(".mat"):
        return load_mat_network(net_file, wt_attr_name)
    return pd.DataFrame({'source': [], 'target': [], wt_attr_name: []})


def map_probes(gs_net: pd.DataFrame, annot_df: pd.DataFrame,
               how_join: str = 'inner') -> pd.DataFrame:
    annot_filter_df = annot_df.loc[annot_df.ID != 'no_match', :]
    gs_net_mapped = gs_net.merge(annot_filter_df, left_on='TF', right_on='ID', how=how_join)
    gs_net_mapped.columns = ['TF', 'TARGET', 'TFID', 'TFPROBE']
    gs_net_mapped = gs_net_mapped.merge(annot_filter_df, left_on='TARGET',
                                        right_on='ID', how=how_join)
    gs_net_mapped.columns = ['TF', 'TARGET', 'TFID', 'TFPROBE', 'TARGETID', 'TARGETPROBE']
    return gs_net_mapped

def map_atid2probes(probe_df: pd.DataFrame, annot_df: pd.DataFrame,
                   how_join: str = 'inner') -> pd.DataFrame:
    annot_filter_df = annot_df.loc[annot_df.ID != 'no_match', :]
    probe_df_mapped = probe_df.merge(annot_filter_df, left_on='ID',
                                     right_on='ID', how=how_join)
    return probe_df_mapped

def map_probes2atid(probe_df: pd.DataFrame, annot_df: pd.DataFrame,
                    how_join: str = 'inner') -> pd.DataFrame:
    annot_filter_df = annot_df.loc[annot_df.ID != 'no_match', :]
    probe_df_mapped = probe_df.merge(annot_filter_df, left_on='PROBE',
                                     right_on='PROBE', how=how_join)
    return probe_df_mapped
