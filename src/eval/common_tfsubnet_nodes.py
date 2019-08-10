from typing import Iterable, Set, Any, List, Tuple
import argparse
import pandas as pd
from data_utils import load_annotation, load_reveng_network, load_gsnetwork, map_probes_cols


def common_network(network_files: Iterable[str]):
    common_tfs: Set[Tuple[str, str]] = set([])
    distinct_tfs: List[Set[Tuple[str, str]]] = []
    for net_file in network_files:
        rv_net: pd.DataFrame = pd.read_csv(net_file, sep='\t')
        src_rcds = rv_net.loc[rv_net.source_ind, ['source_id', 'source_alias']].to_records(index=False)
        tgt_rcds = rv_net.loc[rv_net.target_ind, ['target_id', 'target_alias']].to_records(index=False)
        tf_rcds = list(src_rcds) + list(tgt_rcds)
        print(len(tf_rcds))
        if len(common_tfs) == 0:
            common_tfs = set([(x, y) for x, y in tf_rcds])
        else:
            common_tfs = common_tfs & set([(x, y) for x, y in tf_rcds])
    rmax = len(common_tfs)
    for net_file in network_files:
        rv_net: pd.DataFrame = pd.read_csv(net_file, sep='\t')
        src_rcds = rv_net.loc[rv_net.source_ind, ['source_id', 'source_alias']].to_records(index=False)
        tgt_rcds = rv_net.loc[rv_net.target_ind, ['target_id', 'target_alias']].to_records(index=False)
        tf_rcds = list(src_rcds) + list(tgt_rcds)
        unq_rcs = set([(x, y) for x, y in tf_rcds]) - common_tfs
        if len(unq_rcs) > rmax:
            rmax = len(unq_rcs)
        distinct_tfs.append(unq_rcs)
    cdf = {}
    cdf['COMMON_IDS'] = [x for x, _ in common_tfs] + ["" for _ in range(rmax-len(common_tfs))]
    cdf['COMMON_ALIAS'] = [y for _, y in common_tfs] + ["" for _ in range(rmax-len(common_tfs))]
    for ix, net_file in enumerate(network_files):
        utfs = distinct_tfs[ix]
        cdf['IDS_' + net_file] = [x for x, _ in utfs] + ["" for _ in range(rmax-len(utfs))]
        cdf['ALIAS_' + net_file] = [y for _, y in utfs] + ["" for _ in range(rmax-len(utfs))]
    return pd.DataFrame(cdf)


def main(network_files: Iterable[str], out_file: str) -> None:
    out_df = common_network(network_files)
    out_df.to_csv(out_file, sep="\t", index=False)


if __name__ == "__main__":
    PROG_DESC = """
    Finds a common gold network as the intersection of input networks.
    Network intersection is computed as the intersection of edges of the input networks.
    """
    PARSER = argparse.ArgumentParser(description=PROG_DESC)
    PARSER.add_argument("network_files", nargs="+",
                        help="""sub network build from a reverse engineering methods
                                (currenlty supported: eda, tsv, adj)""")
    PARSER.add_argument("-o", "--out_file",
                        type=argparse.FileType(mode='w'), required=True,
                        help="output file in tab-seperated format")
    ARGS = PARSER.parse_args()
    main(ARGS.network_files, ARGS.out_file)
