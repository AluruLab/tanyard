import argparse
from typing import Iterable, Set, List, Tuple
import pandas as pd
import data_utils as du

def common_unique_genes(annot_file: str, dataset_files: Iterable[str]):
    #
    annot_df = du.load_annotation(annot_file)
    left_file = dataset_files[0]
    right_file = dataset_files[1]
    lxdf = pd.read_csv(left_file, sep="\t", skiprows=3, comment='D', usecols=[0], header=None, names=["PROBE"])
    rxdf = pd.read_csv(right_file, sep="\t", skiprows=3, comment='D', usecols=[0], header=None, names=["PROBE"])
    lonly_df = pd.DataFrame(data={'PROBE' : list(set(lxdf.PROBE) - set(rxdf.PROBE))})
    ronly_df = pd.DataFrame(data={'PROBE' : list(set(rxdf.PROBE) - set(lxdf.PROBE))})
    lonly_df = du.map_probes2atid(lonly_df, annot_df)
    ronly_df = du.map_probes2atid(ronly_df, annot_df)
    #return (lonly_df.loc[:, ['ID']], ronly_df.loc[:, ['ID']])
    return (lonly_df, ronly_df)

def main(annot_file:str, dataset_files: Iterable[str],
         left_file: str, right_file: str) -> None:
    left_df, right_df = common_unique_genes(annot_file, dataset_files)
    left_df.to_csv(left_file, sep="\t", index=False)
    right_df.to_csv(right_file, sep="\t", index=False)


if __name__ == "__main__":
    PROG_DESC = """
    """
    PARSER = argparse.ArgumentParser(description=PROG_DESC)
    PARSER.add_argument("annotation_file",
                        help="""annotation file
                                (a tab seperated file mapping probe to ids)""")
    PARSER.add_argument("dataset_files", nargs="+",
                        help="""sub network build from a reverse engineering methods
                                (currenlty supported: eda, tsv, adj)""")
    PARSER.add_argument("-l", "--left_file",
                        type=argparse.FileType(mode='w'), required=True,
                        help="output file in tab-seperated format")
    PARSER.add_argument("-r", "--right_file",
                        type=argparse.FileType(mode='w'), required=True,
                        help="output file in tab-seperated format")
    ARGS = PARSER.parse_args()
    main(ARGS.annotation_file, ARGS.dataset_files, ARGS.left_file, ARGS.right_file)
