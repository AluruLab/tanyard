import argparse
from typing import List, Dict
import pandas as pd

COMMON_DF_NAME = 'COMMON'

def read_exp_probes(dataset_file: str, col_names: List[str] = None) -> pd.DataFrame:
    if not col_names:
        col_names = ['PROBE', 'ID']
    return pd.read_csv(dataset_file, sep="\t", names=col_names,
                       usecols=[0, 1], skiprows=3)


def common_genes(dataset_names: List[str],
                 dataset_files: List[str]) -> Dict[str, pd.DataFrame]:
    genes_list_dct = {}
    cmn_df = pd.DataFrame(columns=['PROBE', 'ID'])
    for dname in dataset_names:
        genes_list_dct[dname] = pd.DataFrame()
    for _, file_name in zip(dataset_names, dataset_files):
        tcg_df = read_exp_probes(file_name, ['PROBE', 'ID'])
        if cmn_df.empty:
            cmn_df = pd.DataFrame(tcg_df)
        else:
            cmn_df = cmn_df.merge(tcg_df)
    for dname, file_name in zip(dataset_names, dataset_files):
        tcg_df = read_exp_probes(file_name, ['PROBE', 'ID'])
        tcg_df = cmn_df.merge(cmn_df, indicator=True, how='outer')
        genes_list_dct[dname] = tcg_df[tcg_df['_merged'] == 'right_only']
    genes_list_dct[COMMON_DF_NAME] = cmn_df
    return genes_list_dct


def main(dataset_names: str, dataset_files: List[str], out_file: str) -> bool:
    if dataset_names:
        dataset_names = dataset_names.split(",")
        if len(dataset_names) == len(dataset_names):
            combine_df_dct = common_genes(dataset_names, dataset_files)
        else:
            print("Length of dataset names should be equal to length of network files")
            return False
    else:
        ndatasets = len(dataset_files)
        combine_df_dct = common_genes(['data_' + str(x) for x in range(ndatasets)],
                                      dataset_files)
    with pd.ExcelWriter(out_file) as writer:
        for dname, cdf in combine_df_dct.items():
            cdf.to_excel(writer, sheet_name=dname)
    return True

if __name__ == "__main__":
    PROG_DESC = """
    Common and Unique Probes for Tissue/Condition"""
    PARSER = argparse.ArgumentParser(description=PROG_DESC)
    PARSER.add_argument("dataset_files", nargs="+",
                        help="""exp files corresponding to the datasets""")
    PARSER.add_argument("-n", "--dataset_names", type=str,
                        help="""comma seperated names of the datasets;
                                should have as many names as the number of datasets""")
    PARSER.add_argument("-o", "--out_file",
                        type=argparse.FileType(mode='w'), required=True,
                        help="output file in tab-seperated format")
    ARGS = PARSER.parse_args()
    print("""
       ARG : dataset_files : %s
       ARG : dataset_names : %s
       ARG : out_file : %s """ %
          (str(ARGS.dataset_files), str(ARGS.dataset_names),
           str(ARGS.out_file)))
    if not main(ARGS.dataset_names,
                ARGS.dataset_files, ARGS.out_file):
        PARSER.print_usage()
