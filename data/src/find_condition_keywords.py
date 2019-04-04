import pandas as pd
import numpy as np
import argparse

tissues = ["leaf", "flower", "rosette", "seedling1wk", "seedling2wk", "wholeplant", "shoot", "seed", "root"]

condition_keywords = {
    'stress'   : ['cold', 'heat', 'osmotic','salt','drought',
                  'oxidatine', 'uv-b', 'wounding', 'pathogen', 'elicitor' ],
    'hormone'  : ['indoce acetic acid', 'gibberellic acid', 'zeatin',
                  'jasmonic acid', 'abscisic acid', 'cytokinin',
                  'gibberellin', 'auxin', 'ethylene', 'brassinosteroids',
                  'brassinolides', 'salicylic acid'],
    'chemical' : ['phosphorus', 'phosphate', 'sulfur', 'sulphur',
                  'nitrogen', 'calcicum', 'iron', 'potassium',
                  'manganese', 'magnesium', 'zinc', 'copper']
}


infile_pattern = "%s/ae-geo-%s.csv.gz"

def get_accepted(tissue_id):
    ae_geo_df = pd.read_csv(infile_pattern % (tissue_id, tissue_id),
                        encoding='ISO-8859-1')
    #rj_reason_df = pd.read_csv(rjreason_pattern % (tissue_id, tissue_id))
    ae_geo_df = ae_geo_df.loc[ae_geo_df.Accepted == 'Y', :]
    #print(tissue_id, ae_geo_df.shape[0])
    return ae_geo_df


def column_contains(in_df, col_name, contain_string):
    return in_df[col_name].str.contains(contain_string,
                        regex=False, case=False)

def any_column_contains(in_df, col_names, contain_string):
    return pd.DataFrame({
            str(cx) : column_contains(in_df, cx, contain_string) 
                for cx in col_names},
            index=in_df.index).any(axis='columns')

def any_column_contains_map(in_df, col_names, cx_keyword):
    return any_column_contains(in_df, col_names, cx_keyword).apply(
                lambda x: str(cx_keyword) if x else np.nan)

def columns_contain_keywords(in_df, keywords, col_names):
    return pd.DataFrame(data = { 
            str(kx) : any_column_contains_map(in_df, col_names, kx) for kx in keywords
        },
        index=in_df.index)

def join_keywords(col_vals): 
    jx = ",".join(str(y) for y in col_vals if str(y) !='nan')
    return jx if len(jx) > 0 else np.nan

def find_keywords_in(in_df, keywords, col_names = ['SeriesText', 'SampleText']):
    contains_df = columns_contain_keywords(in_df, keywords, col_names)
    return contains_df.agg(lambda x: join_keywords(x), axis=1)

def find_condition_keywords(in_df, conditions):
    rw_dict = {cx : find_keywords_in(in_df, condition_keywords[cx])
                       for cx in conditions}
    rw_dict['FileId'] = in_df['FileId']
    return pd.DataFrame(rw_dict, index=in_df.index)

if __name__ == "__main__":
    prog_desc = """Create a master sheet of accepted datasests;
    Find keywords for conditions;
    The script should be run from data/tables directory
    """
    parser = argparse.ArgumentParser(description=prog_desc)
    parser.add_argument("-o", "out_file", default="Dataset-MasterSheet-Accepted.csv")
    args = parser.parse_args()
    master_df = pd.concat([get_accepted(tissue_id) for tissue_id in tissues], sort=False)
    cond_df = find_condition_keywords(master_df, ['stress', 'hormone', 'chemical'])
    master_df.merge(cond_df).to_csv(args.out_file, index=False)
