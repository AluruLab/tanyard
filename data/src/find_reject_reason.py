import pandas as pd
import argparse

tissues = ["leaf", "flower", "rosette", "seedling1wk", "seedling2wk", "wholeplant", "shoot", "seed", "root"]

infile_pattern = "%s/ae-geo-%s.csv.gz"
fullfile_pattern = "%s/full-list.csv"
fnlcsv_pattern = "%s/final-list.csv"
fulllst_pattern = "%s/full-list.txt"
qcrej_pattern = "%s/qc-reject.txt"
mnrej_pattern = "%s/man-reject.txt"
mnacc_pattern = "%s/man-accept.txt"
gsedup_pattern = "%s/gse-super.txt"
finallst_pattern = "%s/final-list.txt"
rjreason_pattern = "%s/%s-reject-reason.csv"

def load_lst_csv(lst_csv_fname):
    return pd.read_csv(lst_csv_fname ,
                       header=None, names=['SeriesId', 'FileId'])

def load_final_csv(tissue_id):
    return load_lst_csv(fnlcsv_pattern % tissue_id)

def final_file_ids(tissue_id):
    fdf = load_final_csv(tissue_id)
    return set(fdf.FileId.astype(str))

def qcr_reject_ids(tissue_id):
    qcrdf = pd.read_csv(qcrej_pattern % tissue_id,
                        header=None, names=['FileId'])
    return set(qcrdf.FileId.astype(str))

def manual_reject_ids(tissue_id):
    manrdf =pd.read_csv(mnrej_pattern % tissue_id,
                        header=None, names=['FileId'])
    return set(manrdf.FileId.astype(str))

def manual_accept_ids(tissue_id):
    mnacdf = pd.read_csv(mnacc_pattern % tissue_id,
                         header=None, names=['FileId'])
    return set(mnacdf.FileId.astype(str))

def gse_duplicate_ids(tissue_id):
    gseddf = pd.read_csv(gsedup_pattern % tissue_id,
                         header=None, names=['FileId'])
    return set(gseddf.FileId.astype(str))

def rejection_file_ids(tissue_id):
    return gse_duplicate_ids(tissue_id) | (
            (qcr_reject_ids(tissue_id) | 
             manual_reject_ids(tissue_id)) - 
            manual_accept_ids(tissue_id)
           ) 
            
def full_list_file_ids(tissue_id):
    fulldf = pd.read_csv(fulllst_pattern % tissue_id ,
                        header=None, names=['FileId'])
    return set(fulldf.FileId.astype(str))

def accepted_file_ids(tissue_id):
    return full_list_file_ids(tissue_id) - rejection_file_ids(tissue_id)

def final_list_ids(tissue_id):
    finaldf = load_final_csv(tissue_id)
    return set(finaldf.FileId.astype(str))


def add_reject_reason(file_id, final_list, qc_reject, manual_reject, manual_accept, gse_duplicates):
    if file_id in final_list:
        return pd.DataFrame({
            'FileId'   : [file_id],
            'Accepted' : ['Y'],
            'QCNotes'  : ['Manual Acceptance' if file_id in manual_accept else 'NA']
        })
    elif file_id in gse_duplicates:
        return pd.DataFrame({
            'FileId'   : [file_id],
            'Accepted' : ['N'],
            'QCNotes'  : ['Duplicate Sample']
        })
    elif file_id in manual_reject:
        return pd.DataFrame({
            'FileId'   : [file_id],
            'Accepted' : ['N'],
            'QCNotes'  : ['Manual Rejection']
        })
    elif file_id in qc_reject:
        return pd.DataFrame({
            'FileId'   : [file_id],
            'Accepted' : ['N'],
            'QCNotes'  : ['Automatic QC Rejection']
        })
    else:
        return pd.DataFrame({
            'FileId'   : [file_id],
            'Accepted' : ['N'],
            'QCNotes'  : ['Corrupt/Missing CEL']
        })

def rejection_reason_df(tissue_id):
    ae_geo_df = pd.read_csv(infile_pattern % (tissue_id, tissue_id),
                            encoding='ISO-8859-1')
    file_id_df = pd.DataFrame({'FileId' : ae_geo_df.FileId})
    file_id_grouped = file_id_df.groupby('FileId')
    final_list = final_list_ids(tissue_id)
    qc_reject = rejection_file_ids(tissue_id)
    manual_reject = manual_reject_ids(tissue_id)
    manual_accept = manual_accept_ids(tissue_id)
    gse_duplicates = gse_duplicate_ids(tissue_id)
    return file_id_grouped.apply(lambda x: add_reject_reason(x.name, final_list, qc_reject, manual_reject, manual_accept, gse_duplicates))


if __name__ == "__main__":
    prog_desc = """Generate Rejection Reason for long list;
    The script should be run from data/tables directory
    """
    parser = argparse.ArgumentParser(description=prog_desc)
    parser.parse_args()
    for tissue_id in tissues:
        rdf = rejection_reason_df(tissue_id)
        rdf.to_csv(rjreason_pattern % (tissue_id, tissue_id), index=False)