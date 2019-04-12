import argparse
import pandas as pd

TISSUES = ["leaf", "flower", "rosette", "seedling1wk", "seedling2wk",
           "wholeplant", "shoot", "seed", "root"]

TISSUE_FILE_PATTERN = "%s/ae-geo-%s.csv.gz"
TISSUE_FULLFILE_PATTERN = "%s/full-list.csv"
TISSUE_FNLCSV_PATTERN = "%s/final-list.csv"
TISSUE_FULLLST_PATTERN = "%s/full-list.txt"
TISSUE_QCR_PATTERN = "%s/qc-reject.txt"
TISSUE_MNR_PATTERN = "%s/man-reject.txt"
TISSUE_MNAC_PATTERN = "%s/man-accept.txt"
TISSUE_GSEDP_PATTERN = "%s/gse-super.txt"
TISSUE_FNLLST_PATTERN = "%s/final-list.txt"
TISSUE_RJR_PATTERN = "%s/%s-reject-reason.csv"

def load_lst_csv(lst_csv_fname):
    return pd.read_csv(lst_csv_fname,
                       header=None, names=['SeriesId', 'FileId'])

def load_final_csv(tissue_id):
    return load_lst_csv(TISSUE_FNLCSV_PATTERN % tissue_id)

def final_file_ids(tissue_id):
    fdf = load_final_csv(tissue_id)
    return set(fdf.FileId.astype(str))

def qcr_reject_ids(tissue_id):
    qcrdf = pd.read_csv(TISSUE_QCR_PATTERN % tissue_id,
                        header=None, names=['FileId'])
    return set(qcrdf.FileId.astype(str))

def manual_reject_ids(tissue_id):
    manrdf = pd.read_csv(TISSUE_MNR_PATTERN % tissue_id,
                         header=None, names=['FileId'])
    return set(manrdf.FileId.astype(str))

def manual_accept_ids(tissue_id):
    mnacdf = pd.read_csv(TISSUE_MNAC_PATTERN % tissue_id,
                         header=None, names=['FileId'])
    return set(mnacdf.FileId.astype(str))

def gse_duplicate_ids(tissue_id):
    gseddf = pd.read_csv(TISSUE_GSEDP_PATTERN % tissue_id,
                         header=None, names=['FileId'])
    return set(gseddf.FileId.astype(str))

def rejection_file_ids(tissue_id):
    return gse_duplicate_ids(tissue_id) | (
        (qcr_reject_ids(tissue_id) |
         manual_reject_ids(tissue_id)) -
        manual_accept_ids(tissue_id)
        )

def full_list_file_ids(tissue_id):
    fulldf = pd.read_csv(TISSUE_FULLLST_PATTERN % tissue_id,
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
    ae_geo_df = pd.read_csv(TISSUE_FILE_PATTERN % (tissue_id, tissue_id),
                            encoding='ISO-8859-1')
    file_id_df = pd.DataFrame({'FileId' : ae_geo_df.FileId})
    file_id_grouped = file_id_df.groupby('FileId')
    final_list = final_list_ids(tissue_id)
    qc_reject = rejection_file_ids(tissue_id)
    manual_reject = manual_reject_ids(tissue_id)
    manual_accept = manual_accept_ids(tissue_id)
    gse_duplicates = gse_duplicate_ids(tissue_id)
    return file_id_grouped.apply(
        lambda x: add_reject_reason(x.name, final_list, qc_reject,
                                    manual_reject, manual_accept,
                                    gse_duplicates))

def main():
    for tissue_id in TISSUES:
        rdf = rejection_reason_df(tissue_id)
        rdf.to_csv(TISSUE_RJR_PATTERN % (tissue_id, tissue_id), index=False)

if __name__ == "__main__":
    PROG_DESC = """Generate Rejection Reason for long list;
    The script should be run from data/tables directory
    """
    argparse.ArgumentParser(description=PROG_DESC).parse_args()
    main()
