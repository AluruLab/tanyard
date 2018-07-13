import argparse
import os.path
import pathlib
import pandas as pd
import download_utils as DU

def download_single_cels(celdf, data_dir):
    ae_prefix_url = 'https://www.ebi.ac.uk/arrayexpress/files'
    rowid_lst = []
    series_lst = []
    sample_lst = []
    url_lst = []
    file_lst = []
    rcode_lst = []
    for rid, dfrow in celdf.iterrows():
        ext_name = ".cel" if dfrow.SampleFile.lower().endswith(".cel") else ".cel.gz"
        series_dir = os.path.join(data_dir, dfrow.SeriesId)
        cel_fname = dfrow.SeriesId + "_" + str(dfrow.SampleId) + ext_name
        local_cel_file = os.path.join(series_dir, cel_fname)
        if dfrow.SampleFile.startswith('ftp') and dfrow.SeriesId.startswith('E-'):
            fxurl = dfrow.SampleFile
            file_url = ae_prefix_url + fxurl[DU.third_last_of(fxurl, '/'):]
        else:
            file_url = dfrow.SampleFile
        if pathlib.Path(local_cel_file).is_file():
            rcode = 0
        else:
            DU.make_dir_path(series_dir)
            rcode = DU.download_file(file_url, local_cel_file)
        rowid_lst.append(rid)
        series_lst.append(dfrow.SeriesId)
        sample_lst.append(dfrow.SampleId)
        url_lst.append(file_url)
        file_lst.append(local_cel_file)
        rcode_lst.append(rcode)
        print(rid, dfrow.SeriesId, dfrow.SampleId, file_url, local_cel_file, rcode)
    return pd.DataFrame(data={
        'RowId' : rowid_lst,'SeriesId' : series_lst,
        'SampleId' : sample_lst, 'SampleFile': url_lst,
        'SampleCEL' : file_lst, 'ReturnCode': rcode_lst
    })

def print_classification(dx):
    hsdf = dx.loc[DU.has_file, : ]
    vxdf = hsdf.loc[DU.has_no_semicolon, : ]
    txdf = vxdf.loc[DU.ends_with_text, : ]
    ntxdf = vxdf.loc[DU.not_ends_with_text, : ]
    celtxdf = ntxdf.loc[DU.not_ends_with_text, : ]
    cel1txdf = celtxdf.loc[DU.valid_cel, : ]
    hsmtxdf = hsdf.loc[DU.has_semicolon, : ]
    hsmlztxdf = hsmtxdf.loc[lambda df: DU.cel_count_eq_n(df, 0), :]
    hsmlotxdf = hsmtxdf.loc[lambda df: DU.cel_count_eq_n(df, 1), :]
    hsmlttxdf = hsmtxdf.loc[lambda df: DU.cel_count_eq_n(df, 2), :]
    hsmlhtxdf = hsmtxdf.loc[lambda df: DU.cel_count_eq_n(df, 3), :]
    hsmlmtxdf = hsmtxdf.loc[lambda df: DU.cel_count_gt_n(df, 3), :]
    return """{}
   |- {} (NO FILE)
   |- {} (HAS FILE)
        |- {} (NO SEMICOLON)
        |   |- {} (TXT)
        |   |- {} (NO TXT)
        |       |- {} (CONTAINS CEL)
        |           |- {} (ENDS WITH CEL/CEL.GZ)
        |- {} (HAS SEMICOLON)
            |- {} (HAS 0 CEL)
            |- {} (HAS 1 CEL)
            |- {} (HAS 2 CEL)
            |- {} (HAS 3 CEL)
            |- {} (HAS >3 CEL)""".format(
            len(dx), len(dx.loc[DU.empty_file, : ]),
            len(hsdf), len(vxdf), len(txdf), len(ntxdf), 
            len(celtxdf), len(cel1txdf), len(hsmtxdf),
            len(hsmlztxdf), len(hsmlotxdf), len(hsmlttxdf),
            len(hsmlhtxdf), len(hsmlmtxdf))


def main(in_file, data_dir, out_status_file):
    dx = pd.read_csv(in_file)
    hsdf = dx.loc[DU.has_file, : ]
    vxdf = hsdf.loc[DU.has_no_semicolon, : ]
    ntxdf = vxdf.loc[DU.not_ends_with_text, : ]
    celtxdf = ntxdf.loc[DU.not_ends_with_text, : ]
    #celtxaedf = celtxdf.loc[lambda df: DU.starts_with(df, 'SeriesId', 'E-'), : ]
    #celtxgedf = celtxdf.loc[lambda df: starts_with(df, 'SeriesId', 'GSE'), : ]
    #celtxae_ftpdf = celtxaedf.loc[lambda df: starts_with(df, 'SampleFile', 'ftp'), : ]
    #celtxae_httpdf = celtxaedf.loc[lambda df: starts_with(df, 'SampleFile', 'http'), : ]
    hsmtxdf = hsdf.loc[DU.has_semicolon, : ]
    #hsmlztxdf = hsmtxdf.loc[lambda df: DU.cel_count_eq_n(df, 0), :]
    hsmlotxdf = hsmtxdf.loc[lambda df: DU.cel_count_eq_n(df, 1), :]
    oneceldf = pd.concat([celtxdf, hsmlotxdf])
    print(len(oneceldf))
    retdf = download_single_cels(oneceldf, data_dir)
    retdf.to_csv(out_status_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("in_file")
    parser.add_argument("data_dir")
    parser.add_argument("out_status_file")
    args = parser.parse_args()
    main(args.in_file, args.data_dir, args.out_status_file)
