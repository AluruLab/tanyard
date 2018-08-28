import argparse
import os.path
import pathlib
import pandas as pd
import download_utils as DU


def update_nasc_status(celdf, data_dir):
    ae_prefix_url = 'https://www.ebi.ac.uk/arrayexpress/files'
    rowid_lst = []
    series_lst = []
    sample_lst = []
    url_lst = []
    file_lst = []
    rcode_lst = []
    for rid, dfrow in celdf.iterrows():
        if dfrow.SeriesId.startswith("NASC-"):
            ext_name = ".CEL"
            cel_fname = dfrow.SeriesId + "_" + str(dfrow.SampleId) + ext_name
            local_cel_file = os.path.join(series_dir, cel_fname)
            if pathlib.Path(local_cel_file).is_file():
                rcode = 0
            else:
                rcode = 550
            rowid_lst.append(rid)
            series_lst.append(dfrow.SeriesId)
            sample_lst.append(dfrow.SampleId)
            url_lst.append("NA")
            file_lst.append(local_cel_file)
            rcode_lst.append(rcode)
            print(rid, dfrow.SeriesId, dfrow.SampleId, "NA", local_cel_file, rcode)
    return pd.DataFrame(data={
        'RowId' : rowid_lst,'SeriesId' : series_lst,
        'SampleId' : sample_lst, 'SampleFile': url_lst,
        'SampleCEL' : file_lst, 'ReturnCode': rcode_lst
    })


def main(in_file, data_dir, out_status_file):
    dx = pd.read_csv(in_file, encoding = "ISO-8859-1")
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
    retdf = update_nasc_status(oneceldf, data_dir)
    retdf.to_csv(out_status_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("in_file")
    parser.add_argument("data_dir")
    parser.add_argument("out_status_file")
    args = parser.parse_args()
    main(args.in_file, args.data_dir, args.out_status_file)
