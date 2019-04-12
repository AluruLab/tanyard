import time
import shutil
import os
import os.path
import errno
import requests
import pandas as pd

def make_dir_path(dir_path):
    try:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
    except OSError as osex:
        if osex.errno != errno.EEXIST:
            raise

def download_file(url, local_filename):
    rqx = requests.get(url, stream=True)
    with open(local_filename, 'wb') as fptr:
        shutil.copyfileobj(rqx.raw, fptr)
    return rqx.status_code


def download_arrayexpress_metadata(ae_file, meta_dir,
                                   download_status_file):
    aedf = pd.read_csv(ae_file, sep='\t', encoding="ISO-8859-1")
    aexp_url = 'https://www.ebi.ac.uk/arrayexpress/files'
    aejson_url = 'https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}'
    samples_url = 'https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}/samples'
    proto_url = 'https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}/protocols'
    files_url = 'https://www.ebi.ac.uk/arrayexpress/json/v3/files?accession={}'
    access_list = []
    idf_rcodes = []
    sdrf_rcodes = []
    exprts_rcodes = []
    samples_rcodes = []
    proto_rcodes = []
    pfiles_rcodes = []
    ncount = 0
    for axn in aedf["Accession"]:
        # create dir and local files
        axn_dir = os.path.join(meta_dir, axn)
        make_dir_path(axn_dir)
        idf_base = axn + ".idf.txt"
        sdrf_base = axn + ".sdrf.txt"
        idf_fname = os.path.join(axn_dir, axn + ".idf.txt")
        sdrf_fname = os.path.join(axn_dir, axn + ".sdrf.txt")
        # build URLs
        idf_url = "/".join([aexp_url, axn, idf_base])
        sdrf_url = "/".join([aexp_url, axn, sdrf_base])
        access_list.append(axn)
        print(idf_url, idf_fname, end='')
        r_code = download_file(idf_url, idf_fname)
        print(r_code)
        idf_rcodes.append(r_code)
        print(sdrf_url, idf_fname, end='')
        r_code = download_file(sdrf_url, sdrf_fname)
        print(r_code)
        sdrf_rcodes.append(r_code)
        aejurl = aejson_url.format(axn)
        aejfn = os.path.join(axn_dir, axn + ".expts.json")
        r_code = download_file(aejurl, aejfn)
        print(aejurl, aejfn, r_code)
        exprts_rcodes.append(r_code)
        samurl = samples_url.format(axn)
        samfn = os.path.join(axn_dir, axn + ".samples.json")
        r_code = download_file(samurl, samfn)
        print(samurl, samfn, r_code)
        samples_rcodes.append(r_code)
        proturl = proto_url.format(axn)
        protfn = os.path.join(axn_dir, axn + ".protocols.json")
        r_code = download_file(proturl, protfn)
        print(proturl, protfn, r_code)
        proto_rcodes.append(r_code)
        filurl = files_url.format(axn)
        filfn = os.path.join(axn_dir, axn + ".files.json")
        r_code = download_file(filurl, filfn)
        print(filurl, filfn, r_code)
        pfiles_rcodes.append(r_code)
        ncount += 1
        if ncount % 50 == 0:
            print("WAITING 5 SEC... ", end='')
            time.sleep(5)
            print(" ...DONE")
    rdf = pd.DataFrame(data={'Accession' : access_list,
                             'IDF_RETURN': idf_rcodes,
                             'SDRF_RETURN': sdrf_rcodes,
                             'EXPT_RETURN': exprts_rcodes,
                             'SAMPLE_RETURN': samples_rcodes,
                             'PROTO_RETURN': proto_rcodes,
                             'FILES_RETURN': pfiles_rcodes})
    rdf.to_csv(download_status_file)


def main():
    meta_dir = '/home/sriram/work/tanyard/data/meta/AE'
    ae_file = "../tables/AE-ATH1.tsv"
    download_status_file = "ae_download_stats.csv"
    download_arrayexpress_metadata(ae_file, meta_dir,
                                   download_status_file)

if __name__ == '__main__':
    main()
