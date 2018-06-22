import pandas as pd
import requests
import shutil
import os
import errno
import os.path
import time

def make_dir_path(dir_path):
    try:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

def download_file(url, local_filename):
    r = requests.get(url, stream=True)
    with open(local_filename, 'wb') as f:
        shutil.copyfileobj(r.raw, f)
    return r.status_code

        
def download_arrayexpress_metadata(ae_file, meta_dir,
                                   download_status_file):
    aedf = pd.read_csv(ae_file, sep='\t', encoding = "ISO-8859-1")
    aexp_url = 'https://www.ebi.ac.uk/arrayexpress/files'
    access_list = []; idf_rcodes = []; sdrf_rcodes = []
    aejson_url = 'https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}'
    samples_url = 'https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}/samples'
    proto_url = 'https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}/protocols'
    files_url = 'https://www.ebi.ac.uk/arrayexpress/json/v3/files?accession={}'
    exprts_rcodes = []; samples_rcodes = []; proto_rcodes = []; pfiles_rcodes = []
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
        rx = download_file(idf_url, idf_fname)
        print(rx)
        idf_rcodes.append(rx)
        print(sdrf_url, idf_fname, end='')
        rx = download_file(sdrf_url, sdrf_fname)
        print(rx)
        sdrf_rcodes.append(rx)
        aejurl = aejson_url.format(axn)
        aejfn = os.path.join(axn_dir, axn + ".expts.json")
        rx = download_file(aejurl, aejfn)
        print(aejurl, aejfn, rx)
        exprts_rcodes.append(rx)
        samurl = samples_url.format(axn)
        samfn = os.path.join(axn_dir, axn + ".samples.json")
        rx = download_file(samurl, samfn)
        print(samurl, samfn, rx)
        samples_rcodes.append(rx)
        proturl = proto_url.format(axn)
        protfn = os.path.join(axn_dir, axn + ".protocols.json")
        rx = download_file(proturl, protfn)
        print(proturl, protfn, rx)
        proto_rcodes.append(rx)
        filurl = files_url.format(axn)
        filfn = os.path.join(axn_dir, axn + ".files.json")
        rx = download_file(filurl, filfn)
        print(filurl, filfn, rx)
        pfiles_rcodes.append(rx)
        ncount += 1
        if ncount % 50 == 0:
            print("WAITING 5 SEC... ",end='')
            time.sleep(5)
            print(" ...DONE")
    rdf = pd.DataFrame(data = {'Accession' : access_list, 
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
