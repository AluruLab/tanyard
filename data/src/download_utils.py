"""
Utilites for downloading CEL files
"""
import urllib
import shutil
import errno
import os
import os.path
import subprocess
import requests
import pandas as pd
import numpy as np

def last_of(in_txt: str, sub: str):
    "Identify the last occurence of sub in in_txt"
    return in_txt.rfind(sub)

def second_last_of(in_txt: str, sub: str):
    "Identify the penultimate occurence of sub in in_txt"
    first_pos = last_of(in_txt, sub)
    if first_pos > 0:
        return in_txt[:first_pos].rfind(sub)
    return -1

def third_last_of(in_txt: str, sub: str):
    "Identify the third from the last occurence of sub in in_txt"
    second_pos = second_last_of(in_txt, sub)
    if second_pos > 0:
        return in_txt[:second_pos].rfind(sub)
    return -1

def make_dir_path(dir_path: str):
    "Make all the directories in a given path"
    try:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
    except OSError as ose:
        if ose.errno != errno.EEXIST:
            raise


def copy_file(source_file: str, dest_file: str):
    "Copy source_file to dest_file"
    dest_dir = os.path.dirname(dest_file)
    make_dir_path(dest_dir)
    return shutil.copy(source_file, dest_file)


def download_http_file(url: str, local_filename: str):
    """Download input http url file to local_filename;
       Return status of http request"""
    req = requests.get(url, stream=True)
    with open(local_filename, 'wb') as outf:
        shutil.copyfileobj(req.raw, outf)
    return req.status_code

def download_ftp_file(ftp_url: str, local_filename: str):
    """Download input ftp url file to local_filename;
       Return status of 550 if file not found, 200 otherwise"""
    try:
        urllib.request.urlretrieve(ftp_url, filename=local_filename)
    except urllib.error.URLError as uerr:
        if 'No such file or directory' in  uerr.reason:
            return 550
        print(ftp_url)
        raise
    return 200

def download_file(url: str, local_filename: str):
    """Download input http/ftp url file to local_filename ;
    Uses the requests for http; urllib for ftp"""
    if url.startswith('ftp'):
        return download_ftp_file(url, local_filename)
    return download_http_file(url, local_filename)

def download_file_wget(url: str, local_filename: str):
    """Download input htp/ftp url file to local_filename using wget.
       Returns 200 with sucessful; 550 otherwise"""
    r_code = subprocess.run(["wget", url, "-nv", "-O", local_filename],
                            stdout=subprocess.PIPE)
    if r_code == 0:
        return 200
    else:
        return 550

def starts_with(pdf, attr_id, pfx_txt):
    """Download input htp/ftp url file to local_filename using wget"""
    return pd.Series([x.startswith(pfx_txt) for x in pdf[attr_id]],
                     index=pdf.index)

def cel_count(clist):
    """Counts the number of strings that ends with cel or cel.gz"""
    return len([x for x in clist if x.lower().endswith('.cel') or x.lower().endswith('cel.gz')])

def empty_file(pdf: pd.DataFrame):
    """Check if the SampleFile is empty
    Returns a booolean column for the SampleFile column in pdf"""
    return pd.Series([len(str(x)) < 5 for x in pdf.SampleFile])

def has_file(pdf):
    """Check if the SampleFile is not an empty file
    Return a boleean column for the SampleFile column in pdf data frame"""
    return pd.Series([len(str(x)) >= 5 for x in pdf.SampleFile],
                     index=pdf.index)

def has_semicolon(pdf):
    """Check if the SampleFile has ';' (implies multiple CEL files)
    Return a boleean column for the SampleFile column in pdf data frame"""
    return pd.Series([(';' in str(x)) for x in pdf.SampleFile],
                     index=pdf.index)

def has_no_semicolon(pdf):
    """Check if the SampleFile has no ';' (implies multiple CEL files)
    Return a boleean column for the SampleFile column in pdf data frame"""
    return np.logical_not(has_semicolon(pdf))

def ends_with_text(pdf):
    return pd.Series([str(x).endswith('txt') for x in pdf.SampleFile],
                     index=pdf.index)

def not_ends_with_text(pdf):
    return np.logical_not(ends_with_text(pdf))

def contains_cel(pdf):
    return pd.Series(['cel' in x.lower() for x in pdf.SampleFile],
                     index=pdf.index)

def valid_cel(pdf: pd.DataFrame):
    return pd.Series(
        [x.lower().endswith('cel') or x.lower().endswith('cel.gz') for x in pdf.SampleFile],
        index=pdf.index
    )

def cel_count_eq_n(pdf, n_files):
    return pd.Series(
        [cel_count(x.split(';')) == n_files for x in pdf.SampleFile],
        index=pdf.index)

def cel_count_gt_n(pdf, n_files):
    return pd.Series(
        [cel_count(x.split(';')) > n_files for x in pdf.SampleFile],
        index=pdf.index)
