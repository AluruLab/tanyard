import urllib
import shutil
import errno
import os
import os.path
import argparse
import requests
import subprocess
import pandas as pd

def last_of(in_txt, sub):
    return in_txt.rfind(sub)

def second_last_of(in_txt, sub):
    first_pos = last_of(in_txt, sub)
    if first_pos > 0:
        return in_txt[:first_pos].rfind(sub)
    return -1

def third_last_of(in_txt, sub):
    second_pos = second_last_of(in_txt, sub)
    if second_pos > 0:
        return in_txt[:second_pos].rfind(sub)
    return -1

def make_dir_path(dir_path):
    try:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
    except OSError as ose:
        if ose.errno != errno.EEXIST:
            raise

def download_http_file(url, local_filename):
    req = requests.get(url, stream=True)
    with open(local_filename, 'wb') as outf:
        shutil.copyfileobj(req.raw, outf)
    return req.status_code

def download_ftp_file(ftp_url, local_filename):
    try:
        urllib.request.urlretrieve(ftp_url, filename=local_filename)
    except urllib.error.URLError as uerr:
        if 'No such file or directory' in  uerr.reason:
            return 550
        print(ftp_url)
        raise
    return 200

def download_file(url, local_filename):
    if url.startswith('ftp'):
        return download_ftp_file(url, local_filename)
    return download_http_file(url, local_filename)

def download_file_wget(url, local_filename):
    rx = subprocess.run(["wget", url, "-nv", "-O", local_filename], stdout=subprocess.PIPE)
    if rx == 0:
        return 200
    else:
        return 550

def starts_with(pdf, attr_id, pfx_txt):
    return pd.Series([True if x.startswith(pfx_txt) else False for x in pdf[attr_id]],
                     index=pdf.index)

def cel_count(clist):
    return len([x for x in clist if x.lower().endswith('.cel') or x.lower().endswith('cel.gz')])

def empty_file(pdf):
    return pd.Series([True if len(str(x)) < 5 else False for x in pdf.SampleFile])

def has_file(pdf):
    return pd.Series([True if len(str(x)) >= 5 else False for x in pdf.SampleFile],
                     index=pdf.index)

def has_semicolon(pdf):
    return pd.Series([True if (';' in str(x)) else False for x in pdf.SampleFile],
                     index=pdf.index)

def has_no_semicolon(pdf):
    return has_semicolon(pdf) == False

def ends_with_text(pdf):
    return pd.Series([True if str(x).endswith('txt') else False for x in pdf.SampleFile], 
                     index=pdf.index)

def not_ends_with_text(pdf):
    return ends_with_text(pdf) == False

def contains_cel(pdf):
    return pd.Series([True if 'cel' in x.lower() else False for x in pdf.SampleFile],
                     index=pdf.index)

def valid_cel(pdf):
    return pd.Series(
        [True if x.lower().endswith('cel') or x.lower().endswith('cel.gz') else False for x in pdf.SampleFile],
        index=pdf.index
    )

def cel_count_eq_n(pdf, n):
    return pd.Series(
        [True if cel_count(x.split(';')) == n else False for x in pdf.SampleFile],
        index=pdf.index)

def cel_count_gt_n(df, n):
    return pd.Series(
        [True if cel_count(x.split(';')) > n else False for x in df.SampleFile],
        index=df.index)

