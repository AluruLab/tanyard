import argparse
import os.path
import pathlib
import shutil
import pandas as pd
import download_utils as DU

def main(in_file, out_data_dir, out_file, file_fmt):
    dx = pd.read_csv(in_file, encoding = "ISO-8859-1")
    rcode_lst = []
    DU.make_dir_path(out_data_dir)
    for rid, dfrow in dx.iterrows():
        if file_fmt.lower() == "pdf":
           local_file = dfrow.PDF
        elif file_fmt.lower() == "png":
           local_file = dfrow.PNG
        else:
           local_file = dfrow.CEL
        if pathlib.Path(local_file).is_file():
            rx = shutil.copy(local_file, out_data_dir)
            rcode_lst.append(rx)
        else:
            rcode_lst.append(-1)
    if file_fmt.lower() == "png":
        rdx = pd.DataFrame(data={
            'RowId' : dx['RowId'],'CEL' : dx['CEL'],
            'PNG' : dx['PNG'], 'ReturnCode': rcode_lst })
    else:
        rdx = pd.DataFrame(data={
            'RowId' : dx['RowId'],'CEL' : dx['CEL'],
            'PDF' : dx['PDF'], 'ReturnCode': rcode_lst })
    rdx.to_csv(out_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("in_file")
    parser.add_argument("data_dir")
    parser.add_argument("out_file")
    parser.add_argument("file_fmt")
    args = parser.parse_args()
    main(args.in_file, args.data_dir, args.out_file, args.file_fmt)
