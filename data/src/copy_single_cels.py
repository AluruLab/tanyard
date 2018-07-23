import argparse
import os.path
import pathlib
import shutil
import pandas as pd

def main(in_file, out_data_dir, out_file):
    dx = pd.read_csv(in_file, encoding = "ISO-8859-1")
    rcode_lst = []
    for rid, dfrow in dx.iterrows():
        local_cel_file = dfrow.PDF
        if pathlib.Path(local_cel_file).is_file():
            rx = shutil.copy(local_cel_file, out_data_dir)
            rcode_lst.append(rx)
        else:
            rcode_lst.append(-1)
    rdx = pd.DataFrame(data={
        'RowId' : dx['RowId'],'CEL' : dx['CEL'],
        'PDF' : dx['PDF'], 'ReturnCode': rcode_lst })
    rdx.to_csv(out_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("in_file")
    parser.add_argument("data_dir")
    parser.add_argument("out_file")
    args = parser.parse_args()
    main(args.in_file, args.data_dir, args.out_file)
