import argparse
import pathlib
import shutil
import pandas as pd
import download_utils as DU

def main(in_file, out_data_dir, out_file, file_fmt, root_dir):
    if out_data_dir[-1] != '/':
        out_data_dir = out_data_dir + "/"
    cel_df = pd.read_csv(in_file, encoding="ISO-8859-1")
    rcode_lst = []
    DU.make_dir_path(out_data_dir)
    for _, dfrow in cel_df.iterrows():
        if file_fmt.lower() == "pdf":
            local_file = dfrow.PDF
        elif file_fmt.lower() == "png":
            local_file = dfrow.PNG
        else:
            local_file = dfrow.CEL
        if pathlib.Path(local_file).is_file():
            if not root_dir:
                r_code = shutil.copy(local_file, out_data_dir)
                rcode_lst.append(r_code)
            else:
                dest_file = local_file
                dest_file = dest_file.replace(root_dir, out_data_dir)
                print(local_file, dest_file)
                r_code = DU.copy_file(local_file, dest_file)
                # shutil.copy(local_file, out_data_dir)
                rcode_lst.append(r_code)
        else:
            rcode_lst.append(-1)
    if file_fmt.lower() == "png":
        rdx = pd.DataFrame(data={
            'RowId' : cel_df['RowId'], 'CEL' : cel_df['CEL'],
            'PNG' : cel_df['PNG'], 'ReturnCode': rcode_lst})
    elif file_fmt.lower() == "png":
        rdx = pd.DataFrame(data={
            'RowId' : cel_df['RowId'], 'CEL' : cel_df['CEL'],
            'PDF' : cel_df['PDF'], 'ReturnCode': rcode_lst})
    else:
        rdx = pd.DataFrame(data={
            'RowId' : cel_df['RowId'], 'CEL' : cel_df['CEL'],
            'ReturnCode': rcode_lst})
    rdx.to_csv(out_file)


if __name__ == '__main__':
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("in_file", help="Input csv file with cols RowId, CEL, and PDF/PNG")
    PARSER.add_argument("data_dir", help="Destination directory")
    PARSER.add_argument("out_file", help="Output csv fils with status of copying")
    PARSER.add_argument("file_fmt", help="File format ", choices=['png', 'pdf', 'cel'])
    PARSER.add_argument("--root_dir", help="Root dir from the source file to be removed",
                        default="")
    ARGS = PARSER.parse_args()
    print(ARGS)
    main(ARGS.in_file, ARGS.data_dir, ARGS.out_file,
         ARGS.file_fmt, ARGS.root_dir)
