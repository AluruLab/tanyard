import argparse
import numpy as np
import pandas as pd


def main(network_file, output_file):
    net_df = pd.read_csv(network_file, sep="\t")
    th_lst = [float(x)/100 for x in range(1, 100)]
    edge_ct = [(sum(np.sum(abs(net_df) > x, axis=1)) - float(net_df.shape[0]))/2.0
               for x in th_lst]
    out_df = pd.DataFrame({"THRESHOLD": th_lst,
                           "EDGE_COUNT": edge_ct})
    out_df.to_csv(output_file)


if __name__ == "__main__":
    PROG_DESC = """Compute network density vs Threshold"""
    ARGPARSER = argparse.ArgumentParser(description=PROG_DESC)
    ARGPARSER.add_argument("network_file")
    ARGPARSER.add_argument("output_file")
    CMDARGS = ARGPARSER.parse_args()
    main(CMDARGS.network_file, CMDARGS.output_file)
