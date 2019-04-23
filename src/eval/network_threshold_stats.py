import argparse
import numpy as np
import pandas as pd


def main(network_file, output_file):
    net_df = pd.read_csv(network_file, sep="\t", index_col=0)
    th_lst = [float(x)/100 for x in range(1, 100)]
    edge_ct = [0 for _ in range(1, 100)]
    node_ct = [0 for _ in range(1, 100)]
    for i,j in enumerate(range(1, 100)):
        x = j/100 
        edge_inc = np.sum(abs(net_df) > x, axis=1) - 1
        edge_ct[i] = sum(edge_inc)/2.0
        node_ct[i] = sum(edge_inc > 0)
    out_df = pd.DataFrame({"THRESHOLD": th_lst,
                           "EDGE_COUNT": edge_ct,
                           "NODE_COUNT": node_ct})
    out_df.to_csv(output_file)


if __name__ == "__main__":
    PROG_DESC = """Compute network density vs Threshold"""
    ARGPARSER = argparse.ArgumentParser(description=PROG_DESC)
    ARGPARSER.add_argument("network_file")
    ARGPARSER.add_argument("output_file")
    CMDARGS = ARGPARSER.parse_args()
    main(CMDARGS.network_file, CMDARGS.output_file)
