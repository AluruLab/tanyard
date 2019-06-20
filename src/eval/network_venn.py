
import argparse
import venn
from data_utils import load_reveng_network

def main(network_files: str, network_names: str, out_file) -> None:
    if len(network_files) > 4:
        network_files = network_files[0:4]
    if network_names:
        network_names = network_names.split(",")
        if len(network_names) > 4:
            network_names = network_names[0:4]
    if network_names is None:
        network_names = ['net_' + str(ix) for ix in range(len(network_files))]
    network_dfs = [load_reveng_network(nx) for nx in network_files]
    st_lists = [ndf.loc[:, ['source', 'target']].values.to_list() for ndf in network_dfs]
    venn_labels = venn.get_labels(st_lists)
    print(venn_labels)


if __name__ == "__main__":
    PROG_DESC = """Network venn diagram"""
    ARGPARSER = argparse.ArgumentParser(description=PROG_DESC)
    ARGPARSER.add_argument("-n", "--network_names", type=str,
                        help="""comma seperated names of the network;
                                should have as many names as the number of networks""")
    ARGPARSER.add_argument("-o", "--out_file",
                        type=argparse.FileType(mode='w'), required=True,
                        help="output file in tab-seperated format")
    ARGPARSER.add_argument("network_files", nargs="+")
    CMDARGS = ARGPARSER.parse_args()
    print("""
       ARG : network_files : %s
       ARG : network_names : %s
       ARG : out_file : %s """ %
          (str(CMDARGS.network_files), str(CMDARGS.network_names),
           str(CMDARGS.out_file)))
    main(CMDARGS.network_files, CMDARGS.network_names, CMDARGS.out_file)
