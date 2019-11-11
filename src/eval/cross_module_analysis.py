import argparse
import tempfile
import os
import pandas as pd
import compare_gsnetwork as cgs

def cross_module_edges(mod1_genes, mod2_genes):
    rdf = pd.DataFrame(
        [(x, y) for x in mod1_genes for y in mod2_genes],
        columns=['TF', 'TARGET'])
    tfx = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False)
    cstmp = tfx.name
    tfx.close()
    #print(cstmp)
    rdf.to_csv(cstmp, sep="\t", index=False)
    return cstmp

def in_module_genes(mod_genes):
    rdf = pd.DataFrame(
        [(x, y) for x in mod_genes for y in mod_genes if x < y],
        columns=['TF', 'TARGET'])
    tfx = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False)
    cstmp = tfx.name
    tfx.close()
    #print(cstmp)
    rdf.to_csv(cstmp, sep="\t", index=False)
    return cstmp


def eval_network_modpair(annot_file, mod1_file, mod2_file,
                         net_files, wt_attr, max_dist, max_edges,
                         reverse_order):
    with open(mod1_file) as mfx:
        mod1_genes = [lx.strip() for lx in mfx.readlines()]
    with open(mod2_file) as mfx:
        mod2_genes = [lx.strip() for lx in mfx.readlines()]
    gs_cross_fx = cross_module_edges(mod1_genes, mod2_genes)
    gs_mod1_fx = in_module_genes(mod1_genes)
    gs_mod2_fx = in_module_genes(mod2_genes)
    if annot_file is None:
        cross_nhdat = [cgs.eval_network_ids(fx, gs_cross_fx, wt_attr,
                                            max_dist, max_edges, reverse_order)
                       for fx in net_files]
        mod1_nhdat = [cgs.eval_network_ids(fx, gs_mod1_fx, wt_attr,
                                           max_dist, max_edges, reverse_order)
                      for fx in net_files]
        mod2_nhdat = [cgs.eval_network_ids(fx, gs_mod2_fx, wt_attr,
                                           max_dist, max_edges, reverse_order)
                      for fx in net_files]
    else:
        cross_nhdat = [cgs.eval_network_probes(annot_file, fx, gs_cross_fx, wt_attr,
                                               max_dist, max_edges, reverse_order)
                       for fx in net_files]
        mod1_nhdat = [cgs.eval_network_probes(annot_file, fx, gs_mod1_fx, wt_attr,
                                              max_dist, max_edges, reverse_order)
                      for fx in net_files]
        mod2_nhdat = [cgs.eval_network_probes(annot_file, fx, gs_mod2_fx, wt_attr,
                                              max_dist, max_edges, reverse_order)
                      for fx in net_files]
    full_clnames = []
    full_clnames += ['MOD1_GENES', 'MOD2_GENES']
    full_clnames += ['NVRT', 'NEDG', 'NDENS', 'MOD1', 'MOD2']
    full_clnames += ['CRX_EDGES', 'MOD1_EDGES', 'MOD2_EDGES']
    full_clnames += ['MOD1_NETG', 'MOD2_NETG']
    full_clnames += ['CRX_NETEDG', 'MOD1_NETEDG', 'MOD2_NETEDG']
    full_clnames += ['CRX_EDGN'+str(y) for y in range(1, max_dist+1)]
    full_clnames += ['MOD1_EDGN'+str(y) for y in range(1, max_dist+1)]
    full_clnames += ['MOD2_EDGN'+str(y) for y in range(1, max_dist+1)]
    full_clnames += ['CRX_TP', 'MOD1_TP', 'MOD2_TP']
    gs_cmp_data = {str(net_files[x]) :
                   [mod1_file, mod2_file] +
                   [cross_nhdat[x]['NVRT'], cross_nhdat[x]['NEDG'], cross_nhdat[x]['NDENS']] +
                   [cross_nhdat[x]['GSNETID'][0], cross_nhdat[x]['GSNETID'][1]] +
                   [cross_nhdat[x]['GSNETID'][2], mod1_nhdat[x]['GSNETID'][2],
                    mod2_nhdat[x]['GSNETID'][2]] +
                   [cross_nhdat[x]['GRGS'][0][0], cross_nhdat[x]['GRGS'][0][1]] +
                   [cross_nhdat[x]['GRGS'][0][3], mod1_nhdat[x]['GRGS'][0][3],
                    mod2_nhdat[x]['GRGS'][0][3]] +
                   [cross_nhdat[x]['EDGN'][y] for y in range(1, max_dist+1)] +
                   [mod1_nhdat[x]['EDGN'][y] for y in range(1, max_dist+1)] +
                   [mod2_nhdat[x]['EDGN'][y] for y in range(1, max_dist+1)] +
                   [cross_nhdat[x]['FPCTS'][0][1], mod1_nhdat[x]['FPCTS'][0][1],
                    mod2_nhdat[x]['FPCTS'][0][1]] 
                   for x in range(len(net_files))
                   }
    gs_cmp_data_df = pd.DataFrame(data=gs_cmp_data, index=full_clnames)
    print(gs_cmp_data_df.to_csv(sep='\t', index=True))
    os.unlink(gs_cross_fx)
    os.unlink(gs_mod1_fx)
    os.unlink(gs_mod2_fx)


if __name__ == "__main__":
    PROG_DESC = """
    Cross module analysis
    """
    PARSER = argparse.ArgumentParser(description=PROG_DESC)
    PARSER.add_argument("annotation_file",
                        help="""annotation file
                                (a tab seperated file mapping probe to ids)""")
    PARSER.add_argument("mod1_file",
                        help="""genes that should be in module 1""")
    PARSER.add_argument("mod2_file",
                        help="""genes that should be in module 2""")
    PARSER.add_argument("reveng_network_files", nargs="+",
                        help="""network build from a reverse engineering methods
                                (currenlty supported: eda, adj, tsv)""")
    PARSER.add_argument("-d", "--dist", type=int, default=3,
                        help="max. number of hops allowed")
    PARSER.add_argument("-t", "--wt_attr", type=str, default='wt',
                        help="name of weight attribute")
    PARSER.add_argument("-x", "--max_edges", type=int,
                        help="""Maximum number of edges""")
    PARSER.add_argument("-r", "--reverse_order", action='store_true',
                        help="""Order the edges ascending order""")
    ARGS = PARSER.parse_args()
    eval_network_modpair(
        ARGS.annotation_file if ARGS.annotation_file != "-" else None,
        ARGS.mod1_file, ARGS.mod2_file, ARGS.reveng_network_files,
        ARGS.wt_attr, ARGS.dist, ARGS.max_edges,
        ARGS.reverse_order)
