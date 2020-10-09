import argparse
import pandas as pd
import numpy as np
import networkx as nx
import data_utils as du

def shortest_path(net_graph, src, tgt):
    spath = None
    if src in net_graph and tgt in net_graph:
        try:
            spath = (nx.shortest_path(net_graph, src, tgt), src, tgt)
        except nx.NetworkXNoPath:
            spath = (None, src, tgt)
    else:
        if src in net_graph:
            spath = (None, src, None)
        elif tgt in net_graph:
            spath = (None, None, tgt)
        else:
            spath = (None, None, None)
    return spath

def select_edges(net_df: pd.DataFrame, wt_attr_name: str = 'wt',
                 max_edges: int = None,
                 reverse_order: bool = False):
    if max_edges is None or max_edges >= net_df.shape[0]:
        return net_df
    cur_cols = net_df.columns
    maxwt_attr_name = wt_attr_name + '_max'
    net_df[maxwt_attr_name] = net_df[wt_attr_name].abs()
    if reverse_order is True:
        net_df = net_df.nsmallest(n=max_edges, columns=maxwt_attr_name)
    else:
        net_df = net_df.nlargest(n=max_edges, columns=maxwt_attr_name)
    #print(net_df.iloc[0, :])
    return net_df.loc[:, cur_cols]

def zero_div(numerator, denominator):
    if denominator == 0:
        return 0
    return numerator/denominator


def compute_tp_fp(gs_tru_edges, gs_fls_edges, rv_net_graph,
                  wt_attr, reverse_order):
    fp_cts = 0
    tp_cts = 0
    rv_tp_fp_edges = []
    for udx, vdx in gs_fls_edges:
        if udx in rv_net_graph and vdx in rv_net_graph:
            try:
                if nx.shortest_path_length(rv_net_graph, udx, vdx) == 1:
                    rv_tp_fp_edges.append([udx, vdx, False,
                                           rv_net_graph[udx][vdx][wt_attr],
                                           0.0, 0.0])
                    fp_cts += 1
            except nx.NetworkXNoPath:
                pass
    for udx, vdx in gs_tru_edges:
        if udx in rv_net_graph and vdx in rv_net_graph:
            try:
                if nx.shortest_path_length(rv_net_graph, udx, vdx) == 1:
                    rv_tp_fp_edges.append([udx, vdx, True,
                                           rv_net_graph[udx][vdx][wt_attr],
                                           0, 0])
                    tp_cts += 1
            except nx.NetworkXNoPath:
                pass
    prec = (float(tp_cts)/(tp_cts + fp_cts)) if (tp_cts + fp_cts) > 0 else 0
    # recall 
    recall = (float(tp_cts))/len(gs_tru_edges) if len(gs_tru_edges) > 0 else 0
    # aupr
    sorted(rv_tp_fp_edges, key=lambda x: x[3], reverse=reverse_order is False)
    partial_tp_cts = 0
    partial_fp_cts = 0
    for idx, (udx, vdx, flag, _, _, _) in enumerate(rv_tp_fp_edges):
        if flag is True:
            partial_tp_cts += 1
        else:
            partial_fp_cts += 1
        if len(gs_tru_edges) > 0:
            rv_tp_fp_edges[idx][4] = (float(partial_tp_cts))/len(gs_tru_edges)
        if (tp_cts + fp_cts) > 0:
            rv_tp_fp_edges[idx][5] = (float(partial_tp_cts)/(partial_tp_cts + partial_fp_cts))
    # TODO: For each point in rv_tp_fp_edges
    #   compute precision, recall and add in an array
    prc_xcoords = [x for (_, _, _, _, x, _ ) in rv_tp_fp_edges]
    prc_ycoords = [y for (_, _, _, _, _, y ) in rv_tp_fp_edges]
    aupr = np.trapz(prc_ycoords, prc_xcoords)
    print("AUPR", max(prc_ycoords), min(prc_ycoords), max(prc_xcoords),
                  min(prc_xcoords), aupr)
    # call trapz to find area under that array
    return tp_cts, fp_cts, prec, recall, aupr

def shorest_path_graph(gs_net, rv_net_graph, max_dist,
                       tf_col='TFPROBE', tgt_col='TARGETPROBE'):
    gs_spath = [shortest_path(rv_net_graph, x, y)
                for x, y in zip(gs_net.loc[:, tf_col], gs_net.loc[:, tgt_col])]
    #gs_spath.sort()
    dist_histogram = [0 for x in range(max_dist+1)]
    spath_graph_nodes = set(x for _, x, _ in gs_spath if x)
    spath_graph_nodes |= set(y for _, _, y in gs_spath if y)
    for spx, _, _ in gs_spath:
        if spx is not None and (len(spx)-1) <= max_dist:
            dist_histogram[len(spx)-1] += 1
            spath_graph_nodes.update(spx)
        #if spx is not None and (len(spx)-1) <= 1:
        #    print((x, y) if x < y else (y, x))
    spath_graph = nx.Graph()
    spath_graph.add_nodes_from(spath_graph_nodes)
    for spx, _, _ in gs_spath:
        if spx is not None and (len(spx)-1) <= max_dist:
            if len(spx) > 1:
                for s_node, t_node in zip(spx, spx[1:]):
                    spath_graph.add_edge(s_node, t_node)
    return dist_histogram, spath_graph


def eval_network(rv_net, gs_net, max_dist, wt_attr,
                 tf_col, tgt_col, reverse_order):
    rv_net_graph = nx.from_pandas_edgelist(rv_net, edge_attr=wt_attr)
    rv_net_nodes = set(rv_net.source) | set(rv_net.target)
    gs_nedges = gs_net.shape[0]
    gs_tf_nodes = set(gs_net.loc[:, tf_col]) & rv_net_nodes
    gs_tgt_nodes = set(gs_net.loc[:, tgt_col]) & rv_net_nodes
    gs_pos_edges = set((x, y) for x in gs_tf_nodes for y in gs_tgt_nodes if x != y)
    gs_tru_edges = set((x, y) for x, y in zip(gs_net.loc[:, tf_col], gs_net.loc[:, tgt_col]))
    gs_fls_edges = list(gs_pos_edges - gs_tru_edges)
    gs_nodes = gs_tf_nodes | gs_tgt_nodes
    gs_common_nodes = sum((1 if x in rv_net_graph else 0 for x in gs_nodes))
    gs_common_edges = sum((1 if (x in rv_net_graph and y in rv_net_graph) else 0
                           for x, y in zip(gs_net.loc[:, tf_col], gs_net.loc[:, tgt_col])))
    #print(rv_net.shape, len(gs_nodes), gs_nedges, len(gs_tru_edges),
    #    gs_common_nodes, gs_common_edges)
    dist_histogram, spath_graph = shorest_path_graph(gs_net, rv_net_graph,
                                                     max_dist, tf_col, tgt_col)
    tp_cts, fp_cts, prec, recall, aupr = compute_tp_fp(gs_tru_edges, gs_fls_edges,
                                                       rv_net_graph, wt_attr, reverse_order)
    hist_data = {
        'NVRT'  : nx.number_of_nodes(rv_net_graph),
        'NEDG'  : nx.number_of_edges(rv_net_graph),
        'NDENS' : nx.density(rv_net_graph),
        'DIST'  : [x for x in range(max_dist+1)],
        'EDGN'  : dist_histogram,
        'FPCTS' : [(fp_cts, tp_cts, prec, prec*100, 
                    recall, recall*100, aupr) for _ in range(max_dist+1)],
        'PCTGS' : [zero_div(float(x)*100, gs_nedges) for x in dist_histogram],
        'PCTCM' : [zero_div(float(x)*100, gs_common_edges) for x in dist_histogram],
        'PCTSP' : [zero_div(float(x)*100, spath_graph.number_of_edges())
                   for x in dist_histogram],
        'GRSP'  : [(spath_graph.number_of_nodes(), spath_graph.number_of_edges())
                   for _ in range(max_dist+1)],
        'GRGS'  : [(len(gs_tf_nodes), len(gs_tgt_nodes), len(gs_nodes), gs_common_edges)
                   for _ in range(max_dist+1)],
        'GRCM'  : [(str(gs_common_nodes), str(gs_common_edges))
                   for _ in range(max_dist+1)]
    }
    #dist_histogram_df = pd.DataFrame(data=hist_data)
    #print(net_file)
    #print(dist_histogram_df)
    return hist_data


def eval_network_probes(annot_file, net_file, gs_file,
                        wt_attr, max_dist, max_edges, reverse_order):
    annot_df = du.load_annotation(annot_file)
    gs_net = du.load_gsnetwork(gs_file)
    gs_net = gs_net.loc[:, ['TF', 'TARGET']]
    # print("Network size (Before Probe Mapping)", gs_net.shape)
    id_net_shape = [len(set(gs_net.TF)), len(set(x for x in gs_net.TARGET)), gs_net.shape[0]]
    gs_net = du.map_probes(gs_net, annot_df)
    # print("Network size (After Probe Mapping)", gs_net.shape)
    rv_net = select_edges(du.load_reveng_network(net_file, wt_attr_name=wt_attr),
                          wt_attr_name=wt_attr,
                          max_edges=max_edges,
                          reverse_order=reverse_order)
    hist_data = eval_network(rv_net, gs_net, max_dist, wt_attr,
                             'TFPROBE', 'TARGETPROBE', reverse_order)
    hist_data['GSNETID'] = id_net_shape
    return hist_data

def compare_eval_network_probes(annot_file, net_files, gs_file, wt_attr,
                                max_dist, max_edges, reverse_order,
                                prefix):
    nhdat = [eval_network_probes(annot_file, fx, gs_file, wt_attr,
                                 max_dist, max_edges, reverse_order)
             for fx in net_files]
    clnames1 = ['NVRT', 'NEDG', 'NDENS', 'GSTFS', 'GSTARGET', 'GSEDGES'] + [
        'GRGSTF', 'GRGSTGT', 'GTGSV', 'GRGSE', # 'GRCMV', 'GRCME', 'GRSPV', 'GRSPE'
        ] + [
            'EDGN'+str(y) for y in range(1, max_dist+1)] + [
                'FP', 'TP', 'PREC', 'PRECPCT', 'RECALL', 'RECALLPCT', 'AUPR']
    clnames = [ prefix + "-" + cx  for cx in clnames1]
    gs_cmp_data = {str(net_files[x]) :
                   [nhdat[x]['NVRT'], nhdat[x]['NEDG'], nhdat[x]['NDENS']] +
                   [nhdat[x]['GSNETID'][0], nhdat[x]['GSNETID'][1], nhdat[x]['GSNETID'][2]] +
                   [nhdat[x]['GRGS'][0][0], nhdat[x]['GRGS'][0][1], nhdat[x]['GRGS'][0][2],
                    nhdat[x]['GRGS'][0][3]] +
                   #[nhdat[x]['GRCM'][0][0], nhdat[x]['GRCM'][0][1]] +
                   #[nhdat[x]['GRSP'][0][0], nhdat[x]['GRSP'][0][1]] +
                   [nhdat[x]['EDGN'][y] for y in range(1, max_dist+1)] +
                   [nhdat[x]['FPCTS'][0][0], nhdat[x]['FPCTS'][0][1],
                    nhdat[x]['FPCTS'][0][2], nhdat[x]['FPCTS'][0][3],
                    nhdat[x]['FPCTS'][0][4], nhdat[x]['FPCTS'][0][5],
                    nhdat[x]['FPCTS'][0][6]]
                   for x in range(len(net_files))}
    gs_cmp_data_df = pd.DataFrame(data=gs_cmp_data, index=clnames)
    print(gs_cmp_data_df.to_csv(sep='\t', index=True))


def eval_network_ids(net_file, gs_file, wt_attr,
                     max_dist, max_edges, reverse_order):
    gs_net = du.load_gsnetwork(gs_file)
    id_net_shape = [len(set(gs_net.TF)), len(set(gs_net.TARGET)), gs_net.shape[0]]
    rv_net_df = du.load_tsv_network(net_file, wt_attr_name=wt_attr)
    rv_net = select_edges(rv_net_df,
                          wt_attr_name=wt_attr,
                          max_edges=max_edges,
                          reverse_order=reverse_order)
    hist_data = eval_network(rv_net, gs_net, max_dist, wt_attr,
                             'TF', 'TARGET', reverse_order)
    hist_data['GSNETID'] = id_net_shape
    return hist_data

def compare_eval_network_ids(net_files, gs_file, wt_attr,
                             max_dist, max_edges, reverse_order,
                             prefix):
    nhdat = [eval_network_ids(fx, gs_file, wt_attr,
                              max_dist, max_edges, reverse_order)
             for fx in net_files]
    clnames1 = ['NVRT', 'NEDG', 'NDENS', 'GSTFS', 'GSTARGET', 'GSEDGES'] + [
        'GRGSTF', 'GRGSTGT', 'GTGSV', 'GRGSE', # 'GRCMV', 'GRCME', 'GRSPV', 'GRSPE'
        ] + [
            'EDGN'+str(y) for y in range(1, max_dist+1)] + [
                'FP', 'TP', 'PREC', 'PRECPCT']
    clnames = [ prefix + "-" + cx  for cx in clnames1]
    gs_cmp_data = {str(net_files[x]) :
                   [nhdat[x]['NVRT'], nhdat[x]['NEDG'], nhdat[x]['NDENS']] +
                   [nhdat[x]['GSNETID'][0], nhdat[x]['GSNETID'][1], nhdat[x]['GSNETID'][2]] +
                   [nhdat[x]['GRGS'][0][0], nhdat[x]['GRGS'][0][1], nhdat[x]['GRGS'][0][2],
                    nhdat[x]['GRGS'][0][3]] +
                   # [nhdat[x]['GRCM'][0][0], nhdat[x]['GRCM'][0][1]] +
                   # [nhdat[x]['GRSP'][0][0], nhdat[x]['GRSP'][0][1]] +
                   [nhdat[x]['EDGN'][y] for y in range(1, max_dist+1)] +
                   [nhdat[x]['FPCTS'][0][0], nhdat[x]['FPCTS'][0][1],
                    nhdat[x]['FPCTS'][0][2], nhdat[x]['FPCTS'][0][3]]
                   for x in range(len(net_files))}
    gs_cmp_data_df = pd.DataFrame(data=gs_cmp_data, index=clnames)
    print(gs_cmp_data_df.to_csv(sep='\t', index=True))


def compare_eval_network_probe_ranges(annot_file, net_files, gs_file,
                                      wt_attr, max_dist,
                                      eranges, step, reverse_order,
                                      prefix):
    lstart = int(eranges.split(",")[0])
    rend = int(eranges.split(",")[1])
    edge_range = list(range(lstart, rend+1, step))
    print(edge_range)
    nhdat = [eval_network_probes(annot_file, fx, gs_file,
                                 wt_attr, max_dist, nedges, reverse_order)
             for nedges in edge_range for fx in net_files]
    print(len(nhdat))
    # clnames = ['NVRT', 'NEDG', 'NDENS', 'GRGSV', 'GRGSE',
    #     'GRCMV', 'GRCME', 'GRSPV', 'GRSPE'] + [
    #     'EDGN'+str(y) for y in range(max_dist+1)] + [
    #         'FP', 'TP', 'PREC', 'PRECPCT']
    clnames1 = ['NVRT', 'NEDG', 'NDENS', 'GSTFS', 'GSTARGET', 'GSEDGES'] + [
        'EDGN'+str(y) for y in range(max_dist+1)] + [
            'FP', 'TP', 'PREC', 'RCALL']
    clnames = [ prefix + "-" + cx  for cx in clnames1]
    nranges = len(edge_range)
    nfiles = len(net_files)
    nentries = nranges * nfiles
    gs_cmp_data = {str(net_files[0]) + "-" + str(edge_range[x % nranges]):
                   [nhdat[x]['NVRT'], nhdat[x]['NEDG'], nhdat[x]['NDENS']] +
                   [nhdat[x]['GSNETID'][0], nhdat[x]['GSNETID'][1], nhdat[x]['GSNETID'][2]] +
                   # [nhdat[x]['GRGS'][0][0], nhdat[x]['GRGS'][0][1]] +
                   # [nhdat[x]['GRCM'][0][0], nhdat[x]['GRCM'][0][1]] +
                   # [nhdat[x]['GRSP'][0][0], nhdat[x]['GRSP'][0][1]] +
                   [nhdat[x]['EDGN'][y] for y in range(max_dist+1)] +
                   [nhdat[x]['FPCTS'][0][0], nhdat[x]['FPCTS'][0][1],
                    nhdat[x]['FPCTS'][0][2],
                    float(nhdat[x]['FPCTS'][0][1])/float(nhdat[x]['GSNETID'][2])]
                   for x in range(nentries)}
    gs_cmp_data_df = pd.DataFrame(data=gs_cmp_data, index=clnames)
    print(gs_cmp_data_df.to_csv(sep='\t', index=True))


if __name__ == "__main__":
    PROG_DESC = """
    Compares the given networks against the gold standard networks and reports
    the number of hops in the input networks between two nodes that are present as
    edge in the gold standard network.
    """
    PARSER = argparse.ArgumentParser(description=PROG_DESC)
    PARSER.add_argument("annotation_file",
                        help="""annotation file
                                (a tab seperated file mapping probe to ids)""")
    PARSER.add_argument("gs_network_file",
                        help="""gold standard network
                                (tab seperated file of TF-TARGET interactions)""")
    PARSER.add_argument("reveng_network_files", nargs="+",
                        help="""network build from a reverse engineering methods
                                (currenlty supported: eda, adj, tsv)""")
    PARSER.add_argument("-d", "--dist", type=int, default=3,
                        help="max. number of hops allowed")
    PARSER.add_argument("-t", "--wt_attr", type=str, default='wt',
                        help="name of weight attribute")
    PARSER.add_argument("-x", "--max_edges", type=int,
                        help="""Maximum number of edges""")
    PARSER.add_argument("-n", "--range_edges", type=str,
                        help="""Range of edges""")
    PARSER.add_argument("-s", "--range_steps", type=int, default=100000,
                        help="""Steps of edges""")
    PARSER.add_argument("-r", "--reverse_order", action='store_true',
                        help="""Order the edges ascending order""")
    PARSER.add_argument("-p", "--prop_prefix", type=str, default="PROP",
                        help="""Prefix String for Properties""")
    ARGS = PARSER.parse_args()
    if ARGS.range_edges is None:
        if ARGS.annotation_file == "-":
            compare_eval_network_ids(ARGS.reveng_network_files,
                                     ARGS.gs_network_file, ARGS.wt_attr,
                                     ARGS.dist, ARGS.max_edges,
                                     ARGS.reverse_order,
                                     ARGS.prop_prefix)
        else:
            compare_eval_network_probes(ARGS.annotation_file,
                                        ARGS.reveng_network_files,
                                        ARGS.gs_network_file,
                                        ARGS.wt_attr,
                                        ARGS.dist,
                                        ARGS.max_edges,
                                        ARGS.reverse_order,
                                        ARGS.prop_prefix)
    else:
        if ARGS.annotation_file == "-":
            #compare_eval_network_ids(ARGS.reveng_network_files,
            #                         ARGS.gs_network_file, ARGS.wt_attr,
            #                         ARGS.dist, ARGS.max_edges)
            print("Option Not supported Yet!")
        else:
            compare_eval_network_probe_ranges(ARGS.annotation_file,
                                              ARGS.reveng_network_files,
                                              ARGS.gs_network_file,
                                              ARGS.wt_attr,
                                              ARGS.dist,
                                              ARGS.range_edges,
                                              ARGS.range_steps,
                                              ARGS.reverse_order,
                                              ARGS.prop_prefix)
