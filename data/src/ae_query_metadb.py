import os
import os.path
import json
import pprint
import pandas as pd

def ae_raw_files(fxnd):
    return [x['url'] for x in fxnd['file'] if x['kind'] == 'raw']

def load_json_file(fname):
    with open(fname) as jsnf:
        return json.load(jsnf)
    
def ae_expt_accession(expt_node):
    return expt_node['accession']

def ae_expt_name(expt_node):
    return expt_node['name']

def ae_expt_type(expt_node):
    return ";".join(expt_node['experimenttype'])

def ae_expt_desc(expt_node):
    return ";".join([x['text'] for x in expt_node['description']])

def ae_expt_sub_date(expt_node):
    return expt_node['releasedate']

def ae_expt_upd_date(expt_node):
    return expt_node['lastupdatedate']

def ae_expt_pubmed(expt_node):
    if 'bibliography' in expt_node:
        return ";".join([str(x['accession']) if 'accession' in x else ''
                         for x in expt_node['bibliography']])
    else:
        ''

def ae_expt_organism(expt_node):
    return ";".join(expt_node['organism'])

def ae_sample_name(smnx):
    return smnx['source']['name']

def ae_sample_attr_dict(smnx):
    var_dict = {}
    if 'variable' in smnx:
        var_dict = { x['name'] : str(x['value']) for x in smnx['variable']  }
    char_dict = {}
    if 'characteristic' in smnx:
        char_dict = { x['category'] : x['value'] for x in smnx['characteristic'] }
    return {**var_dict, **char_dict}

def ae_sample_attrs(sm_attrs):
    return ";".join([ "{}: [{}]".format(x, y) 
                     for x, y in sm_attrs.items() ])

def ae_sample_file(smnx, fxnd):
    raw_files = ae_raw_files(fxnd)
    fn_list = []
    if 'file' in smnx:
        if len(raw_files) == 1:
            for x in smnx['file']:
                try:
                    if x['type'] == 'data' and x['name'] is not None:
                        fn_list.append(raw_files[0] + "/" + x['name'])
                except TypeError:
                    print(x)
        if len(fn_list) == 0:
            for x in smnx['file']:
                if 'url' in x:
                    fn_list.append(x['url'])
    return  ";".join(fn_list)


def ae_expts_data(expts_json_file, samples_json_file, files_json_file):
    pp = pprint.PrettyPrinter(indent=4)
    ex = load_json_file(expts_json_file)
    sx = load_json_file(samples_json_file)
    fx = load_json_file(files_json_file)
    exnd = ex['experiments']['experiment'][0]
    fxnd = fx['files']['experiment']
    smxnd = sx['experiment']['sample']
    expts_data = []
    try:
        acc_id = ae_expt_accession(exnd)
        series_data = {
         'SeriesId' :  acc_id,
         'SeriesTitle' :  ae_expt_name(exnd),
         'SeriesExperimentType' :  ae_expt_type(exnd),
         'SeriesDescription' :  ae_expt_desc(exnd),
         'SeriesSubmissionDate' :  ae_expt_sub_date(exnd),
         'SeriesUpdateDate' :  ae_expt_upd_date(exnd),
         'SeriesPubMedID' :  ae_expt_pubmed(exnd),
         'SeriesLink' : 'https://www.ebi.ac.uk/arrayexpress/experiments/{}'.format(
            acc_id),
         'SeriesOrganism' :  ae_expt_organism(exnd)}
    except KeyError as kerr:
        pp.pprint(exnd)
        raise kerr
    for smx in smxnd:
        try:
            sm_attrs = ae_sample_attr_dict(smx)
            smp_data = {'SampleId' :  ae_sample_name(smx),
                        'SampleOrganism' :  sm_attrs['Organism'] if 'Organism' in sm_attrs else '',
                        'SampleAttributes' :  ae_sample_attrs(sm_attrs),
                        'SampleFile' :  ae_sample_file(smx, fxnd)}
            expts_data.append({**series_data, **smp_data})
        except KeyError as kerr:
            print(acc_id)
            pp.pprint(smx)
            pp.pprint(exnd)
            raise kerr
    return expts_data

def main(ae_file, meta_ae_dir, df_out_file):
    aedf = pd.read_csv(ae_file, sep='\t', encoding = "ISO-8859-1")
    sum(aedf.Assays)
    expts_json_file_format = os.path.join(meta_ae_dir, "{}", "{}.expts.json")
    samples_json_file_format = os.path.join(meta_ae_dir, "{}", "{}.samples.json")
    files_json_file_format = os.path.join(meta_ae_dir, "{}", "{}.files.json")
    ae_expts_full = [
        ae_expts_data(
            expts_json_file_format.format(aeid, aeid),
            samples_json_file_format.format(aeid, aeid),
           files_json_file_format.format(aeid, aeid)
        ) for aeid in aedf['Accession']]
    pdf = pd.DataFrame([item for sublist in ae_expts_full for item in sublist])
    sel_cols = [
         'SeriesId',
         'SeriesTitle',
         'SeriesExperimentType',
         'SeriesDescription',
         'SeriesSubmissionDate',
         'SeriesUpdateDate',
         'SeriesPubMedID',
         'SeriesLink',
         'SeriesOrganism',
         'SampleId',
         'SampleOrganism',
         'SampleAttributes',
         'SampleFile'
            ]
    pdf2 = pdf.drop_duplicates()
    pdf3 = pdf2[sel_cols].sort_values(by=["SeriesId"])   

    pdf3.to_csv(df_out_file)




if __name__ == '__main__':
    ae_file = "../tables/AE-ATH1.tsv"
    meta_ae_dir = "../meta/AE/"
    df_out_file = "../tables/AE-META-ATH1.csv"
    main(ae_file, meta_ae_dir, df_out_file)

