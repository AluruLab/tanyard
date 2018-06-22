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

def ae_sample_name(smnx): return smnx['source']['name']

def ae_sample_attr_dict(smnx): 
    return {**{ x['name'] : str(x['value']) for x in smnx['variable']  },
            **{ x['category'] : x['value'] for x in smnx['characteristic'] }}

def ae_sample_attrs(sm_attrs):
    return ";".join([ "{}: [{}]".format(x, y) 
                     for x, y in sm_attrs.items() ])

def ae_sample_file(smnx): return ";".join([x['url'] for x in smnx['file']])

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
    except KeyError as ke:
        pp.pprint(exnd)
        raise
    for smx in smxnd:
        try:
            sm_attrs = ae_sample_attr_dict(smx)
            smp_data = {'SampleId' :  ae_sample_name(smx),
                        'SampleOrganism' :  sm_attrs['Organism'] if 'Organism' in sm_attrs else '',
                        'SampleAttributes' :  ae_sample_attrs(sm_attrs),
                        'SampleFile' :  ae_sample_file(smx)}
            expts_data.append({**series_data, **smp_data})
        except KeyError as ke:
            print(acc_id)
            print(ae_raw_files(fxnd))
            pp.pprint(smx)
            pp.pprint(exnd)
            raise
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
    pdf.to_csv(df_out_file)




ae_file = "../tables/AE-ATH1.tsv"
meta_ae_dir = "../meta/AE/"

