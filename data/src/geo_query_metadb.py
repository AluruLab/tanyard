import argparse
import sqlite3
import pandas as pd

def run_sql_query(cxn, sql_query):
    tresult = []
    try:
        cxcr = cxn.cursor()
        cxcr.execute(sql_query)
        tresult = cxcr.fetchall()
    except sqlite3.Error as ser:
        print(str(ser))
    try:
        cxcr.close()
    except sqlite3.Error as ser:
        print(str(ser))
    return tresult

def describe_table(cxn, table_name):
    tdesc = []
    try:
        cxcr = cxn.cursor()
        cxcr.execute('PRAGMA table_info({})'.format(table_name))
        tdesc = cxcr.fetchall()
    except sqlite3.Error as ser:
        print(str(ser))
    try:
        cxcr.close()
    except sqlite3.Error as ser:
        print(str(ser))
    return tdesc

def list_tables(cxn):
    return pd.read_sql_query("SELECT name FROM sqlite_master WHERE type='table'", cxn)

def describe_table_df(cxn, tname):
    return pd.DataFrame(describe_table(cxn, tname))[[1, 2]]


# Query GPL18's charactestics
#
# rqry = """select
#  distinct gsm.gpl, gpl.manufacturer, gpl.organism, gpl.title
# from
#  gse, gse_gpl, gse_gsm, gsm , gpl
# where
#  gsm.gpl != "GPL198" and
#  gse_gpl.gpl = "GPL198" and
#  gpl.gpl = gsm.gpl and
#  gse.gse = gse_gpl.gse and
#  gse_gsm.gse = gse_gpl.gse and
#  gse.gse = gse_gsm.gse and
#  gsm.gsm = gse_gsm.gsm
# """
# pd.read_sql_query(rqry, conn)

# Generic query over metadata
#
# rqry = """select
#  gse.gse as GSE_ID, gse.title as GSE_TITLE,
#  gse.status as GSE_STATUS,
#  gse.submission_date as GSE_SUBMISSION_DATE,
#  gse.last_update_date as GSE_UPDATE_DATE,
#  gse.pubmed_id as GSE_PUBMED,
#  gse.summary as GSE_SUMMARY,
#  gse.type as GSE_TYPE,
#  gse.web_link as GSE_LINK,
#  gse.overall_design as GSE_DESIGN,
#  gse.contact as GSE_CONTACT,
#  gse.supplementary_file as GSE_SFILE,
#  gsm.title as GSM_TITLE,
#  gsm.gsm as GSM_ID,
#  gsm.submission_date as GSM_SUBMISSION_DATE,
#  gsm.last_update_date as GSM_UPDATE_DATE,
#  gsm.type as GSM_TYPE,
#  gsm.source_name_ch1 as GSM_SOURCE,
#  gsm.organism_ch1 as GSM_ORGANISM,
#  gsm.characteristics_ch1 as GSM_CHARACTERSITCS,
#  gsm.molecule_ch1 as GSM_MOLECULE,
#  gsm.label_ch1 as GSM_LABEL,
#  gsm.treatment_protocol_ch1 as GSM_TREATMENT,
#  gsm.extract_protocol_ch1 as GSM_EXTRACTION,
#  gsm.label_protocol_ch1 as GSM_PROTOCOL_LABEL,
#  gsm.hyb_protocol as GSM_HYBRIDIZATION,
#  gsm.description as GSM_DESCRIPTION,
#  gsm.data_processing as GSM_PROCESSING,
#  gsm.supplementary_file as GSM_SFILE,
#  gsm.data_row_count as GSM_NROW
# from
#  gse, gse_gsm, gsm
# where
#  gsm.gpl = "GPL198" and
#  gse.gse = gse_gsm.gse and
#  gsm.gsm = gse_gsm.gsm
# """
# rdf = pd.read_sql_query(rqry, conn)

# print(len(rdf))

def query_meta_data_since(out_data_file, since_date):
    rqry = """select
    gse.gse as SeriesId,
    gse.title as SeriesTitle,
    gse.type as SeriesExperimentType,
    gse.summary as SeriesDescription,
    gse.submission_date as SeriesSubmissionDate,
    gse.last_update_date as SeriesUpdateDate,
    gse.status as SeriesStatus,
    gse.overall_design as SeriesOverallDesign,
    gse.pubmed_id as SeriesPubMedID,
    ('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' || gse.gse) as SeriesLink,
    gse.supplementary_file as SeriesFile,
    gsm.gsm as SampleId,
    gsm.title as SampleTitle,
    gsm.type as SampleType,
    gsm.description as SampleDescription,
    gsm.submission_date as SampleSubmissionDate,
    gsm.last_update_date as SampleUpdateDate,
    gsm.organism_ch1 as SampleOrganism,
    gsm.molecule_ch1 as SampleMolecule,
    ('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' || gsm.gsm) as SampleLink,
    (
    'SOURCE : [' || gsm.source_name_ch1 || ']; ' ||
    'CHARACTESTICS : [ ' || gsm.characteristics_ch1 || ']; ' ||
    'LABEL : [ ' || gsm.label_ch1 || ']; ' ||
    'TREATMENT : [ ' || gsm.treatment_protocol_ch1 || ']; ' ||
    'EXTRACTION : [ ' || gsm.extract_protocol_ch1 || ']; ' ||
    'PROTOCOL_LABEL : [ ' || gsm.label_protocol_ch1 || ']; ' ||
    'HYBRIDIZATION_LABEL : [ ' || gsm.hyb_protocol || ']; ' ||
    'GSM_PROCESSING : [ ' || gsm.data_processing || ']; '
    ) as SampleAttributes,
    gsm.supplementary_file as GSM_FILE
    from
    gse, gse_gsm, gsm
    where
    gsm.gpl = "GPL198" and
    gse.gse = gse_gsm.gse and
    gsm.gsm = gse_gsm.gsm and
    date(gsm.submission_date) > '{0}'
    """
    conn = sqlite3.connect("GEOmetadb.sqlite")
    rdf = pd.read_sql_query(rqry.format(since_date), conn)
    rdf.to_csv(out_data_file)
    conn.close()


def query_meta_data(out_data_file):
    rqry = """select
    gse.gse as SeriesId,
    gse.title as SeriesTitle,
    gse.type as SeriesExperimentType,
    gse.summary as SeriesDescription,
    gse.submission_date as SeriesSubmissionDate,
    gse.last_update_date as SeriesUpdateDate,
    gse.status as SeriesStatus,
    gse.overall_design as SeriesOverallDesign,
    gse.pubmed_id as SeriesPubMedID,
    ('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' || gse.gse) as SeriesLink,
    gse.supplementary_file as SeriesFile,
    gsm.gsm as SampleId,
    gsm.title as SampleTitle,
    gsm.type as SampleType,
    gsm.description as SampleDescription,
    gsm.submission_date as SampleSubmissionDate,
    gsm.last_update_date as SampleUpdateDate,
    gsm.organism_ch1 as SampleOrganism,
    gsm.molecule_ch1 as SampleMolecule,
    ('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' || gsm.gsm) as SampleLink,
    (
    'SOURCE : [' || gsm.source_name_ch1 || ']; ' ||
    'CHARACTESTICS : [ ' || gsm.characteristics_ch1 || ']; ' ||
    'LABEL : [ ' || gsm.label_ch1 || ']; ' ||
    'TREATMENT : [ ' || gsm.treatment_protocol_ch1 || ']; ' ||
    'EXTRACTION : [ ' || gsm.extract_protocol_ch1 || ']; ' ||
    'PROTOCOL_LABEL : [ ' || gsm.label_protocol_ch1 || ']; ' ||
    'HYBRIDIZATION_LABEL : [ ' || gsm.hyb_protocol || ']; ' ||
    'GSM_PROCESSING : [ ' || gsm.data_processing || ']; '
    ) as SampleAttributes,
    gsm.supplementary_file as GSM_FILE
    from
    gse, gse_gsm, gsm
    where
    gsm.gpl = "GPL198" and
    gse.gse = gse_gsm.gse and
    gsm.gsm = gse_gsm.gsm
    """
    conn = sqlite3.connect("GEOmetadb.sqlite")
    rdf = pd.read_sql_query(rqry, conn)
    rdf.to_csv(out_data_file)
    conn.close()


if __name__ == '__main__':
    PROG_DESC = """Query GEOmetadb.sqlite for the GEO datasets
    If datsets are required since a start date, provide the
    -s or --since_date argument in the format YYYY-MM-DD
    """
    PARSER = argparse.ArgumentParser(description=PROG_DESC)
    PARSER.add_argument("out_file")
    PARSER.add_argument("-s", "--since_date")
    ARGS = PARSER.parse_args()
    if ARGS.since_date:
        query_meta_data_since(ARGS.out_file,
                              ARGS.since_date)
    else:
        query_meta_data(ARGS.out_file)
