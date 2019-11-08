import intermine.webservice as imweb
import pandas as pd
import argparse

ARAPORT_URL = "https://apps.araport.org/thalemine/service"

def enrich_entry(mod_item, mod_id, ngenes, enrich_type):
    mod_item['MODULE_ID'] = mod_id
    mod_item['MODULE_GENES'] = ngenes
    mod_item['ENRICHMENT'] = enrich_type
    return mod_item

def main(araport_token, module_file, output_file):
    module_df = pd.read_csv(module_file, sep="\t")
    service = imweb.Service(ARAPORT_URL, token=araport_token)
    num_modules = max(module_df.MODULE_ID)
    enrich_data = []
    for mod_id in range(1, num_modules+1):
        mod_genes = list(module_df.loc[module_df.MODULE_ID == mod_id,
                                       'GENE_ID'])
        lman = service.list_manager()
        mod_list = lman.get_list("Module List")
        if mod_list is not None:
            lman.delete_lists([mod_list])
        mod_list = lman.create_list(content=mod_genes, list_type="Gene",
                                    name="Module List")
        ngenes = len(mod_genes)
        enrich_data += [enrich_entry(item, mod_id, ngenes, 'Pathway Enrichment')
                        for item in mod_list.calculate_enrichment("pathway_enrichment")]
        enrich_data += [enrich_entry(item, mod_id, ngenes, 'GO Enrichment')
                        for item in mod_list.calculate_enrichment("go_enrichment_for_gene")]
        lman.delete_lists([mod_list])
    out_df = pd.DataFrame(enrich_data)
    out_df.to_csv(output_file, sep="\t")


if __name__ == "__main__":
    PROG_DESC = """
    Analyze modules with ARAPORT
    """
    PARSER = argparse.ArgumentParser(description=PROG_DESC)
    PARSER.add_argument("araport_token",
                        help="""Token for loggin in ARAPORT""")
    PARSER.add_argument("module_file",
                        help="""annotation file
                                (a tab seperated file mapping probe to ids)""")
    PARSER.add_argument("output_file", 
                        help="""output file""")
    ARGS = PARSER.parse_args()
    main(ARGS.araport_token, ARGS.module_file, ARGS.output_file)
