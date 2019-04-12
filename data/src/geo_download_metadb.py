import argparse
import rpy2.robjects.packages
import rpy2.robjects.vectors


def main():
    r_base = rpy2.robjects.packages.importr('base')
    r_geomdb = rpy2.robjects.packages.importr('GEOmetadb')
    if not r_base.file_exists("GEOmetadb.sqlite")[0]:
        r_geomdb.getSQLiteFile()

if __name__ == '__main__':
    PROG_DESC = """Downloads the latest GEOmetadb file
    using the GEOmetadb R package"""
    argparse.ArgumentParser(description=PROG_DESC).parse_args()
    main()
