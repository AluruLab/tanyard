import argparse
import rpy2.robjects.packages
import rpy2.robjects.vectors


def main():
    r_base = rpy2.robjects.packages.importr('base')
    r_geomdb = rpy2.robjects.packages.importr('GEOmetadb')
    if r_base.file_exists("GEOmetadb.sqlite")[0] == False:
        r_geomdb.getSQLiteFile()

if __name__ == '__main__':
    prog_desc="""Downloads the latest GEOmetadb file
    """
    parser = argparse.ArgumentParser(description=prog_desc)
    parser.parse_args()
    main()