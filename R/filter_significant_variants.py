#! /usr/bin/env python3

import argparse
import gzip

GZIP_EXTENSION = ".gz" # We need to know if input is compressed to make it compatible
READ_AS_TEXT = "rt"  # As well giving them the correct read
WRITE_TO_TEXT = "wt" # and write arguments

def parseargs():
    """
        We need 5 arguments to filter the original file, of which only 3 are mandatory
        Input file, output file, and which column hold pvalues

        in: The input file (path)
        out: The output file, where filtered variants will be sabed (results overwritten). By default,
            goes to stdout
        pcol: column index that hold pvalues
        pval: threshold, max p-val to be kept. By default, its 0.05
        sep: separator bewteen cells in the input file. The same one will be used in the output file
            Incorrect separator will produce an error. By default it's '\t' (tab)
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-pval", type=float, action="store", dest="pval", default=0.05)
    parser.add_argument("-in", type=str, action="store", dest="infile")
    parser.add_argument("-out", type=str, action="store", dest="outfile", default="/dev/stdout")
    parser.add_argument("-sep", type=str, action="store", dest="sep", default="\t")
    parser.add_argument("-pcol", type=int, action="store", dest="pcol")
    return parser.parse_args()

def loadcsv(file:str, sep:str):
    """
        A generator that reads one line at a time from the file, and splits 
        it according to the separator given. 
        
        Uses the correct function whether it's compressed in gz or not. 
        It's considered gzipped if has the extension present in GZIP_EXTENSION

        It's read as text according to the variable READ_AS_TEXT
    """
    open_function:function = gzip.open if file.endswith(GZIP_EXTENSION) else open
    for line in open_function(file, READ_AS_TEXT):
        yield line.strip("\n").split(sep)
    pass

def filter_by_pvalue(x:float, max_pval:float) -> bool:
    """
        Determines if a pvalue should be filtered or not

        Just checks both are floats and compares them
    """
    pval = float(x)
    return pval < max_pval

def write_significant_variables(significant, outfile, sep) -> None:
    """
        Writes the registries in the output file using the separator sep
        as a gzipped file
    """
    fhandler = gzip.open(outfile, WRITE_TO_TEXT)
    text:list = []
    for line in significant:
        fhandler.write(sep.join(line)+"\n")
    fhandler.close()

def main():
    args = parseargs()
    print("Loading variants")
    variants:list = loadcsv(args.infile, args.sep)
    print("Filtering variants")
    fhandler = gzip.open(args.outfile, WRITE_TO_TEXT)
    header = next(variants)
    fhandler.write(args.sep.join(header)+"\n")
    for v in variants: # Skip the header, just from 1 to end
        if (float(v[args.pcol]) < args.pval):
            fhandler.write(args.sep.join(v)+"\n") 
            #Quería hacerlo con la función filter, pero esta no admite generadores.
        # De esta forma es más eficiente en memoria pero no está tan funcionalizado
    fhandler.close()
    pass

if __name__ == "__main__":
    main()
