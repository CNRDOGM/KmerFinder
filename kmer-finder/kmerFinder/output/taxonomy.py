#!/usr/bin/env python3
'''  '''
from argparse import ArgumentParser
import subprocess
import sys
import os

def getTaxonomy():
    '''  '''
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", dest="infile", help="inputfile")
    parser.add_argument("-t", "--taxfile", dest="taxfile", help="taxonomy file")
    parser.add_argument("-b", "--bacteria", dest="bacteria",action="store_true", help="is it a bacteria DB?")
    parser.add_argument("-c", "--organism", dest="organism",action="store_true", help="is it an organism DB?")
    parser.add_argument("-o", "--outfile", dest="outfile", help="outputfile")
    args = parser.parse_args()
    # -------------------------------------------------------------
    #     parse commandline options
    # -------------------------------------------------------------
    
    # get infile:
    if args.infile is not None:
        infile = open(args.infile, "r")
    else:
        sys.stderr.write("Please specify inputfile!\n")
        sys.exit(2)
    
    # get infile:
    if args.taxfile is not None:
        taxfile = open(args.taxfile, "r")
    else:
        sys.stderr.write("Please specify path to taxonomy file!\n")
        sys.exit(2)
    
    # get infile:
    if args.outfile is not None:
        outfile = open(args.outfile, "w")
    else:
        sys.stderr.write("Please specify outputfile!\n")
        sys.exit(2)
    #------------------------------------------------------------
    #  store taxonomy
    #------------------------------------------------------------
    
    tax={}
    for l in taxfile:
        l=l.strip()
        if l == '': continue
        if l[0] == '#': continue
        l=l.split("\t")
        if args.organism == True:
            tax[l[1]] = "\t".join(l[2:])
        else:
            tax[l[0]] = "\t".join(l[2:])
    
    taxfile.close()
    
    #-------------------------------------------------------------
    #  parse kmerfinder output
    #-------------------------------------------------------------
    
    for l in infile:
        l = l.strip()
        if l == '': continue
        if l[0] == "#":
            if args.bacteria == True:
                outfile.write(l + "\tTAXID\tTaxonomy\tTAXID Species\tSpecies\n")
            else:
                outfile.write(l + "\tTAXID\tTaxonomy\n")
        else:
            tmp=l.split("\t")
            tmp=tmp[0]
            tmp=tmp.strip(">")
            if tmp[0:3] == "gi|":
                tmp=tmp.split("|")[1]
            print(tmp, tax.keys())
            if tmp in tax:
                t=tax[tmp]
            else:
                if args.bacteria == True:
                    t="unknown\tunknown\tunknown\tunknown"
                else:
                    t="unknown\tunknown"
            outfile.write(l + "\t" + t + "\n")
    
    infile.close()
    outfile.close()

if __name__ == "__main__":
    getTaxonomy()
