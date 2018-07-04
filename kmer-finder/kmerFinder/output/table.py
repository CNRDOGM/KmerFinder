#!/usr/bin/env python3
'''  '''
from argparse import ArgumentParser
import subprocess
import sys
import os

def createTSV():
    '''  '''
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", dest="infile", help="inputfile")
    parser.add_argument("-d", "--database", dest="database",
                        help="which database was searched?")
    parser.add_argument("-o", "--outfile", dest="outfile", help="outputfile")
    args = parser.parse_args()
    
    # -------------------------------------------------------------
    #     parse commandline options
    # -------------------------------------------------------------
    
    # get infile:
    if args.infile != None:
        infile = open(args.infile, "r")
    else:
        sys.stderr.write("Please specify inputfile!\n")
        sys.exit(2)


    # get database:
    if args.database != None:
        database = args.database
    else:
        sys.stderr.write("Please specify database!\n")
        sys.exit(2)
    
    # open outputfile:
    if args.outfile != None:
        outfile = open(args.outfile, "w")
    else:
        sys.stderr.write("Please specify outputfile!\n")
        sys.exit(2)
    #------------------------------------------------------------
    #   parse kmerfinder + taxonomy results
    #------------------------------------------------------------
    w = False
    # if database.split("_")[0] == "bacteria":
    for l in infile:
        l = l.strip()
        # print header line:
        if l[0] == "#":
            l = l.split("\t")
            if "total depth" in l:
                outfile.write("#Hit\tScore\tz-score\tQuery Coverage [%]\t"
                              "Template Coverage [%]\tDepth\tTotal Query "
                              "Coverage [%]\tTotal Template Coverage [%]\tTotal"
                              " Depth\n")
                w = True
            else:
              outfile.write("#Hit\tScore\tz-score\tQuery Coverage [%]\t"
                            "Template Coverage [%]\tDepth\n")
        # parse hits:
        else:
            l = l.split("\t")
            # create the first column:
            taxlist = l[-3].split("; ")
            if l[-1].strip("\n") in taxlist and l[-1].strip("\n").lower() != "unknown":
                idx = taxlist.index(l[-1].strip("\n"))
                hit = ", ".join(taxlist[idx:len(taxlist) - 1])
            else:
                hit = l[0]
            sys.stderr.write(hit)
            #-------------------------------------------------------------------------
            #  CREATE LINK:
            #-------------------------------------------------------------------------
            
            # Bacteria
            if database.split("_")[0] == "bacteria":
                link = (" <a href='ftp://ftp.ncbi.nlm.nih.gov/genomes/"
                        "Bacteria/%s'> get sequence </a>")%(l[0])
            
            # Plasmids
            elif database == "plasmids":
                link = (" <a href='ftp://ftp.ncbi.nlm.nih.gov/genomes/"
                        "Plasmids/%s'> get sequence </a>")%(l[0])
            
            # Type strains:
            elif database == "type_strains":
              link = ("<a href='http://www.microbial-earth.org/%s'> get "
                      "sequence </a>")%(l[0])
            
            # Fungi
            elif database == "fungi":
                # NCBI complete genome:
                if l[-5][0:3] == "gi|":
                    link = (" <a href='ftp://ftp.ncbi.nlm.nih.gov/genomes/"
                            "Fungi/%s'> get sequence </a>")%(l[0])
                # NCBI assembly:
                elif l[0].split("_")[0] == "GCA":
                    #a_id = l[0].split("_")
                    #a_id = a_id[-2] + "_" + a_id[-1].strip(".fsa")
                    #a_id = a_id[0] + "_" + a_id[1]
                    a_id = l[0]
                    link = (" <a href='ftp://ftp.ncbi.nih.gov/genomes/all/%s'> "
                            "get sequence </a>")%(a_id)
                else:
                    link = (" <a href='ftp://ftp.ensemblgenomes.org/pub/fungi/"
                            "current/fasta/'> get sequence </a>")
            
            # Protozoa:
            elif database == "protists":
                # NCBI complete genome:
                if l[-5][0:3] == "gi|":
                    link = (" <a href='ftp://ftp.ncbi.nlm.nih.gov/genomes/"
                            "Protozoa/'> get sequence </a>")%(l[0])
                # NCBI assembly:
                elif l[0].split("_")[0] == "GCA":
                    #a_id = l[0].split("_")
                    #a_id = a_id[-2] + "_" + a_id[-1].strip(".fsa")
                    a_id = l[0]
                    link = (" <a href='ftp://ftp.ncbi.nih.gov/genomes/all/%s'>"
                            " get sequence </a>")%(a_id)
                else:
                    link = (" <a href='ftp://ftp.ensemblgenomes.org/pub/"
                            "protists/current/fasta/'> get sequence </a>")
            
            # Viruses:
            elif database.split("_")[0] == "virus":
                link = (" <a href='ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/'"
                        "> get sequence </a>")%(l[0])
            
            # ResFinder:
            elif database == "resfinder":
                link = ""
            hit = hit + link
            
            # print to outputfile:
            if w == True:
                outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(
                    hit, l[1], l[3], l[5], l[6], l[7], l[8], l[9], l[10]
                    ))
            else:
                outfile.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(
                    hit, l[1], l[3], l[5], l[6], l[7]
                    ))
    
    # NCBI complete has gi| in description:
    #
    # ncbi assembly:
    # "ftp://ftp.ncbi.nih.gov/genomes/all/" + a_id
    #
    # ENSEMBL:
    # "ftp://ftp.ensemblgenomes.org/pub/fungi/current/fasta/"
    
    infile.close()
    outfile.close()

if __name__ == "__main__":
    createTSV(infile, database, outfile)
