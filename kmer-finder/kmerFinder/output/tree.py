#!/usr/bin/env python3
#
# Import libraries
#
import sys
import time
import os
from math import sqrt
from argsparse import ArgumentParser
from operator import itemgetter
import re
#import bsddb3 as bsddb
import numpy
#import anydbm
#import shelve

if sys.version_info < (3, 0):
    import cPickle as pickle
else:
    import pickle
    xrange = range

#
# Functions
#
# Construct the reverse complemant from a sequence
#


def reversecomplement_old(seq):
    '''Reverse complement'''
    comp = ''
    for s in seq:
        if s == 'A':
            comp = comp + 'T'
        elif s == 'T':
            comp = comp + 'A'
        elif s == 'C':
            comp = comp + 'G'
        elif s == 'G':
            comp = comp + 'C'
        else:
            comp = comp + s
    return comp[::-1]


def reversecomplement(seq):
    '''Reverse complement'''
    return seq.translate(maketrans("ATGC", "TACG"))[::-1]


def makeTree():
    #
    # Parse command line options
    #
    parser = ArgumentParser()
    parser.add_argument("-t", "--templatefile", dest="templatefilename",
                        help="read from TEMFILE", metavar="TEMFILE")
    parser.add_argument(
        "-o", "--outputfile", dest="outputfilename", help="write to OUTFILE", metavar="OUTFILE")
    parser.add_argument(
        "-x", "--prefix", dest="prefix", help="prefix", metavar="_id")
    parser.add_argument("-a", "--printall", dest="printall", action="store_true",
                        help="Print matches to all templates in templatefile unsorted")
    parser.add_argument("-p", "--pickleinput", dest="pickleinput",
                        action="store_true", help="use pickle input")
    parser.add_argument("-d", "--dbinput", dest="dbinput",
                        action="store_true", help="use database input")
    parser.add_argument("-c", "--columns", dest="columns",
                        action="store_true", help="write output in column format")
    args = parser.parse_args()
    #
    # set up prefix filtering
    #
    if args.prefix != None:
        prefix = args.prefix
    else:
        prefix = ''
    prefixlen = len(prefix)
    #
    # Open files
    #
    t0 = time.time()
    #
    if args.templatefilename != None:
        if args.pickleinput == True:
            templatefile = open(args.templatefilename, "rb")
        # elif args.dbinput == True:
        #  #templates = anydbm.open(args.templatefilename, 'r')
        #  templates = bsddb.db.DB()
        #  templates.set_cachesize(1, 0)  # 1Gb, 0b
        #  dbtype = bsddb.db.DB_BTREE
        #  # or
        #  # dbtype = bsddb.db.DB_HASH
        #  templates.open(args.templatefilename, None, flags = bsddb.db.DB_RDONLY)
        else:
            templatefile = open(args.templatefilename, "r")
    else:
        sys.exit("No template file specified")
    #
    if args.outputfilename != None:
        outputfile = open(args.outputfilename, "w")
    else:  # If no output filename choose the same as the input filename
        outputfilename = os.path.splitext(args.inputfilename)[0]
        outputfile = open(outputfilename, "w")
    #
    # Make database of nmers
    #
    if not args.dbinput == True:
        templates = {}
    templateentries = {}
    templatebackground = {}
    templatebackgroundtot = 0
    Ntemplates = 0
    seq = ""
    name = "None"
    oligolen = 16
    #
    # Read Template file
    #
    # outputfile.write("%s\n" % ("# Reading database of templates"))
    if args.pickleinput == True:
        templates = pickle.load(templatefile)
        #templates = shelve.open(args.templatefilename, flag='c', protocol=2, writeback=False)
    elif not args.dbinput == True:
        for line in templatefile:
            fields = line.split()
            templates[fields[1]] = fields[2]
            #
            # Find number of Nmers in each template, and in total
            #
            matches = fields[2].split(",")
            for match in matches:
                if match in templatebackground:
                    templatebackground[match] += 1
                else:
                    templatebackground[match] = 1
                    Ntemplates += 1
                templatebackgroundtot += 1
                # print templatebackgroundtot
        oligolen = len(fields[1])
    #
    # Search for matches
    #
    # outputfile.write("%s\n" % ("# making distance matrix"))
    Nhits = 0
    numbers = {}
    names = {}
    #nkmers = {}
    nnumbers = 0
    for submer in templates:
            # print submer
        matches = templates[submer].split(",")
        for match in matches:
            #      if match in nkmers:
            #        nkmers[match] +=1
            #      else:
            #        nkmers[match] = 1
            if not match in numbers:
                numbers[match] = nnumbers
        names[nnumbers] = match
        nnumbers += 1
    #
    #
    #
    # for num in names:
    #  print num, names[num]

    #mat = [[0]*nnumbers]*nnumbers
    mat = numpy.zeros((nnumbers, nnumbers))
    # print mat
    #nkmers = [0]*nnumbers
    nkmers = numpy.zeros((nnumbers))
    # max=0
    for submer in templates:
        # print submer
        matches = templates[submer].split(",")
        # print submer
        for match1 in matches:
            # print match1
            for match2 in matches:
                # print match2
                mat[numbers[match1]][numbers[match2]] += 1
        if match1 == match2:
            nkmers[numbers[match1]] += 1
        # print submer, match1, match2, mat[numbers[match1]][numbers[match2]], mat[numbers[match2]][numbers[match1]], nkmers[numbers[match1]], nkmers[numbers[match2]]
            #mat[numbers[match2]][numbers[match1]] +=1
        # if mat[numbers[match2]][numbers[match1]] > max:
        #  max = mat[numbers[match2]][numbers[match1]]
    #
    #
    #

    # print mat

    # print nkmers
    #
    # for nam in numbers:
    #  print nam, numbers[nam], nkmers[numbers[nam]]
        # print nam, numbers[nam]
    #
    # write in column format
    #
    if args.columns == True:
        for i in xrange(0, nnumbers):
            mynamei = names[i][0:11]
            for j in xrange(0, nnumbers):
                mynamej = names[j][0:11]
                #lmax = max(nkmers[names[i]],nkmers[names[j]])
                #lmax = nkmers[names[i]]
                if nkmers[i] > nkmers[j]:
                    lmax = nkmers[i]
                else:
                    lmax = nkmers[j]
                outputfile.write("%s %s %s %s %s %s\n" % (
                    mynamei, mynamej, mat[i][j], nkmers[i], nkmers[j], lmax))
    #
    # Write in neighbor format
    #
    else:
        outputfile.write("%s\n" % (nnumbers))
        for i in xrange(0, nnumbers):
            mynamei = names[i][0:10]
            outputfile.write("%-10s " % (mynamei))
            for j in xrange(0, nnumbers):
                mynamej = names[i][0:10]
                #lmax = max(nkmers[names[i]],nkmers[names[j]])
                #lmax = nkmers[names[i]]
                if nkmers[i] > nkmers[j]:
                    lmax = nkmers[i]
                else:
                    lmax = nkmers[j]
                outputfile.write("%0.8f" % (1.0 - mat[i][j] / lmax))
                if ((j + 1) % 6 == 0) or (j == nnumbers - 1):
                    outputfile.write("\n")
                else:
                    outputfile.write(" ")
    #
    # Close files
    #
    #t1 = time.time()
    # sys.stdout.write("\r# %s entries (%s entries / s). Total time used %s seq" % ("{:,}".format(nnumbers), "{:,}".format(nnumbers / (t1-t0)),int(t1-t0)))
    # sys.stdout.flush()
    # sys.stdout.write("\n")
    # print "# Closing files"
    # inputfile.close()
    # templatefile.close()
    # outputfile.close()

if __name__ == '__main__':
    makeTree()
