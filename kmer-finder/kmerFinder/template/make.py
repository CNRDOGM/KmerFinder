#!/usr/bin/env python3

# Copyright (c) 2014, Ole Lund, Technical University of Denmark
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#
# Import libraries
#
import sys
import time
import os
from argparse import ArgumentParser
from operator import itemgetter
import re
from string import maketrans

if sys.version_info < (3, 0):
    import cPickle as pickle
else:
    import pickle

#################################################################
# FUNCTIONS:
#################################################################


# --------------------------------------
# reverse complement sequence:
# --------------------------------------
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


# ---------------------------------------------------------
# check homology
# ---------------------------------------------------------

def check_homology(inputseqsegments, inputs):
    global etta

    # Make list of unique k-mers in entry
    queryindex = {}
    uquerymers = 0

    for seq in inputseqsegments:
        seqlen = len(seq)
        for qseq in[seq, reversecomplement(seq)]:
            for j in range(0, seqlen - kmersize + 1):
                submer = qseq[j:j + kmersize]
                if prefix == qseq[j:j + prefixlen]:
                    if submer in queryindex:
                        queryindex[submer] += 1
                    else:
                        queryindex[submer] = 1
                        uquerymers += 1

    # Search for matches:
    mincoverage = 1
    templateentries = {}

    for submer in queryindex:
        if submer in inputs:
            if queryindex[submer] >= mincoverage:
                matches = inputs[submer].split(",")
                # get unique list:
                umatches = list(set(matches))
                for match in umatches:
                    # count matches:
                    if match in templateentries:
                        templateentries[match] += 1
                    else:
                        templateentries[match] = 1

    # calculate frac_q:
    frac_q = 0.0
    hitname = ""
    score = 0
    sortedlist = sorted(
        templateentries.items(), key=itemgetter(1), reverse=True)
    for template, score in sortedlist:
        # multiplication by 2 because of reverse complement
        frac_q = score / (float(uquerymers) + etta)
        hitname = template
        break

    del templateentries
    del queryindex

    return (hitname, frac_q, score)

# -----------------------------------------------------------
# update database:
# -----------------------------------------------------------


def update_database(inputseqsegments, inputname, filters):

    # define global variables:
    global inputs, lengths, ulengths, descriptions, desc
    global Nstored, Nstored_old, Nustored, Nustored_old
    global kmer_count, kmersize, prefixlen, stepsize, t0, t1, printfreq
    global filterfilename, homthres

    # Start of database update
    for s in inputseqsegments:
        for seq in[s, reversecomplement(s)]:
            start = 0
            while start < len(seq) - kmersize + 1:
                submer = seq[start:start + kmersize]

                if prefix == seq[start:start + prefixlen]:
                    if (filterfilename is not None and submer not in filters
                       or filterfilename is None):
                        Nstored += 1
                        if submer in inputs:
                            if (inputs[submer].find(inputname) == -1):
                                Nustored += 1
                            inputs[submer] = inputs[submer] + "," + inputname
                        else:
                            inputs[submer] = inputname
                            Nustored += 1
                kmer_count += 1
                start += stepsize
    # update nr of kmers:
    lengths[inputname] = Nstored - Nstored_old
    # update nr of unique kmers:
    ulengths[inputname] = Nustored - Nustored_old
    # update descriptions:
    descriptions[inputname] = desc
    # print inputname,lengths[inputname], Nstored, Nstored_old,"i: ",i,
    # len(inputseq)
    Nstored_old = Nstored
    Nustored_old = Nustored


# ------------------------------------------------------------
# process entry
# ------------------------------------------------------------
def process_entry(old_inputname):

    global inputseqsegments, inputs, homthres, filters

    sys.stdout.write("%s %s\n" % ("# Entry read", old_inputname))

    # check homology:
    if homthres is not None:
        sys.stdout.write("%s\n" % ("# Checking for homology"))

        (hitname, frac_q, score) = check_homology(inputseqsegments, inputs)
        sys.stdout.write(
            "# Max frac_q similarity of %s to %s frac_q: %s Score: %s\n" %
            (old_inputname, hitname, frac_q, score))

    # decide to exclude / include entry:
    if (homthres is not None and frac_q >= homthres):
        # exclude entry:
        sys.stdout.write(
            "# Skipping entry: %s in database due to similarity"
            "to %s frac_q: %s\n" % (old_inputname, hitname, frac_q)
        )
    if (homthres is None or (homthres is not None and frac_q < homthres)):
        # include entry -> update database:
        sys.stdout.write("%s %s\n" % ("# Including entry: ", old_inputname))
        update_database(inputseqsegments, old_inputname, filters)

def makeTemplateDB():
    # #########################################################################
    # PARSE COMMAND LINE OPTIONS:
    # #########################################################################

    parser = ArgumentParser()
    parser.add_argument("-i", "--inputfile", dest="inputfilename",
                      help="read from INFILE", metavar="INFILE")
    parser.add_argument("-l", "--inputfilelist", dest="inputfilelist",
                      help="read a list of fatsa file locations",
                      metavar="INFILELIST")
    parser.add_argument("-f", "--filterfile", dest="filterfilename",
                      help="filter (ignore) K-mers present in FILTERFILE",
                      metavar="FILTERFILE")
    parser.add_argument("-o", "--outputfile", dest="outputfilename",
                      help="write to OUTFILE", metavar="OUTFILE")
    parser.add_argument("-k", "--kmersize", dest="kmersize", help="Size of KMER",
                      metavar="KMERSIZE")
    parser.add_argument("-t", "--homthres", dest="homthres",
                      help="Threshold for homology reduction", metavar="HOMTHRES")
    parser.add_argument("-s", "--stepsize", dest="stepsize",
                      help="Size of step between K-mers", metavar="STEPSIZE")
    parser.add_argument("-x", "--prefix", dest="prefix", help="type of prefix",
                      metavar="PREFIX")
    parser.add_argument("-a", "--templatefile", dest="templatefilename",
                      help="add to database TEMFILE", metavar="TEMFILE")
    parser.add_argument("-c", "--organismlist", dest="organismlistname",
                      help="provide organism list to replace IDs ORGLIST",
                      metavar="ORGLIST")
    args = parser.parse_args()

    ##########################################################################
    # DEFINE GLOBAL VARIABLES:
    ##########################################################################

    global inputs, lengths, ulengths, descriptions, desc, inputseqsegments
    global Nstored, Nstored_old, Nustored, Nustored_old
    global kmer_count, kmersize, prefix, prefixlen, stepsize, t0, t1, printfreq
    global filterfilename, homthres, filters, organismlist, etta

    # Open file for input sequence with kmers to save in database:
    if args.inputfilename is not None:
        if args.inputfilename == "--":
            inputfile = sys.stdin
        else:
            # inputfile = open(args.inputfilename,"r")
            inputfile = args.inputfilename

    # Open list of FASTA file locations:
    if args.inputfilelist is not None:

        if args.inputfilelist == "--":
            inputfilelist = sys.stdin
        else:
            inputfilelist = open(args.inputfilelist, "r")
    elif inputfile != "":
        inputfilelist = [inputfile]

    # Open file to filter on (kmers not to save in database):
    if args.filterfilename is not None:
        filterfile = open(args.filterfilename, "r")
        filterfilename = optins.filterfilename
    else:
        filterfilename = None

    # Open templatefile to add to already existing databse:
    if args.templatefilename is not None:
        templatefile = open(args.templatefilename + ".p", "rb")
        templatefile_lengths = open(args.templatefilename + ".len.p", "rb")
        try:
            templatefile_ulengths = open(
                args.templatefilename + ".ulen.p", "rb")
        except:
            # do nothing
            two = 2
            templatefile_ulengths = None
            templatefile_descriptions = open(
                args.templatefilename + ".desc.p", "rb")


    # Harcode this to be true so I do not need to use the -p option:

    if args.outputfilename is not None:
        outputfile = open(args.outputfilename + ".p", "wb")
        outputfile_lengths = open(args.outputfilename + ".len.p", "wb")
        outputfile_ulengths = open(args.outputfilename + ".ulen.p", "wb")
        outputfile_descriptions = open(args.outputfilename + ".desc.p", "wb")
    else:
        outputfile = sys.stdout

    # get kmer-size:
    if args.kmersize is not None:
        kmersize = int(args.kmersize)
    else:
        kmersize = 16


    # get size of step when looking for K-mers in the sequence:
    if args.stepsize is not None:
        stepsize = int(args.stepsize)
    else:
        stepsize = 1

    # Homology threshold for when to include entry in database:
    if args.homthres is not None:
        homthres = float(args.homthres)
    else:
        homthres = None

    # Prefix to use fro filtering sequences:
    if args.prefix is not None:
        prefix = args.prefix
        prefixlist = [prefix]
        prefixlen = len(prefixlist[0])
    else:
        prefix = ''
        prefixlist = [prefix]
        prefixlen = len(prefixlist[0])


    # get organism list:
    if args.organismlistname is not None:
        organismlist = open(args.organismlistname, "r")
    else:
        organismlist = None

    ##################################################################
    # INITIALIZE STATISTICS
    ##################################################################

    # # of kmers
    kmer_count = 0
    # Start time to keep track of progress
    t0 = time.time()
    # Print progress
    printfreq = 100000
    # frequenct to save sorted list to db
    dbsavefreq = 30000000

    etta = 0.0001

    ###################################################################
    # READ SEQUENCES FROM FILTERFILE AND SAVE KMERS
    ###################################################################

    filterseqsegments = []
    filters = {}
    Nfilters = 0
    i = 0
    t1 = time.time()
    if args.filterfilename is not None:
        sys.stdout.write("%s\n" % ("# Reading filterfile"))
        for line in filterfile:
            line = line.rstrip('\n')
            fields = line.split()
            if len(line) >= 1:
                if fields[0][0] == ">":
                    if (i > 0):
                        #
                        # Fasta entry read
                        #
                        filterseq = ''.join(filterseqsegments)
                        for seq in [filterseq, reversecomplement(filterseq)]:
                            start = 0
                            while start < len(seq) - kmersize:
                                submer = seq[start:start + kmersize]
                                if prefix == seq[start:start + prefixlen]:
                                    if submer not in filters:
                                        filters[submer] = filtername
                                kmer_count += 1
                                if kmer_count % printfreq == 0:
                                    t1 = time.time()
                                    sys.stdout.write(
                                        "\r%s kmers (%s kmers / s)" %
                                        ("{:,}".format(kmer_count),
                                         "{:,}".format(kmer_count / (t1 - t0)))
                                    )
                                    sys.stdout.flush()
                                start += stepsize
                    del filterseqsegments
                    filterseqsegments = []
                    i = 0
                    filterseq = ""
                    filtername = fields[0][1:]
                else:
                    filterseqsegments.append("")
                    filterseqsegments[i] = fields[0]
                    i += 1
        filterseq = ''.join(filterseqsegments)
        for seq in [filterseq, reversecomplement(filterseq)]:
            start = 0
            while start < len(seq) - kmersize:
                submer = seq[start:start + kmersize]
                if prefix == seq[start:start + prefixlen]:
                    if submer not in filters:
                        filters[submer] = filtername
                kmer_count += 1
                if kmer_count % printfreq == 0:
                    t1 = time.time()
                    sys.stdout.write(
                        "\r%s kmers (%s kmers / s)" %
                        ("{:,}".format(kmer_count),
                            "{:,}".format(kmer_count / (t1 - t0)))
                    )
                    sys.stdout.flush()
                start += stepsize
        #
        # Print final statistics for filterfile
        #
        t1 = time.time()
        sys.stdout.write(
            "\r%s kmers (%s kmers / s)" %
            ("{:,}".format(kmer_count), "{:,}".format(kmer_count / (t1 - t0)))
        )
        sys.stdout.flush()
        # sys.stdout.write("\n")
    del filterseqsegments


    ##########################################################################
    # READ TEMPLATEFILE
    ##########################################################################

    inputs = {}
    lengths = {}
    ulengths = {}
    descriptions = {}

    Nstored = 0
    Nstored_old = Nstored
    # count only one per sequence:
    Nustored = 0
    Nustored_old = Nustored

    if args.templatefilename is not None:
        sys.stdout.write("%s\n" % ("# Reading database of templates"))
        inputs = pickle.load(templatefile)
        lengths = pickle.load(templatefile_lengths)
        try:
            ulengths = pickle.load(templatefile_ulengths)
        except:
            sys.stderr.write('No ulen.p file found for database')
            SystemExit()
        descriptions = pickle.load(templatefile_descriptions)

        # Count number of k-mers and number of unique k-mers:
        for name in lengths:
            Nstored += 1
            if name in ulengths:
                Nustored += 1
        Nstored_old = Nstored
        Nustored_old = Nustored


    ##########################################################################
    # READ ORGANISM LIST
    ##########################################################################


    if args.organismlistname is not None:
        sys.stdout.write("%s\n" % ("# Reading organism list"))

        organism = {}

        for l in organismlist:
            l = l.strip()
            l = l.split("\t")
            organism[l[0]] = l[1]

    #################   ######################################################
    #   PROCESS A LIST OF FASTA FILE LOCATIONS
    #######################################################################

    t1 = time.time()

    sys.stdout.write("%s\n" % ("# Reading inputfile(s)"))

    kmer_count = 0
    inputseqsegments = []
    inputseqsegments.append("")
    i = 0
    desc = ''
    old_inputname = ''
    original_inputname = ''


    for l in inputfilelist:
        l = l.strip()
        fastafile = open(l, "r")
        for line in fastafile:
            line = line.rstrip('\n')
            fields = line.split()
            if len(line) >= 1:
                if fields[0][0] == ">":

                    # get input name and translate if necessary:
                    inputname = fields[0][1:]
                    if organismlist is not None:
                        original_inputname = inputname
                        inputname = organism[inputname]

                    # new entry:
                    if inputname != old_inputname:
                        if old_inputname != '':

                            # process entry (homology check and include in
                            # database):
                            process_entry(old_inputname)

                            # prepare sequence variable:
                            del inputseqsegments
                            inputseqsegments = []
                            i = 0
                            inputseqsegments.append("")

                        # update old_inputname:
                        old_inputname = inputname

                        # store description:
                        if organismlist is not None:
                            desc = original_inputname
                        else:
                            desc = ' '.join(fields[1:len(fields)])
                        kmer_count_old = kmer_count

                    elif inputname == old_inputname:
                        # append description:
                        if organismlist is not None:
                            desc = desc + ", " + original_inputname
                        else:
                            desc = desc + ", " + ' '.join(fields[1:len(fields)])

                        # prepare sequnece variable:
                        inputseqsegments.append("")
                        i = i + 1

                else:
                    # read sequnece:
                    inputseqsegments[i] = inputseqsegments[i] + fields[0].upper()

        fastafile.close()

    # process last entry:
    process_entry(inputname)

    del inputseqsegments

    ############################################################
    # PRINT DATABASE OF KMERS
    ############################################################


    pickle.dump(inputs, outputfile, 2)
    pickle.dump(lengths, outputfile_lengths, 2)
    pickle.dump(ulengths, outputfile_ulengths, 2)
    pickle.dump(descriptions, outputfile_descriptions, 2)

    ###########################################################
    # PRINT FINAL STATISTICS
    ###########################################################
    t2 = time.time()
    sys.stdout.write(
        "\r%s kmers (%s kmers / s)\n" %
        ("{:,}".format(kmer_count), "{:,}".format(kmer_count / (t2 - t1)))
    )
    sys.stdout.flush()
    sys.stdout.write("\n")
    sys.stdout.write("# Total time used: %s s\n" % (t2 - t0))
    #
    # Done
    #

if __name__ == '__main__':
    makeTemplateDB()
