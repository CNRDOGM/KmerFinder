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
import sys, time
import os
from argparse import ArgumentParser
from operator import itemgetter
import re

if sys.version_info < (3, 0):
    import cPickle as pickle
else:
    import pickle

def makeorganismDB():
    ############################################################
    # Parse command line options
    ############################################################
    parser = ArgumentParser()
    parser.add_argument("-i", "--inputfile", dest="inputfilename", help="read from INFILE", metavar="INFILE")
    parser.add_argument("-o", "--outputfile", dest="outputfilename", help="write to OUTFILE", metavar="OUTFILE")
    parser.add_argument("-t", "--templatefile", dest="templatefilename", help="add to database TEMFILE", metavar="TEMFILE")
    #parser.add_argument("-p", "--pickleoutput", dest="pickleoutput",action="store_true", help="use pickle output")
    args = parser.parse_args()
    #
    # Open file for input sequence with kmers to save in database
    #
    if args.inputfilename != None:
      if args.inputfilename == "--":
        inputfile = sys.stdin
      else:
        inputfile = open(args.inputfilename,"r")


    #
    # open templatefile (already existing databse)
    #
    if args.templatefilename != None:
      templatefile = open(args.templatefilename+".p", "rb" )
      templatefile_lengths = open(args.templatefilename+".len.p", "rb" )
      try:
        templatefile_ulengths = open(args.templatefilename+".ulen.p", "rb" )
      except:
        # do nothing
        two=2
      templatefile_descriptions = open(args.templatefilename+".desc.p", "rb" )

    #
    # Harcode this to be true so I do not need to use the -p option
    #
    args.pickleoutput = True
    if args.outputfilename != None:
      if args.pickleoutput == True:
        outputfile = open(args.outputfilename+".p", "wb")
        outputfile_lengths = open(args.outputfilename+".len.p", "wb")
        outputfile_ulengths = open(args.outputfilename+".ulen.p", "wb")
        outputfile_descriptions = open(args.outputfilename+".desc.p", "wb")
      else:
        outputfile = open(args.outputfilename,"w")
    else:
      outputfile = sys.stdout



    #################################################################################
    #	Read inputfile
    #################################################################################

    org={}
    new_descriptions={}

    for l in inputfile:
        l=l.strip()
        l=l.split("\t")

        org[l[0]]=l[1]

        # get descriptions:
        if l[1] in new_descriptions:
          new_descriptions[l[1]] = new_descriptions[l[1]]+", "+ l[0]
        else:
          new_descriptions[l[1]] = l[0]

    #################################################################################
    # Read Template file
    #################################################################################

    sys.stdout.write("%s\n" % ("# Reading database of templates"))
    inputs = pickle.load(templatefile)
    lengths = pickle.load(templatefile_lengths)
    ulengths = pickle.load(templatefile_ulengths)
    descriptions = pickle.load(templatefile_descriptions)


    #################################################################################
    #	Change kmer-entries and counts
    #################################################################################

    sys.stdout.write("%s\n" % ("# Updating DB"))

    new_inputs={}
    new_lengths={}
    new_ulengths={}
    
    if sys.version_info < (3, 0):
        iter_inputs = inputs.iteritems()
    else:
        iter_inputs = inputs.items()
    for key, value in iter_inputs:
      # replace acc with organism name:
      value=value.split(",")
      #for i in range(0,len(value)):
      #  value[i]=org[value[i]]
      #
      #  # update counts:
      #  if value[i] in new_lengths:
      #    new_lengths[value[i]] = new_lengths[value[i]] + 1
      #  else:
      #    new_lengths[value[i]] = 1
      count = 0

      for i in value:

        value[count]=org[i]
        count = count +1

        # update counts:
        if org[i] in new_lengths:
          new_lengths[org[i]] = new_lengths[org[i]] + 1
        else:
          new_lengths[org[i]] = 1

      # update unique counts:
      uvalue= list(set(value))
      for i in uvalue:
        if i in new_ulengths:
          new_ulengths[i] = new_ulengths[i] + 1
        else:
          new_ulengths[i] = 1

      # update the kmer-dictionnary
      nvalue = ','.join(value)
      new_inputs[key]=nvalue


    ################################################################################
    #	Print new database
    ################################################################################


    pickle.dump(new_inputs, outputfile,2)
    pickle.dump(new_lengths, outputfile_lengths,2)
    pickle.dump(new_ulengths, outputfile_ulengths,2)
    pickle.dump(new_descriptions, outputfile_descriptions,2)

if __name__ == '__main__':
    makeorganismDB()
