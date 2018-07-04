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

def readPrintKmerDB():
    ############################################################
    # Parse command line options
    ############################################################

    parser = ArgumentParser()
    parser.add_argument("-t", "--templatefile", dest="templatefilename", help="add to database TEMFILE", metavar="TEMFILE")
    parser.add_argument("-u", "--ulengths", dest="print_ulengths", action="store_true", help="print amount of unique kmers found in each template sequence")
    parser.add_argument("-l", "--lengths", dest="print_lengths", action="store_true", help="print amount of kmers found in each template sequence")
    parser.add_argument("-i", "--inputs", dest="print_inputs", action="store_true", help="print hash table of kmers")
    parser.add_argument("-d", "--descriptions", dest="print_descriptions", action="store_true", help="print descriptions for each template sequence")
    args = parser.parse_args()


    #################################################################################
    # open templatefile (already existing databse)
    #################################################################################

    if args.templatefilename != None:

      if args.print_ulengths == True:
        templatefile_ulengths = open(args.templatefilename+".ulen.p", "rb" )

      if args.print_lengths ==True:
        templatefile_lengths = open(args.templatefilename+".len.p", "rb" )

      if args.print_inputs == True:
        templatefile = open(args.templatefilename+".p", "rb" )

      if args.print_descriptions == True:
        templatefile_descriptions = open(args.templatefilename+".desc.p", "rb" )
    else:
      sys.stderr.write("Please specify database!\n")
      sys.exit(2)


    #################################################################################
    # Read Template file
    #################################################################################

    sys.stdout.write("%s\n" % ("# Reading database of templates"))


    if args.print_ulengths == True:
      ulengths = pickle.load(templatefile_ulengths)

    if args.print_lengths ==True:
      lengths = pickle.load(templatefile_lengths)

    if args.print_inputs == True:
      inputs = pickle.load(templatefile)

    if args.print_descriptions == True:
      descriptions = pickle.load(templatefile_descriptions)


    #################################################################################
    #	Change kmer-entries and counts
    #################################################################################

    sys.stdout.write("%s\n" % ("# Printing DB entries:"))

    if args.print_ulengths == True:
      print("ulengths:")
      for key, value in ulengths.iteritems():
        print(key + " " + str(value))

    if args.print_lengths ==True:
      print("lengths:")
      for key, value in lengths.iteritems():
        print(key + " " + str(value))

    if args.print_inputs == True:
      print("hash table:")
      for key, value in inputs.iteritems():
        print(key + " " + str(value))

    if args.print_descriptions == True:
      print("descriptions:")
      for key, value in descriptions.iteritems():
        print(key + " " + str(value))

if __name__ == '__main__':
    readPrintKmerDB()
