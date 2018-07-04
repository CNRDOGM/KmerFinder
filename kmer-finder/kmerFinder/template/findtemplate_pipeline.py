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
from math import sqrt, pow
from argparse import ArgumentParser
from operator import itemgetter
import re

if sys.version_info < (3, 0):
    import cPickle as pickle
else:
    import pickle
#
# Functions

#
# Conservative two sided p-value from z-score
#
def fastp(z):
  '''Conservative two sided p-value from z-score'''
  if   z > 10.7016: p =  1e-26
  elif z > 10.4862: p =  1e-25
  elif z > 10.2663: p =  1e-24
  elif z > 10.0416: p =  1e-23
  elif z > 9.81197: p =  1e-22
  elif z > 9.5769:  p =  1e-21
  elif z > 9.33604: p =  1e-20
  elif z > 9.08895: p =  1e-19
  elif z > 8.83511: p =  1e-18
  elif z > 8.57394: p =  1e-17
  elif z > 8.30479: p =  1e-16
  elif z > 8.02686: p =  1e-15
  elif z > 7.73926: p =  1e-14
  elif z > 7.4409:  p =  1e-13
  elif z > 7.13051: p =  1e-12
  elif z > 6.8065:  p =  1e-11
  elif z > 6.46695: p =  1e-10
  elif z > 6.10941: p =  1e-9
  elif z > 5.73073: p =  1e-8
  elif z > 5.32672: p =  1e-7
  elif z > 4.89164: p =  1e-6
  elif z > 4.41717: p =  1e-5
  elif z > 3.89059: p =  1e-4
  elif z > 3.29053: p =  1e-3
  elif z > 2.57583: p =  0.01
  elif z > 1.95996: p =  0.05
  elif z > 1.64485: p =  0.1
  else            : p =  1.0
  return p
#
# Conservative two sided p-value from z-score
#
def z_from_two_samples(r1,n1,r2,n2):
  '''Comparison of two fractions, Statistical methods in medical research, Armitage et al. p. 125'''
  #
  # r1: positives in sample 1
  # n1: size of sample 1
  # r2: positives in sample 2
  # n2: size of sample 2
  #
  #r1=38
  #r2=102
  #n1=100
  #n2=200
  p1= float(r1)/(float(n1)+etta)
  p2= float(r2)/(float(n2)+etta)
  q1=1-p1
  q2=1-p2
  p=(float(r1)+r2)/(n1+n2+etta)
  q=1-p
  z=(p1-p2)/sqrt(p*q*(1/(n1+etta)+1/(n2+etta))+etta)
  #print r1,n1,r2,n2,p1,p2,q1,q2,p,q,z
  return z
#
#
#
def reversecomplement(seq):
    '''Reverse complement'''
    comp = ''
    for s in seq:
        if s == 'A': comp = comp + 'T'
        elif s == 'T': comp = comp + 'A'
        elif s == 'C': comp = comp + 'G'
        elif s == 'G': comp = comp + 'C'
        else: comp = comp + s
    return comp[::-1]

#
# Parse command line options
#
parser = ArgumentParser()
parser.add_argument("-i", "--inputfile", dest="inputfilename", help="read from INFILE", metavar="INFILE")
parser.add_argument("-t", "--templatefile", dest="templatefilename", help="read from TEMFILE", metavar="TEMFILE")
parser.add_argument("-o", "--outputfile", dest="outputfilename", help="write to OUTFILE", metavar="OUTFILE")
parser.add_argument("-k", "--kmersize", dest="kmersize", help="Size of k-mer, default 16", metavar="KMERSIZE")
parser.add_argument("-x", "--prefix", dest="prefix", help="prefix, e.g. ATGAC, default none", metavar="_id")
parser.add_argument("-a", "--printall", dest="printall",action="store_true", help="Print matches to all templates in templatefile unsorted")
parser.add_argument("-w", "--winnertakesitall", dest="wta",action="store_true", help="kmer hits are only assigned to most similar template")
parser.add_argument("-e", "--evalue", dest="evalue", help="Maximum E-value", metavar="EVALUE")
parser.add_argument("-r", "--round_output", dest="round_out",action="store_true", help="Print rounded output 1=round 0= do not round")
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
#
#
if args.evalue != None:
  evalue = float(args.evalue)
else:
  evalue = 0.05

#if args.round_out != None:
#  round_out = int(args.round_out)
#else:
#  round_out = 0


#
# Open files
#
t0 = time.time()
if args.inputfilename != None:
  if args.inputfilename == "--":
    inputfile = sys.stdin
  else:
    inputfile = open(args.inputfilename,"r")
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
else:
  sys.exit("No template file specified")
#
if args.outputfilename != None:
  outputfile = open(args.outputfilename,"w")
else: # If no output filename choose the same as the input filename
  outputfilename = os.path.splitext(args.inputfilename)[0]
  outputfile = open(outputfilename,"w")
#
# Size of K-mer
#
if args.kmersize != None:
  kmersize = int(args.kmersize)
else:
  kmersize = 16
#
# Make database of nmers
#
templates = {}
templateentries = {}
templateentries_tot = {}
templatebackground = {}
templatebackgroundtot = 0
Ntemplates =0
oligolen=kmersize
#
# Read Template file
#
sys.stderr.write("%s\n" % ("# Reading database of templates"))
templates = pickle.load(templatefile)
templates_lengths = pickle.load(templatefile_lengths)
try:
  templates_ulengths = pickle.load(templatefile_ulengths)
except:
   sys.stderr.write('No ulen.p file found for database')
   SystemExit()
templates_descriptions = pickle.load(templatefile_descriptions)
#
# Count number of k-mers, and number of unique k-mers
#
template_tot_len = 0
template_tot_ulen = 0
Ntemplates = 0
#length added
for name in templates_lengths:
  template_tot_len  += templates_lengths[name]
  template_tot_ulen += templates_ulengths[name]
  Ntemplates += 1

#print template_tot_len
#
# Read inputfile
#
#queryseq = []
queryseq = ""
queryseqsegments = []
#consensusseq = []
#queryname = []
#querydesc = []
#queryname = ""
#querydesc = ""
Nquerys=0
queryindex = {}
qtotlen=0
querymers=0
uquerymers=0
i=0
if args.inputfilename != None:
  sys.stderr.write("%s\n" % ("# Reading inputfile"))
  for line in inputfile:
    fields=line.split()
    if len(fields)>=1:
      if fields[0][0] == ">":
        if (i>0):
          #queryseq[-1] = ''.join(queryseqsegments)
          queryseq = ''.join(queryseqsegments)

          # Update dictionary of kmers
          #
          seqlen = len(queryseq)
          qtotlen += seqlen;
          #for qseq in [queryseq[i],reversecomplement(queryseq[i])]:
          for qseq in [queryseq]:
            for j in range(0, seqlen-oligolen+1):
              submer = qseq[j:j+oligolen]
              if prefix == qseq[j:j+prefixlen]:
                if submer in queryindex:
                  if submer in templates:
                    queryindex[submer] += 1
                  querymers += 1
                else:
                  if submer in templates:
                    queryindex[submer] = 1
                  querymers += 1
                  uquerymers += 1
        del queryseqsegments
        queryseqsegments = []
        i=0
        #queryseq.append("")
        #consensusseq.append("")
        #queryname.append(fields[0][1:])
        #querydesc.append(re.sub(r"^[^\s]+\s","",line.strip()))
      elif fields[0][0] == "@":
        # Fastq file
        if (i>0):
          #queryseq[-1] = ''.join(queryseqsegments)
          queryseq = ''.join(queryseqsegments)
          #
          # Update dictionary of kmers
          #
          seqlen = len(queryseq)
          qtotlen += seqlen;
          #for qseq in [queryseq[i],reversecomplement(queryseq[i])]:
          for qseq in [queryseq]:
            for j in range(0, seqlen-oligolen+1):
              submer = qseq[j:j+oligolen]
              if prefix == qseq[j:j+prefixlen]:
                if submer in queryindex:
                  if submer in templates:
                    queryindex[submer] += 1
                  querymers += 1
                else:
                  if submer in templates:
                    queryindex[submer] = 1
                  querymers += 1
                  uquerymers += 1
        del queryseqsegments
        queryseqsegments = []
        i=0
        #queryseq.append("")
        queryseq = ""
        #consensusseq.append("")
        #queryname.append(fields[0][1:])
        #querydesc.append(re.sub(r"^[^\s]+\s","",line.strip()))
        try:
          line = next(inputfile)
          fields=line.split()
          queryseqsegments.append("")
          queryseqsegments[i] = fields[0]
          i+=1
          line = next(inputfile)
          line = next(inputfile)
        except:
          break
      else:
        queryseqsegments.append("")
        queryseqsegments[i] = fields[0].upper()
        i+=1
  #queryseq[-1] = ''.join(queryseqsegments)
  queryseq = ''.join(queryseqsegments)
  #
  # Update dictionary of K-mers
  #
  seqlen = len(queryseq)
  qtotlen += seqlen;
  #for qseq in [queryseq[i],reversecomplement(queryseq[i])]:
  for qseq in [queryseq]:
    for j in range(0, seqlen-oligolen+1):
      submer = qseq[j:j+oligolen]
      if prefix == qseq[j:j+prefixlen]:
        if submer in queryindex:
          if submer in templates:
            queryindex[submer] += 1
          querymers += 1
        else:
          if submer in templates:
            queryindex[submer] = 1
          querymers += 1
          uquerymers += 1
del queryseqsegments
#
# Search for matches
#
sys.stderr.write("%s\n" % ("# Searching for matches of input in template"))
mincoverage = 1
Nhits=0
for submer in queryindex:
  if submer in templates:
    if queryindex[submer] >= mincoverage:
      matches = templates[submer].split(",")
      matches = list(set(matches))
      for match in matches:
        #
        # Update counts for the templates containing the k-mer
        #
        Nhits += 1
        if match in templateentries:
          templateentries[match] += 1
        else:
          templateentries[match] = 1
      #
      # Make list of unique matches (by converting to a dict and then back to a vector)
      # A match may occur more than once if the k-mer is found in more than one position
      # in that template
      #
      umatches = list(set(matches))
      for match in umatches:
        #
        # Add the number of times the k-mer is found in the input sequence to the
        # templateentries_tot dictionary
        #
        if match in templateentries_tot:
          templateentries_tot[match] += queryindex[submer]
        else:
          templateentries_tot[match] = queryindex[submer]
#
# Print best scoring entries sorted
#
minscore = 0
#
# etta is a small number to avoid division by zero
#
etta = 1.0e-8
sys.stderr.write("%s\n" % ("# Search statistics"))
sys.stderr.write("%s\n" % ("# Total number of hits: %s") % (Nhits))
sys.stderr.write("%s\n" % ("# Total number of kmers in templates : %s") % (template_tot_len))
sys.stderr.write("%s\n" % ("# Minimum number of k-mer hits to report template: %s") % (minscore))
sys.stderr.write("%s\n" % ("# Maximum multiple testing corrected E-value to report match : %s") % (evalue))
sys.stderr.write("%s\n" % ("# Printing best matches"))

if (args.round_out==True) and (args.wta!=True):
  outputfile.write("#Template\tScore\tExpected\tz\tp_value\tfrac_q\tfrac_d\tcoverage\tKmers in Template\tDescription\n")
elif (args.round_out==True) and (args.wta==True):
  outputfile.write("#Template\tScore\tExpected\tz\tp_value\tfrac_q\tfrac_d\tcoverage\ttotal frac_q\ttotal frac_d\ttotal coverage\tKmers in Template\tDescription\n")
elif (args.round_out!=True) and (args.wta!=True):
  outputfile.write("Template\tScore\tExpected\tz\tp_value\tfrac_q\tfrac_d\tcoverage\tunique_Kmers_in_Template\tDescription\n")
elif (args.round_out!=True) and (args.wta==True):
  outputfile.write("Template\tScore\tExpected\tz\tp_value\tfrac_q\tfrac_d\tcoverage\ttotal_frac_q\ttotal_frac_d\ttotal_coverage\tunique_Kmers_in_Template\tDescription\n")


#if args.pickleinput == True:
if not args.wta == True:
  #
  # Winner takes alle option is false, just do normal statistics
  #
  sortedlist= sorted(templateentries.items(), key = itemgetter(1), reverse=True)
  for template,score in sortedlist:
    if score > minscore:
      expected = float(Nhits)*float(templates_ulengths[template])/float(template_tot_ulen)
      #z = (score - expected)/sqrt(score + expected+etta)
      #p  = fastp(z)
      #
      # If expected < 1 the above poisson approximation is a poor model
      # Use instead: probabilyty of seing X hits is p**X if probability
      # of seing one hit is p (like tossing a dice X times)
      #
      #if expected <1:
      #  p = expected**score
      #
      # Comparison of two fractions, Statistical methods in medical research, Armitage et al. p. 125
      #
      z = z_from_two_samples(score,templates_ulengths[template],Nhits,template_tot_ulen)
      p  = fastp(z)
      #
      # Correction for multiple testing
      #
      p_corr = p*Ntemplates
      frac_q = score/(float(uquerymers)+etta)
      frac_d = score/(templates_ulengths[template]+etta)
      coverage = 2*templateentries_tot[template]/float(templates_lengths[template])
      if p_corr <= evalue:
        if args.round_out==True:
          outputfile.write("%-12s\t%8s\t%8s\t%8.1f\t%4.1e\t%4.1e\t%4.1e\t%4.1e\t%8d\t%s\n" %
            (template, score, int(round(expected)), round(z,1), p_corr, frac_q, frac_d, coverage, templates_ulengths[template], templates_descriptions[template].strip()))
        else:
          outputfile.write("%-12s\t%8s\t%8s\t%8.1f\t%4.1e\t%4.1e\t%4.1e\t%4.1e\t%8d\t%s\n" %
            (template, score, expected, z, p_corr, frac_q, frac_d, coverage, templates_ulengths[template], templates_descriptions[template].strip()))
else:
  #
  # New version of algorithm
  #
  w_templateentries = templateentries.copy()
  w_templateentries_tot = templateentries_tot.copy()
  w_Nhits = Nhits

  maxhits=100
  hitcounter=1
  stop = False
  while not stop:
    hitcounter += 1
    if hitcounter > maxhits:
      stop= True
    sortedlist= sorted(w_templateentries.items(), key = itemgetter(1), reverse=True)
    for template,score in sortedlist:
      if score > minscore:
        expected = float(w_Nhits)*float(templates_ulengths[template])/float(template_tot_ulen)
        #z = (score - expected)/sqrt(score + expected+etta)
        #p  = fastp(z)
        #
        # If expected < 1 the above poisson approximation is a poor model
        #
        #if expected <1:
        #  p = expected**score
        #
        # Comparison of two fractions, Statistical methods in medical research, Armitage et al. p. 125
        #
        z = z_from_two_samples(score,templates_ulengths[template],w_Nhits,template_tot_ulen)
        p  = fastp(z)
        #
        # correction for multiple testing
        #
        p_corr = p*Ntemplates
        #print score,float(uquerymers),etta
        frac_q = score/(float(uquerymers)+etta)
        frac_d = score/(templates_ulengths[template]+etta)
        coverage = 2*w_templateentries_tot[template]/float(templates_lengths[template])

        #calculate total values:
        tot_frac_q = templateentries[template]/(float(uquerymers)+etta)
        tot_frac_d = templateentries[template]/(templates_ulengths[template]+etta)
        tot_coverage = 2*templateentries_tot[template]/float(templates_lengths[template])

        if p_corr <= evalue:
          if args.round_out==True:
            outputfile.write("%-12s\t%8s\t%8s\t%8.1f\t%4.1e\t%4.1e\t%4.1e\t%4.1e\t%4.1e\t%4.1e\t%4.1e\t%8d\t%s\n" %
            (template, score, int(round(expected)), round(z,1), p_corr, frac_q, frac_d, coverage, tot_frac_q, tot_frac_d, tot_coverage, templates_ulengths[template], templates_descriptions[template].strip()))
          else:
            outputfile.write("%-12s\t%8s\t%8s\t%8.1f\t%4.1e\t%4.1e\t%4.1e\t%4.1e\t%4.1e\t%4.1e\t%4.1e\t%8d\t%s\n" %
              (template, score, expected, z, p_corr, frac_q, frac_d, coverage, tot_frac_q, tot_frac_d, tot_coverage, templates_ulengths[template], templates_descriptions[template].strip()))
          #
          # remove all kmers in best hit from queryindex
          #
          for submer in queryindex:
            if submer in templates:
              matches = templates[submer].split(",")
              #matches = list(set(templates[submer].split(",")))
              if template in matches:
                #querymers -= queryindex[submer]
                #uquerymers -= 1
                #if (uquerymers<0):
                #  uquerymers = 0
                queryindex[submer] = 0
          #
          # find best hit like before
          #
          w_Nhits=0
          del w_templateentries
          del w_templateentries_tot
          w_templateentries = {}
          w_templateentries_tot = {}
          for submer in queryindex:
            if submer in templates:
              if queryindex[submer] >= mincoverage:
                matches = templates[submer].split(",")
                matches = list(set(matches))
                for match in matches:
                  #
                  # Update counts for the templates containing the k-mer
                  #
                  w_Nhits += 1
                  if match in w_templateentries:
                    w_templateentries[match] += 1
                  else:
                    w_templateentries[match] = 1
                #
                # Make list of unique matches (by converting to a dict and then back to a vector)
                # A match may occur more than once if the k-mer is found in more than one position
                # in that template
                #
                umatches = list(set(matches))
                for match in umatches:
                  #
                  # Add the number of times the k-mer is found in the input sequence to the
                  # templateentries_tot dictionary
                  #
                  if match in w_templateentries_tot:
                    w_templateentries_tot[match] += queryindex[submer]
                  else:
                    w_templateentries_tot[match] = queryindex[submer]

        else:
          stop = True
          break
#
# Close files
#
t1 = time.time()
if args.round_out==True:
  sys.stderr.write("#info\t %s kmers (%s kmers / s). Total time used: %s sec" % ("{:,}".format(querymers), "{:,}".format(querymers / (t1-t0)),int(t1-t0)))
else:
  sys.stderr.write("\r# %s kmers (%s kmers / s). Total time used: %s sec" % ("{:,}".format(querymers), "{:,}".format(querymers / (t1-t0)),int(t1-t0)))
sys.stdout.flush()
sys.stderr.write("\n")
sys.stderr.write("DONE!")
sys.stderr.write("# Closing files\n")
