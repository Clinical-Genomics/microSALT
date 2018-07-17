"""Performs sample comparison analysis on a custom set of files
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import os

class Meta_Analyzer():

  def __init__(self, file_list, config, log):
    self.config = config
    self.logger = log
    self.listfile = file_list
    self.distancefile = 'cgMLST_distance.txt'

    self.files = list()
    self.names = list()

    with open (self.listfile, 'r') as infile:
      for line in infile:
        if line != '':
          self.files.append(line.rstrip())
    self.names = list()
    for entry in self.files:
      self.names.append(entry.split("/")[-3])


  def calc_dist(self):
    """Calculates distances between different typing profiles"""
    #Load a profile dictionary
    fp = dict()
    for sample in self.names:
      fp[sample] = dict()
      fn = self.files[self.names.index(sample)]
      with open(fn, 'r') as currfile:
        for line in currfile:
          line = line.rstrip().split(':')
          if line[0] in fp[sample]:
            self.logger.error("You dun goofed!")
          fp[sample][line[0]] = line[1]
    #Create a score dictionary
    score = dict()
    for sample1 in self.names:
      for locus, allele in fp[sample1].items():
        for sample2 in self.names:
          if not '{}-{}'.format(sample1, sample2) in score:
            score['{}-{}'.format(sample1, sample2)] = 0
          if locus in fp[sample2]:
            if allele != fp[sample2][locus]:
              score['{}-{}'.format(sample1, sample2)] += 1

    #Write distance file information
    f = open(self.distancefile, 'w+')
    header = ""
    for sample in self.names:
      header += "\t{}".format(sample)
    f.write("{}\n".format(header))
    for sample1 in self.names:
      row="{}".format(sample1)
      for sample2 in self.names:
        row +="\t{}".format(score['{}-{}'.format(sample1, sample2)]) 
      f.write("{}\n".format(row))
    f.close()
    self.logger.info("Distance file written to {}/{}".format(os.getcwd(), self.distancefile))
