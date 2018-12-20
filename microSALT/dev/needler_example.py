import pdb
import re
import subprocess
import os


class Meta_Set():
  def __init__(self):
    self.ref = "/mnt/hds/proj/bioinfo/MICROBIAL/MLST/references/gene_set/escherichia_coli.fna"
    self.found = list()

  def init(self):
    #Load dictionary
    dropped = 0
    with open(self.ref, 'r') as infile:
      lastkey = ""
      codons = ['ATG', 'TAG', 'TGA', 'TAA']
      rv_codons = ['TAC', 'ATC', 'ACT', 'ATT']
      for line in infile:
        line = line.rstrip()
        if '>' in line:
          lastkey = line
        #Length and codon filters
        elif len(line) >= 50 and (line[0:3] in codons and line[-3:] in codons) or (line[0:3] in rv_codons and line[-3:] in rv_codons):
          self.found.append(dict())
          self.found[-1] = (lastkey , line)
        else:
          #print("Dropped sequence {}. Length: {}, Start: {}, Stop: {}".format(lastkey.split(' ')[0], len(line), line[0:3], line[-3:]))
          dropped += 1

    print("Sequences initially: {} items".format(dropped+len(self.found)))
    print("Sequences remaining after length and codon filter: {} items".format(len(self.found)))

  def filter_codon(self):
    dropped = 0
    for line in self.found:
      lastkey = ""
      codons = ['ATG', 'TAG', 'TGA', 'TAA']
      rv_codons = ['TAC', 'ATC', 'ACT', 'ATT']
      #Length and codon filters
      if len(line[1]) >= 50 and (line[1][0:3] in codons and line[1][-3:] in codons) or (line[1][0:3] in rv_codons and line[1][-3:] in rv_codons):
        self.found.append(dict())
        self.found[-1] = (lastkey , line)
      else:
        dropped += 1

    print("Sequences remaining after length and codon filter: {} items".format(len(self.found)))

  def filter_id(self):  
    #Apply identity criteria
    found = self.found
    refined = list()
    while len(found) > 1:
      a = found.pop(0)
      with open('a', 'w') as temporary:
        temporary.write("{}\n{}".format(a[0],a[1]))
      with open('b', 'w') as temporary:
        for item in found:
          temporary.write("{}\n{}\n".format(item[0], item[1]))
      cmd = "water -asequence a -bsequence b -gapopen 10.0 -gapextend 0.5 -outfile out.water -datafile EDNAFULL"
      subprocess.check_output(cmd.split(), stderr=subprocess.DEVNULL)
      out = open("{}/out.water".format(os.getcwd()), 'r')
      content = out.readlines()
      #Read output
      similar = False
      for outl in content:
        if 'Length' in outl:
          length = int(re.search('\d+', outl).group(0))
        elif 'Identity' in outl:
          id = float(re.search('\(\ *(\d+.\d+)%\)', outl).group(1))
          #Identity filter when id found (always after length)
          if id >= 90 and length >= 100:
            similar = True
            break
          #Reset
          else:
            length = 0
            id = 0
      if not similar:
        refined.append(dict())
        refined[-1]=(a[0], a[1])

    print("Sequences after similarity filter: {} items".format(len(refined)))
    self.found = refined

  def filter_corner(self):
    #Apply corner criteria
    #Unnecessary due to codon filter. Use water and base on filter_id if considering.


metaset=Meta_Set()
metaset.init()
metaset.filter_id()
#metaset.filter_corner()
