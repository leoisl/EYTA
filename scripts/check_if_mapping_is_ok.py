#usage: python check_if_mapping_is_ok.py <output_folder> <k_value>
import glob
import sys
from Bio.Seq import Seq

def rc(read):
  my_seq = Seq(read)
  return str(my_seq.reverse_complement())

def parseLine(line):
  unitigInfos=[]
  splitLine = line.split()
  for i in xrange(0, len(splitLine), 3):
    unitigInfos.append((int(splitLine[i]), splitLine[i+1], int(splitLine[i+2])))
  return unitigInfos

def checkIfMappingIsOk(read, unitigInfos, unitigs, k):
  assembly=""
  for currentUnitig in unitigInfos:
    #get the correct kmer to add
    toAdd=unitigs[currentUnitig[0]]
    if currentUnitig[1]=='R':
      toAdd=rc(toAdd)
    toAdd = toAdd[currentUnitig[2]:currentUnitig[2]+k]

    if len(assembly)==0:
      assembly+=toAdd
    else:
      assembly+=toAdd[k-1:]

  if read != assembly:
    print "Fatal error: " + read + " != " + assembly
    sys.exit(1)

k = int(sys.argv[2])
print "Reading unitigs..."
with  open(sys.argv[1]+"/graph_merged.colored.nodes") as unitigFile:
  unitigs=[line.rstrip().split()[2] for line in unitigFile]

for mappingFilename in glob.glob(sys.argv[1]+"/*mapping*"):
  print "Processing " + mappingFilename
  with open(mappingFilename) as mappingFile:
    mapping = mappingFile.readlines()
  for i in xrange(0, len(mapping), 3):
    read = mapping[i].split()[2]
    unitigInfos = parseLine(mapping[i+1])
    checkIfMappingIsOk(read, unitigInfos, unitigs, k)

    unitigInfos = parseLine(mapping[i+2])
    checkIfMappingIsOk(rc(read), unitigInfos, unitigs, k)
  print "All ok!"
