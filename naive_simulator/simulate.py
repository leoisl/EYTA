#This script is a very naive simulator
#It simulates perfect transcripts and perfect short reads
#Perfect transcripts are simply the full transcripts with no errors
#Perfect short reads are short reads from the transcripts with no errors
#For each gene, we simulate x% (at least 1) as perfect transcripts and 100% as short reads

#Input:
#1. inputFolder: A folder with one fasta file per gene
#Current input folder is downloaded from ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
inputFolder = "input"
#2. percentagePerfectTranscripts: percentage of the transcripts of a gene that will be perfect transcripts (the nb of transcripts from each gene to simulate is at minimum 1)
percentagePerfectTranscripts = 0.1
#3. shortReadLength: length of the short read
shortReadLength = 101
#4. shortReadSpacing: the i-th short read begins at position (i*shortReadSpacing)
shortReadSpacing = 30
#5. outputFolder: Where the results are going to be stored
outputFolder = "output"

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from pprint import pprint
import numpy as np

def printToOutput(outputShortReads, header, seq, reverse):
  if reverse:
    seq = str(Seq(seq).reverse_complement())
  outputShortReads.write(header + "\n")
  outputShortReads.write(seq + "\n")


for filename in os.listdir(inputFolder):
  #some logging
  fullFilename = inputFolder+"/"+filename

  #save only the sufficiently long records to the allRecords var
  allRecords=[]
  for record in SeqIO.parse(fullFilename, "fasta"):
    if len(record.seq)>=shortReadLength:
      allRecords.append(record)


  #skip trivial genes
  if (len(allRecords)<=1):
    print "Skipping " + fullFilename + " because it has <=1 sequences larger than " + str(shortReadLength)
    continue

  #here, things are fine
  print "Simulating " + fullFilename

  #get the perfect transcripts
  nbOfPerfectTranscripts = max(1, int(percentagePerfectTranscripts*len(allRecords)))
  perfectTranscriptsIndexes = np.array(np.random.choice(len(allRecords), nbOfPerfectTranscripts, False))
  perfectTranscripts = []
  for index in perfectTranscriptsIndexes:
    perfectTranscripts.append(allRecords[index])

  #create output folder
  if not os.path.exists(outputFolder):
    os.makedirs(outputFolder)
  
  #print the perfect transcripts
  with open (outputFolder+"/"+filename+".transcripts.fa", "w") as outputPerfectTranscripts:
    SeqIO.write(perfectTranscripts, outputPerfectTranscripts, "fasta")

  #simulate short reads for all transcripts
  with open (outputFolder+"/"+filename+".short_reads.fa", "w") as outputShortReads:
    for record in allRecords:
      #simulation
      initialPos=0
      readNumber=0
      while (initialPos+shortReadLength<len(record.seq)):
        reverse = bool(np.random.randint(2))
        header = ">{}:{}[{}:{}]:reverse={}".format(fullFilename,record.id, initialPos, initialPos+shortReadLength, str(reverse))
        seq = str(record.seq)[initialPos:initialPos+shortReadLength]
        printToOutput(outputShortReads, header, seq, reverse)
        initialPos+=shortReadSpacing
        readNumber+=1