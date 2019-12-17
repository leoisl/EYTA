#Separate the Ensembl transcripts into genes according to the gene field in the header
#E.g.:
#>ENST00000434970.2 cdna chromosome:GRCh38:14:22439007:22439015:1 gene:ENSG00000237235.2 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene gene_symbol:TRDD2 description:T-cell receptor delta diversity 2 [Source:HGNC Symbol;Acc:HGNC:12255]
#CCTTCCTAC
#Will be save to a file called ENSG00000237235.2.fa

#Input:
inputFile = "rawInput/Homo_sapiens.GRCh38.cdna.all.fa"
outputFolder = "input" #input for the simulate.py script

import os
from Bio import SeqIO
from Bio.Seq import Seq
from pprint import pprint

#create output folder
if not os.path.exists(outputFolder):
  os.makedirs(outputFolder)

#process input file
gene2Record={}
for record in SeqIO.parse(inputFile, "fasta"):
  print "Processing " + record.id
  gene = record.description.split()[3].split(":")[1]
  if gene not in gene2Record:
    gene2Record[gene]=[]
  gene2Record[gene].append(record)

#output
for gene in gene2Record:
  with open (outputFolder+"/"+gene+".fa", "w") as outputGene:
    print "Writing " + gene
    SeqIO.write(gene2Record[gene], outputGene, "fasta")


