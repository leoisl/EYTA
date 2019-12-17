#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script checks if two fasta files are equal (disregarding the headers, order and checking reverse-complement also)
"""
import sys

if len(sys.argv)!=3:
	print "Usage:", sys.argv[0], "<EYTA_file_1.fa> <EYTA_file_2.fa>"
	sys.exit(1)


def getSequencesFromFasta(file):
	seqs=[]
	with open(file) as fin:
		for line in fin:
			line=line.rstrip()
			if line[0]!=">" and len(line)>0:
				seqs.append(line.upper())
	seqs.sort()
	return seqs


file1=sys.argv[1]
file2=sys.argv[2]

seqsFile1=getSequencesFromFasta(file1)
seqsFile2=getSequencesFromFasta(file2)

equal=True
for seq in seqsFile1:
	if seq not in seqsFile2:
		print "Sequence %s from file %s is not in file %s"%(seq, file1, file2)
		equal=False

for seq in seqsFile2:
	if seq not in seqsFile1:
		print "Sequence %s from file %s is not in file %s"%(seq, file2, file1)
		equal=False

if equal:
	print "THE FILES ARE EQUAL"
else:
	print "THE FILES ARE DIFFERENT"