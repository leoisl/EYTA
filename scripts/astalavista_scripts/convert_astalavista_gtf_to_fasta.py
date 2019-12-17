#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Example usage:
python convert_astalavista_gtf_to_fasta.py --genome /data2/leandro/EYTA/ASTALAVISTA/ASTALAVISTA_run_on_GRCh38.p10_Ensembl90/GRCh38/ --gtf Homo_sapiens.GRCh38.90.gtf --astalavista Homo_sapiens.GRCh38.90.gtf_astalavista.gtf --output Homo_sapiens.GRCh38.90.gtf_astalavista.gtf.all_pairwise_events.fa

This script converts astalavista GTF to FASTA
Let's say we ran ASTALAVISTA in a transcriptome as:
astalavista -t asta -c <genome> e [ASI] -d 2 -i <transcriptome> -a [SEQ]

The lines of the GTF file will look like this:
1       Undefined       as_event        12613   13221   .       +       .       transcript_id "ENST00000450305,ENST00000456328"; gene_id "1:11869-14409W"; flanks "12613-,13221-"; structure "1^3-4^,2^"; splice_chain "12697^12975-13052^,12721^"; sources "Undefined,Undefined"; NMorMRNA "null"; degree "4"; dimension "2_2";

We have to transform this into sequences, so we can compare with EYTA output.

Some of the lines are not transformable. The difference between EYTA and ASTALAVISTA is that EYTA uses sequences and not the annotations, so it needs flanking exons to have at least k bases in common.
ASTALAVISTA uses annotations, so it does not need this, but it needs the same splicing site as flanking.
Sometimes, flanking exons do not have k bases in common, and in this case we will discard this event (hopefully, we won't have many)

"""

import sys
import os
import re
import pprint
import copy
import argparse

class ExonError(Exception):
	pass
class FlankingExonsDoNotCorrespond(Exception):
	pass
class FlankingExonsTooShort(Exception):
	pass
class FlankingExonDoNotExist(Exception):
	pass


class ErrorsLogger:
	"""
	Just holds some global variables related to errors
	"""
	exonErrors = 0
	flankingExonsDoNotCorrespond = 0
	flankingExonsTooShort = 0
	flankingExonDoNotExist = 0
	successful = 0
	

	@staticmethod
	def report():
		print "# of errors logged: "
		print "exonErrors = ", ErrorsLogger.exonErrors
		print "flankingExonsDoNotCorrespond = ", ErrorsLogger.flankingExonsDoNotCorrespond
		print "flankingExonsTooShort = ", ErrorsLogger.flankingExonsTooShort
		print "FlankingExonDoNotExist = ", ErrorsLogger.flankingExonDoNotExist
		print "successful = ", ErrorsLogger.successful

def loadGenome(filenames):
	"""
	load the sequences in filenames
	@param: filenames - a list of .fa files
	"""
	print "Loading genome: ", filenames
	chr2Seq={}
	for filename in filenames:
		with open(filename) as file:
			name = os.path.basename(filename)[:-3]

			#read the sequence
			seqLines=[]
			for line in file:
				if line[0]==">":
					continue
				seqLines.append(line.rstrip())

			chr2Seq[name]="".join(seqLines)

	print "Done!"
	return chr2Seq

class GTFExons:
	"""
	Represents the exons of each transcript in the GTF
	"""
	def __init__(self):
		self.transcriptInfo={}

	def addExon(self, lineSplit):
		"""
		Add an exon to this index
		line has the following format:
		1       havana  exon    11869   12227   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2";
		exon_number "1"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana";
		transcript_biotype "processed_transcript"; exon_id "ENSE00002234944"; exon_version "1"; tag "basic"; transcript_support_level "1";
		"""
		#get what interest us
		chromossome = lineSplit[0]
		strand = lineSplit[6]
		acceptor = int(lineSplit[3])
		donor = int(lineSplit[4])

		if acceptor > donor: #this is weird, we should inverse them
			acceptor, donor = donor, acceptor

		transcript = lineSplit[lineSplit.index("transcript_id")+1][1:-2] #[1:-2] is to remove the begin and the end stuff that we don't care, i.e. transcript_id "ENST00000456328";
		geneName = lineSplit[lineSplit.index("gene_name")+1][1:-2]
		transcriptName = lineSplit[lineSplit.index("transcript_name") + 1][1:-2]

		if transcript not in self.transcriptInfo:
			self.transcriptInfo[transcript]={"geneName": geneName, "transcript": transcript, "transcriptName": transcriptName, "chr": chromossome, "exons":[]}

		self.transcriptInfo[transcript]["exons"].append({"start": acceptor, "end": donor}) #I put start since it can be an acceptor or transcriptional start or w/e... Same for end
		self.transcriptInfo[transcript]["exons"].sort(cmp=lambda x, y : x["start"]-y["start"])


	def __getExonNbContainingThisSS(self, transcript, spliceSiteCoordinate):
		"""
		Given a spliceSiteCoordinate (e.g.: ('12613', 'start') or ('13052', 'end') ), returns the exon of it
		"""
		#check the transcript
		if transcript not in self.transcriptInfo:
			raise Exception('Transcript ' + transcript + " not found")
		exonList = self.transcriptInfo[transcript]["exons"]

		#look for the given splice site
		for exonNumber, exon in enumerate(exonList):
			if exon[spliceSiteCoordinate[1]] == spliceSiteCoordinate[0]:
				return exonNumber

		raise Exception("Splice site %s not found in transcript %s"%(spliceSiteCoordinate, transcript))

	def getPreviousExonBeforeThisSS(self, transcript, spliceSiteCoordinate):
		"""
		Given a spliceSiteCoordinate (e.g.: ('12613', 'start') or ('13052', 'end') ), search the previous exon before this splice site.
		If the splice site is a end site, then the previous exon is the exon of this end site
		Otherwise, if it is an start site, then the previous exon is the exon that comes before the exon of this splice site
		If the splicesite is not found in the list of exons, we raise an exception
		"""
		exonNumber = self.__getExonNbContainingThisSS(transcript, spliceSiteCoordinate)

		if spliceSiteCoordinate[1]=="start":
			if exonNumber-1<0:
				raise FlankingExonDoNotExist("No exon before 0")
				#I don't think we should do this, raising the exception is better
				#previousExon = self.transcriptInfo[transcript]["exons"][exonNumber]
			else:
				previousExon = self.transcriptInfo[transcript]["exons"][exonNumber-1]
		else:
			previousExon = self.transcriptInfo[transcript]["exons"][exonNumber]
		return previousExon

	def getNextExonAfterThisSS(self, transcript, spliceSiteCoordinate):
		"""
		Given a spliceSiteCoordinate (e.g.: ('12613', 'start') or ('13052', 'end') ), search the next exon after this splice site.
		If the splice site is a end site, then the next exon is the next exon after the exon of this end site
		Otherwise, if it is an start site, then the next exon is the exon of this start
		If the splicesite is not found in the list of exons, we raise an exception
		"""
		exonNumber = self.__getExonNbContainingThisSS(transcript, spliceSiteCoordinate)

		if spliceSiteCoordinate[1] == "start":
			nextExon = self.transcriptInfo[transcript]["exons"][exonNumber]
		else:
			if exonNumber+1 >= len(self.transcriptInfo[transcript]["exons"]):
				raise FlankingExonDoNotExist("No exon after end")
				#I don't think we should do this, raising the exception is better
				#nextExon = self.transcriptInfo[transcript]["exons"][exonNumber]
			else:
				nextExon = self.transcriptInfo[transcript]["exons"][exonNumber+1]
		return nextExon

	def getExonOfTheseBothSS(self, transcript, spliceSiteCoordinateStart, spliceSiteCoordinateEnd):
		"""
		Given two spliceSiteCoordinates (e.g.: ('12613', 'start') or ('13052', 'end') ), search for the exon with these both splice sites
		"""
		exonStartNb = self.__getExonNbContainingThisSS(transcript, spliceSiteCoordinateStart)
		exonEndNb = self.__getExonNbContainingThisSS(transcript, spliceSiteCoordinateEnd)
		if exonStartNb == exonEndNb:
			return self.transcriptInfo[transcript]["exons"][exonStartNb]
		else:
			raise ExonError("Start %s and end %s do not have the same exons in transcript %s"%(spliceSiteCoordinateStart, spliceSiteCoordinateEnd, transcript))


def loadTranscriptome(transcritomeFilename):
	"""
	Load the transcriptome into a GTFExons object
	"""
	print "Loading transcriptome: ", transcritomeFilename
	gtf = GTFExons()
	with open(transcritomeFilename) as transcritomeFile:
		for line in  transcritomeFile:
			#skip comments
			if line[0]=="#":
				continue

			#process only exons
			lineSplit = line.rstrip().split()
			if lineSplit[2]!="exon":
				continue

			#add the exon
			gtf.addExon(lineSplit)

	print "Done!"
	return gtf


ssCoordinateType = re.compile("(\d+)([\^\-\(\)\[\]]?)")

def toStartOrEnd(ASTALAVISTASymbol):
	if ASTALAVISTASymbol in ["-", "(", "["]:
		return "start"
	elif ASTALAVISTASymbol in ["^", ")", "]"]:
		return "end"
	else:
		raise Exception("Unrecognized ASTALAVISTA symbol %s"%ASTALAVISTASymbol)

def processASTALAVISTAField(field):
	"""
	Process a field like '"5584810-5585063^,5579670-5579808^5584810-5584975^";' and returns:
	[[[5584810, 'start'], [5585063, 'end']],
     [[5579670, 'start'], [5579808, 'end'], [5584810, 'start'], [5584975, 'end']]]
	"""
	fieldSplit = field.replace('"', '').replace(";", "").split(",")
	ASCodesFormatted=[]
	for ASCode in fieldSplit:
		#get the regex matches
		regexMatches = re.findall(ssCoordinateType, ASCode)
		formatted = []
		for ASCodeAlmostFormatted in regexMatches:
			formatted.append([int(ASCodeAlmostFormatted[0]), toStartOrEnd(ASCodeAlmostFormatted[1])])
		ASCodesFormatted.append(formatted)
	return ASCodesFormatted
#test
#pprint.pprint(processASTALAVISTAField('"5584810-5585063^,5579670-5579808^5584810-5584975^";'))
#sys.exit(1)

def getFlankingExonList(transcript, ASCode, gtf):
	"""
	Get the flanking exon list as an ASCode given as [('12613', 'start'), ('12697', 'end'), ('12975', 'start'), ('13052', 'end'), ('13221', 'start')]
	Flanking means that if the first SS is a start, then we'll get the exon before it
	If the last SS is an end, then we'll get the exon after it
	"""
	exonList=[]

	#we put the first exon in the list - the one coming before the first SS
	exonList.append(gtf.getPreviousExonBeforeThisSS(transcript, ASCode[0]))

	#startIndexInternal is the start of the exon after exonList[0] (the left flanking exon)
	#which can be ASCode[0] or ASCode[1]
	#if ASCode[0] is a start, then startIndex is 0, else is 1
	startIndexInternal = 0 if ASCode[0][1]=="start" else 1

	#similar for endIndex
	endIndexInternal = len(ASCode) if ASCode[-1][1]=="end" else len(ASCode)-1
	for i in range(startIndexInternal, endIndexInternal, 2):
		exonList.append(gtf.getExonOfTheseBothSS(transcript, ASCode[i], ASCode[i+1]))

	#we put the last exon in the list - the one after the last SS
	exonList.append(gtf.getNextExonAfterThisSS(transcript, ASCode[-1]))

	#we can have duplicates here - we remove them
	exonListUnique = []
	for exon in exonList:
		if exon not in exonListUnique:
			exonListUnique.append(exon)
		else:
			print "[ERROR]: duplicated exons in : %s"%(str(ASCode))
			print "I do not really know why and when this happens, so I am printing it and if you find this error, please report"
			sys.exit(1)

	return exonListUnique


def getSeqFromExonList(exonList, transcript, gtf, genome):
	# we get the sequences themselves
	chromossome = gtf.transcriptInfo[transcript]["chr"]
	seqList = [genome[chromossome][exon['start'] - 1:exon['end']] \
			   for exon in exonList]  # exon['start']-1 to deal with the 1-based gtf coordinates (should do for the 'end' also, but in the end we don't need since gtf coordinates are inclusive)
	return "".join(seqList)


def checkIfExonHasAtLeastKBases(exon, k):
	if exon['end']-exon['start']+1 < k:
		raise FlankingExonsTooShort("Flanking exon with < k = %d bases - %s" % (k, exon))

def reverseSpliceChain(spliceChain):
	spliceChain=spliceChain[::-1]
	for SS in spliceChain:
		if SS[1]=="start":
			SS[1]="end"
		else:
			SS[1]="start"
	return spliceChain

bubbleId=0
def fromASTALAVISTALineToFastaBubble(line, genome, gtf, k, outputFile):
	"""
	MAIN FUNCTION - transforms ASTALAVISTA GTF line to sequences
	@param line: astalavista GTF line, e.g.:
	1       Undefined       as_event        12613   13221   .       +       .       transcript_id "ENST00000450305,ENST00000456328"; gene_id "1:11869-14409W"; flanks "12613-,13221-"; structure "1^3-4^,2^"; splice_chain "12697^12975-13052^,12721^"; sources "Undefined,Undefined"; NMorMRNA "null"; degree "4"; dimension "2_2";
	"""
	lineSplit = line.split()

	#get the transcripts
	transcripts1, transcripts2 = lineSplit[lineSplit.index("transcript_id")+1].replace('"', '').replace(';', '').split(',')
	transcripts1 = transcripts1.split("/")
	transcripts2 = transcripts2.split("/")

	#get the left and the right flanks
	flanks = lineSplit[lineSplit.index("flanks")+1]
	leftFlank, rightFlank = processASTALAVISTAField(flanks)

	#get the splice_chain
	spliceChain = lineSplit[lineSplit.index("splice_chain")+1]
	spliceChain1, spliceChain2 = processASTALAVISTAField(spliceChain)

	#check if it is in the minus strand - we should inverse
	if lineSplit[6]=="-":
		"""
		DEBUG code, but seems right!
		print "it is on the minus strand:"
		print "Before reversing:"
		print "line: %s\nleftFlank: %s\nrightFlank: %s\nspliceChain1: %s\nspliceChain2: %s"%(line, str(leftFlank), str(rightFlank), str(spliceChain1), str(spliceChain2))
		"""
		leftFlank, rightFlank = reverseSpliceChain(rightFlank), reverseSpliceChain(leftFlank)
		spliceChain1=reverseSpliceChain(spliceChain1)
		spliceChain2 = reverseSpliceChain(spliceChain2)
		"""
		DEBUG code, but seems right!
		print "after reversing:"
		print "line: %s\nleftFlank: %s\nrightFlank: %s\nspliceChain1: %s\nspliceChain2: %s"%(line, str(leftFlank), str(rightFlank), str(spliceChain1), str(spliceChain2))
		"""



	#get the two local variations
	ASCode1=leftFlank+spliceChain1+rightFlank
	ASCode2=leftFlank+spliceChain2+rightFlank

	#transform each local variation into a sequence
	global bubbleId
	for transcript1 in transcripts1:
		for transcript2 in transcripts2:
			try:
				exonList1 = copy.deepcopy(getFlankingExonList(transcript1, ASCode1, gtf)) #we make a deep copy because we will change the exons here
				exonList2 = copy.deepcopy(getFlankingExonList(transcript2, ASCode2, gtf)) #we make a deep copy because we will change the exons here

				#we transform the flanking exons to have k chars, at most
				checkIfExonHasAtLeastKBases(exonList1[0], k)
				checkIfExonHasAtLeastKBases(exonList1[-1], k)
				checkIfExonHasAtLeastKBases(exonList2[0], k)
				checkIfExonHasAtLeastKBases(exonList2[-1], k)

				exonList1[0]["start"] = exonList1[0]["end"] - k + 1
				exonList1[-1]["end"] = exonList1[-1]["start"] + k - 1
				exonList2[0]["start"] = exonList2[0]["end"] - k + 1
				exonList2[-1]["end"] = exonList2[-1]["start"] + k - 1

				#do the checks - the flanking exons must be exactly the same
				if exonList1[0]!=exonList2[0]:
					raise FlankingExonsDoNotCorrespond("LEFT Flanking exons errors:\n%s\n%s" % (exonList1, exonList2))
				if exonList1[-1]!=exonList2[-1]:
					raise FlankingExonsDoNotCorrespond("RIGHT Flanking exons errors:\n%s\n%s" % (exonList1, exonList2))



				# we build the headers
				transcriptInfo1 = gtf.transcriptInfo[transcript1]
				header1 = ">Bubble_%d_Gene_%s_Transcript_%s_TranscriptId_%s_Splicing_1_Flanks_%s_SpliceChain_%s" % (bubbleId, transcriptInfo1["geneName"], transcriptInfo1["transcriptName"], transcriptInfo1["transcript"], flanks, spliceChain)
				transcriptInfo2 = gtf.transcriptInfo[transcript2]
				header2 = ">Bubble_%d_Gene_%s_Transcript_%s_TranscriptId_%s_Splicing_2_Flanks_%s_SpliceChain_%s" % (bubbleId, transcriptInfo2["geneName"], transcriptInfo2["transcriptName"], transcriptInfo2["transcript"], flanks, spliceChain)
				bubbleId+=1

				#get the sequences
				seq1 = getSeqFromExonList(exonList1, transcript1, gtf, genome)
				seq2 = getSeqFromExonList(exonList2, transcript2, gtf, genome)

				outputFile.write("%s\n%s\n%s\n%s\n"%(header1, seq1, header2, seq2))


				ErrorsLogger.successful+=1
			except ExonError as e:
				#print e
				ErrorsLogger.exonErrors += 1
			except FlankingExonsDoNotCorrespond as e:
				#print e
				ErrorsLogger.flankingExonsDoNotCorrespond += 1
			except FlankingExonsTooShort as e:
				#print e
				ErrorsLogger.flankingExonsTooShort += 1
			except FlankingExonDoNotExist as e:
				#print e
				ErrorsLogger.flankingExonDoNotExist += 1


#test
#fromASTALAVISTALineToFastaBubble('1       Undefined       as_event        12613   13221   .       +       .       transcript_id "ENST00000450305,ENST00000456328"; gene_id "1:11869-14409W"; flanks "12613-,13221-"; structure "1^3-4^,2^"; splice_chain "12697^12975-13052^,12721^"; sources "Undefined,Undefined"; NMorMRNA "null"; degree "4"; dimension "2_2";', None, None)




if __name__=="__main__":
	#parse the args
	parser = argparse.ArgumentParser(description='Convert ASTALAVISTA GTF to fasta sequences')
	parser.add_argument('--genome', dest="genomeFolder", required=True, help='a path to a folder in which each file is a chromosome')
	parser.add_argument('--gtf', dest="transcritomeFile", required=True, help='a path to the gtf file that was used as input to ASTALAVISTA')
	parser.add_argument('--astalavista', dest="astalavistaGTFFilename", required=True, help='a path to ASTALAVISTA gtf file')
	parser.add_argument('--output', dest="outputFilename", help='output file', default="output.fa")
	parser.add_argument('-k', dest="k", type=int, help='kmer value - should be the one you used to run EYTA', default=31)
	args = parser.parse_args()


	#load the genome
	onlyFAfiles = [os.path.join(args.genomeFolder, f) for f in os.listdir(args.genomeFolder) if os.path.isfile(os.path.join(args.genomeFolder, f)) and os.path.basename(f)[-3:]==".fa"]
	genome = loadGenome(onlyFAfiles)

	#load the transcriptome
	gtf = loadTranscriptome(args.transcritomeFile)

	#load the astalavista GTF file
	with open(args.astalavistaGTFFilename) as astalavistaGTFFile, open(args.outputFilename, "w") as outputFile:
		for line in astalavistaGTFFile:
			#skip comments
			if line[0]=="#":
				continue

			#transform each line into sequences
			fromASTALAVISTALineToFastaBubble(line.rstrip(), genome, gtf, args.k, outputFile)

	#report the errors and success
	ErrorsLogger.report()