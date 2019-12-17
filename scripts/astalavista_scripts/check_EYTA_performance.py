"""
Executing it:
see README
"""

"""
Reading Homo_sapiens.GRCh38.90.gtf.cov_30x.short_reads.fastq we know the transcripts that were expressed in the short reads
Reading Homo_sapiens.GRCh38.90.gtf.shallow_transcriptome_0.1.fa  we know the transcripts that are in the shallow transcriptome
Then, we read Homo_sapiens.GRCh38.90.gtf_astalavista.gtf.all_pairwise_events.fa to get all bubbles found by ASTALAVISTA, but keeping only the ones
	such that one transcript is expressed in SR+LR and the other only in the SR.
"""
import sys
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import numpy as np
import re

getTranscriptPattern = re.compile("Transcript=([a-zA-Z0-9]+);")

class HeaderSeq:
	"""
	This class stores the header and the seq of a fasta line.
	We need to define this as class and implement __hash__ and __eq__ in order for objects to be hashable and be able to be dict's keys
	"""
	def __init__(self, header, seq):
		self.__header = header
		self.__seq = seq

	def __eq__(self, other):
		return (self.__header, self.__seq) == (other.__header, other.__seq)

	def __ne__(self, other):
		return not self.__eq__(other)

	def __str__(self):
		return ">%s\n%s"%(self.__header, self.__seq)

	def reverseComplement(self):
		"""
		Returns a HeaderSeq which is the RC of this one (the header keeps the same)
		"""
		return HeaderSeq(self.__header, self.__seq.reverse_complement())

	def getHeader(self):
		return self.__header

	def getSeq(self):
		return self.__seq

	def isNovel(self):
		return "Transcript=Novel;" in self.getHeader()

	def getTranscript(self):
		return getTranscriptPattern.search(self.getHeader()).group(1)

class Bubble:
	"""
	Represents a bubble: a pair of HeaderSeq
	"""
	def __addKmers(self, seq):
		k=31
		for i in range(len(seq)-k+1):
			self.kmers.add(seq[i:i+k])

	def __init__(self, upperPath, lowerPath, isRC=False):
		self.__upperPath = upperPath
		self.__lowerPath = lowerPath
		
		#add all k-mers
		self.kmers=set()
		self.__addKmers(self.__upperPath.getSeq())
		self.__addKmers(self.__lowerPath.getSeq())

		#precompute the RC
		if isRC==False:
			self.reverseComplementBubble = Bubble(self.__upperPath.reverseComplement(), self.__lowerPath.reverseComplement(), True)


	def __eq__(self, other):
		"""
		A bubble is equal to the other just sequence-wise, and we don't also care about the strand
		"""
		reverseComplementedBubble = self.reverseComplement()
		return (self.__upperPath.getSeq() == other.__upperPath.getSeq() and self.__lowerPath.getSeq() == other.__lowerPath.getSeq()) or \
			   (self.__lowerPath.getSeq() == other.__upperPath.getSeq() and self.__upperPath.getSeq() == other.__lowerPath.getSeq()) or \
			   (reverseComplementedBubble.__upperPath.getSeq() == other.__upperPath.getSeq() and reverseComplementedBubble.__lowerPath.getSeq() == other.__lowerPath.getSeq()) or \
			   (reverseComplementedBubble.__lowerPath.getSeq() == other.__upperPath.getSeq() and reverseComplementedBubble.__upperPath.getSeq() == other.__lowerPath.getSeq())

	def __ne__(self, other):
		return not self.__eq__(other)

	def __str__(self):
		return "%s\n%s"%(str(self.__upperPath), str(self.__lowerPath))

	def reverseComplement(self):
		#return Bubble(self.__upperPath.reverseComplement(), self.__lowerPath.reverseComplement())
		return self.reverseComplementBubble

	def getUpperPath(self):
		return self.__upperPath

	def getLowerPath(self):
		return self.__lowerPath

	def shouldBeFoundByEYTA(self, transcriptsExpressedInShortReads, shallowTranscriptome):
		"""
		For this to be true, one path should be in the transcriptome + SR, and the other only in the SR
		"""
		def isExpressedInShortReads(path):
			"""
			Returns true if the transcript is expressed in the short reads
			It might not mean that the whole transcript is covered, but in our case, yes
			"""
			def getTranscriptIdFromASTALAVISTA(header):
				headerSplit = header.split("_")
				return headerSplit[headerSplit.index("TranscriptId")+1]
			return getTranscriptIdFromASTALAVISTA(path.getHeader()) in transcriptsExpressedInShortReads

		def isExpressedInShallowTranscriptome(path):
			for transcript in shallowTranscriptome:
				if path.getSeq() in transcript.getSeq() or path.getSeq().reverse_complement() in transcript.getSeq():
					return True
			return False


		def checkPaths(path1, path2):
			return	isExpressedInShallowTranscriptome(path1) and \
					isExpressedInShortReads(path1) and \
					not isExpressedInShallowTranscriptome(path2) and \
					isExpressedInShortReads(path2)

		return checkPaths(self.__upperPath, self.__lowerPath) or checkPaths(self.__lowerPath, self.__upperPath)

	@staticmethod
	def __computeEditDistanceOfAShortStringInsideALongerStringCore(longString, shortString):
		len1 = len(longString); len2=len(shortString)
		col=[0]*(len2+1)
		prevCol=[i for i in range(len2+1)]
		editDistance = prevCol[-1]
		bestIndex=0

		for i in range(len1):
			col[0]=0
			for j in range(len2):
				col[j+1] = min(prevCol[j+1] + 1, col[j] + 1, prevCol[j] + (0 if longString[i]==shortString[j] else 1))

			if col[-1] < editDistance:
				editDistance=col[-1]
				bestIndex=i

			col,prevCol=prevCol,col

		return editDistance

	@staticmethod
	def __getLocalAlignment(str1, str2):
		return pairwise2.align.localms(str1, str2, 2, -1, -10000, -10000, one_alignment_only=True)[0]

	def isInsideExact(self, other):
		#first try an exact match with and 1-off base
		reverseComplementedBubble = self.reverseComplement()
		return (self.__upperPath.getSeq() in other.__upperPath.getSeq() and self.__lowerPath.getSeq() in other.__lowerPath.getSeq()) or \
			   (self.__lowerPath.getSeq() in other.__upperPath.getSeq() and self.__upperPath.getSeq() in other.__lowerPath.getSeq()) or \
			   (reverseComplementedBubble.__upperPath.getSeq() in other.__upperPath.getSeq() and reverseComplementedBubble.__lowerPath.getSeq() in other.__lowerPath.getSeq()) or \
			   (reverseComplementedBubble.__lowerPath.getSeq() in other.__upperPath.getSeq() and reverseComplementedBubble.__upperPath.getSeq() in other.__lowerPath.getSeq()) or \
			   (self.__upperPath.getSeq()[1:-1] in other.__upperPath.getSeq() and self.__lowerPath.getSeq() in other.__lowerPath.getSeq()) or \
			   (self.__lowerPath.getSeq()[1:-1] in other.__upperPath.getSeq() and self.__upperPath.getSeq() in other.__lowerPath.getSeq()) or \
			   (reverseComplementedBubble.__upperPath.getSeq()[1:-1] in other.__upperPath.getSeq() and reverseComplementedBubble.__lowerPath.getSeq() in other.__lowerPath.getSeq()) or \
			   (reverseComplementedBubble.__lowerPath.getSeq()[1:-1] in other.__upperPath.getSeq() and reverseComplementedBubble.__upperPath.getSeq() in other.__lowerPath.getSeq())

	def __shareEnoughKmers(self, other):
		threshold=10 #what would be enough?
		return len(self.kmers & other.kmers)>=threshold or len(self.kmers & other.reverseComplement().kmers)>=threshold

	def getInexactScore(self, other):
		'''
		@returns: the score of the best local alignment and the local alignment
		'''
		if self.__shareEnoughKmers(other):
			reverseComplementedBubble = self.reverseComplement()

			#computes all the local alignments
			LA={i:{} for i in range(4)}
			LA[0][0] = Bubble.__getLocalAlignment(self.__upperPath.getSeq(), other.__upperPath.getSeq())
			LA[0][1] = Bubble.__getLocalAlignment(self.__lowerPath.getSeq(), other.__lowerPath.getSeq())

			LA[1][0] = Bubble.__getLocalAlignment(self.__lowerPath.getSeq(), other.__upperPath.getSeq())
			LA[1][1] = Bubble.__getLocalAlignment(self.__upperPath.getSeq(), other.__lowerPath.getSeq())

			LA[2][0] = Bubble.__getLocalAlignment(reverseComplementedBubble.__upperPath.getSeq(), other.__upperPath.getSeq())
			LA[2][1] = Bubble.__getLocalAlignment(reverseComplementedBubble.__lowerPath.getSeq(), other.__lowerPath.getSeq())

			LA[3][0] = Bubble.__getLocalAlignment(reverseComplementedBubble.__lowerPath.getSeq(), other.__upperPath.getSeq())
			LA[3][1] = Bubble.__getLocalAlignment(reverseComplementedBubble.__upperPath.getSeq(), other.__lowerPath.getSeq())

			#get the arg of the best score
			bestScoreIndex = np.argmax(np.array([LA[i][0][2] + LA[i][1][2] for i in range(4)]))

			#returns the score and the best alignment
			return LA[bestScoreIndex][0][2]+LA[bestScoreIndex][1][2], LA[bestScoreIndex]

		return -1.0, []







def getTranscriptsExpressedInShortReads(srFilename):
	print  >> sys.stderr, "Getting Transcripts Expressed In Short Reads..."
	transcriptsExpressedInShortReads = set()
	for record in SeqIO.parse(srFilename, "fastq"):
		transcript = record.id.split("_")[0]
		transcriptsExpressedInShortReads.add(transcript)
	print >> sys.stderr,  "Done!"
	return transcriptsExpressedInShortReads



def getTranscriptome(transcriptomeFilename):
	print  >> sys.stderr, "Parsing Transcriptome %s..."%transcriptomeFilename
	transcriptome = [HeaderSeq(record.id, record.seq) for record in SeqIO.parse(transcriptomeFilename, "fasta")]
	print  >> sys.stderr, "Done!"
	return transcriptome
#test:
#print getTranscriptIdFromASTALAVISTA('Gene_FP565260.3_Transcript_FP565260.3-201_TranscriptId_ENST00000612610_Splicing_1_Flanks_"5025049^,5026630^";_SpliceChain_"5026280-,5026432-";')
#sys.exit(0)


def getASTALAVISTABubblesToBeFound(transcriptsExpressedInShortReads, shallowTranscriptome, ASTALAVISTA_AS_bubbles_filename):
	"""
	Returns a set where each element is a Bubble containing the ASTALAVISTA bubbles to be found by EYTA
	"""
	print >> sys.stderr,  "Getting ASTALAVISTA bubbles to be found..."
	ASTALAVISTABubblesToBeFound=[]

	allASTALAVISTASeqs=[HeaderSeq(record.id, record.seq) for record in SeqIO.parse(ASTALAVISTA_AS_bubbles_filename, "fasta")]

	totalBubbles = 0
	repeatedBubbles = 0
	for i in xrange(0, len(allASTALAVISTASeqs), 2):
		print  >> sys.stderr, "Checking ASTALAVISTA bubble %d..."%i
		bubble = Bubble(allASTALAVISTASeqs[i], allASTALAVISTASeqs[i+1])
		if bubble.shouldBeFoundByEYTA(transcriptsExpressedInShortReads, shallowTranscriptome):
			totalBubbles+=1
			if bubble not in ASTALAVISTABubblesToBeFound:
				#print  >> sys.stderr, "[getASTALAVISTABubblesToBeFound]: ADDING:\n%s"%str(bubble)
				ASTALAVISTABubblesToBeFound.append(bubble)
			else:
				repeatedBubbles+=1
				#print  >> sys.stderr, "[getASTALAVISTABubblesToBeFound]: NO NEED TO ADD:\n%s"%str(bubble)
	print "We had %d ASTALAVISTA bubbles to be found by EYTA, but %d are repeated (represent the same event in different transcripts or strands)"%(totalBubbles, repeatedBubbles)
	return ASTALAVISTABubblesToBeFound


def getASTALAVISTABubblesFoundByEYTA(ASTALAVISTABubblesToBeFound, EYTA_AS_bubbles_filename):
	EYTA_AS_seqs = [ HeaderSeq(record.id, record.seq) for record in SeqIO.parse(EYTA_AS_bubbles_filename, "fasta") ]
	EYTA_AS_bubbles = []
	totalBubbles = 0
	repeatedBubbles = 0
	for i in xrange(0, len(EYTA_AS_seqs), 2):
		print  >> sys.stderr, "Checking EYTA bubble %d..."%i
		bubble = Bubble(EYTA_AS_seqs[i], EYTA_AS_seqs[i+1])
		totalBubbles+=1
		if bubble not in EYTA_AS_bubbles:
			#print  >> sys.stderr, "[getASTALAVISTABubblesFoundByEYTA]: ADDING:\n%s"%str(bubble)
			EYTA_AS_bubbles.append(bubble)
		else:
			repeatedBubbles+=1
			#print  >> sys.stderr, "[getASTALAVISTABubblesFoundByEYTA]: NO NEED TO ADD:\n%s"%str(bubble)
	print "We had %d EYTA bubbles, but %d are repeated (represent the same event in different transcripts or strands)"%(totalBubbles, repeatedBubbles)
	print "%d EYTA bubbles considered"%(totalBubbles - repeatedBubbles)

	ASTALAVISTABubble2EYTABubblesMappingIntoIt={index: [] for index, ASTALAVISTABubble in enumerate(ASTALAVISTABubblesToBeFound)}
	EYTABubble2ASTALAVISTABubbleItMapsTo = {index: [] for index, EYTA_AS_bubble in enumerate(EYTA_AS_bubbles)}

	for ASTALAVISTABubbleIndex, ASTALAVISTABubble in enumerate(ASTALAVISTABubblesToBeFound):
		print >> sys.stderr, "Processing ASTALAVISTA and EYTA bubbles... %.2f"%(float(ASTALAVISTABubbleIndex)/len(ASTALAVISTABubblesToBeFound)*100.0)
		for EYTA_AS_bubbleIndex, EYTA_AS_bubble in enumerate(EYTA_AS_bubbles):
			if ASTALAVISTABubble.isInsideExact(EYTA_AS_bubble) or EYTA_AS_bubble.isInsideExact(ASTALAVISTABubble):
				ASTALAVISTABubble2EYTABubblesMappingIntoIt[ASTALAVISTABubbleIndex].append(EYTA_AS_bubble)
				EYTABubble2ASTALAVISTABubbleItMapsTo[EYTA_AS_bubbleIndex].append(ASTALAVISTABubble)

		if len(ASTALAVISTABubble2EYTABubblesMappingIntoIt[ASTALAVISTABubbleIndex])==0:
			#try inexact mapping - get the best inexact mapping score
			bestInexactMapping={"score": 0.0, "ASTALAVISTABubbleIndex": None, "EYTA_AS_bubbleIndex": None}
			for EYTA_AS_bubbleIndex, EYTA_AS_bubble in enumerate(EYTA_AS_bubbles):
				score, alignment = ASTALAVISTABubble.getInexactScore(EYTA_AS_bubble)
				if score > bestInexactMapping["score"]:
					bestInexactMapping = {"score": score, "ASTALAVISTABubbleIndex": ASTALAVISTABubbleIndex, "EYTA_AS_bubbleIndex": EYTA_AS_bubbleIndex}

			#add it
			if bestInexactMapping["ASTALAVISTABubbleIndex"] is not None:
				ASTALAVISTABubble2EYTABubblesMappingIntoIt[bestInexactMapping["ASTALAVISTABubbleIndex"]].append(EYTA_AS_bubbles[bestInexactMapping["EYTA_AS_bubbleIndex"]])
				EYTABubble2ASTALAVISTABubbleItMapsTo[bestInexactMapping["EYTA_AS_bubbleIndex"]].append(ASTALAVISTABubblesToBeFound[bestInexactMapping["ASTALAVISTABubbleIndex"]])


	#goes through all EYTA bubbles that have no mapping and do inexact mapping on them
	#TODO: unsure to do this one...
	'''
	with open("inexactMappingEYTABubblesDebug", "w") as fileDebugInexactMapping :
		for EYTA_AS_bubbleIndex, EYTA_AS_bubble in enumerate(EYTA_AS_bubbles):
			print >> sys.stderr, "Processing unmmaped EYTA bubbles... %.2f"%(float(EYTA_AS_bubbleIndex)/len(EYTA_AS_bubbles)*100.0)
			if len(EYTABubble2ASTALAVISTABubbleItMapsTo[EYTA_AS_bubbleIndex])==0:
				#try inexact mapping - get the best inexact mapping score
				fileDebugInexactMapping.write("Trying inexact mapping on EYTA bubbles...\n")
				bestInexactMapping={"score": 0.0, "ASTALAVISTABubbleIndex": None, "EYTA_AS_bubbleIndex": None}
				bestAlignment=None
				for ASTALAVISTABubbleIndex, ASTALAVISTABubble in enumerate(ASTALAVISTABubblesToBeFound):
					score,alignment = ASTALAVISTABubble.getInexactScore(EYTA_AS_bubble)
					if score > bestInexactMapping["score"]:
						bestInexactMapping = {"score": score, "ASTALAVISTABubbleIndex": ASTALAVISTABubbleIndex, "EYTA_AS_bubbleIndex": EYTA_AS_bubbleIndex}
						bestAlignment = alignment

				#add it
				if bestInexactMapping["ASTALAVISTABubbleIndex"] is not None:
					#TODO: not sure about this one...
					#ASTALAVISTABubble2EYTABubblesMappingIntoIt[bestInexactMapping["ASTALAVISTABubbleIndex"]].append(EYTA_AS_bubbles[bestInexactMapping["EYTA_AS_bubbleIndex"]])
					EYTABubble2ASTALAVISTABubbleItMapsTo[bestInexactMapping["EYTA_AS_bubbleIndex"]].append(ASTALAVISTABubblesToBeFound[bestInexactMapping["ASTALAVISTABubbleIndex"]])
					fileDebugInexactMapping.write("Worked:\n")
					fileDebugInexactMapping.write(format_alignment(*bestAlignment[0])+"\n")
					fileDebugInexactMapping.write(format_alignment(*bestAlignment[1])+"\n")
				else:
					fileDebugInexactMapping.write("DID NOT Work\n")
	'''

	ASTALAVISTABubblesFoundByEYTA = [ASTALAVISTABubblesToBeFound[key] for key,value in ASTALAVISTABubble2EYTABubblesMappingIntoIt.items() if len(value)>0]
	ASTALAVISTABubblesNotFoundByEYTA = [ASTALAVISTABubblesToBeFound[key] for key,value in ASTALAVISTABubble2EYTABubblesMappingIntoIt.items() if len(value)==0]
	EYTABubblesThatMapToNoASTALAVISTABubbles = [EYTA_AS_bubbles[key] for key,value in EYTABubble2ASTALAVISTABubbleItMapsTo.items() if len(value)==0]
	EYTABubblesThatMapToExactlyOneASTALAVISTABubble = [EYTA_AS_bubbles[key] for key,value in EYTABubble2ASTALAVISTABubbleItMapsTo.items() if len(value)==1]
	EYTABubblesThatMapToMoreThanOneASTALAVISTABubbles = [EYTA_AS_bubbles[key] for key,value in EYTABubble2ASTALAVISTABubbleItMapsTo.items() if len(value)>1]

	return ASTALAVISTABubblesFoundByEYTA, ASTALAVISTABubblesNotFoundByEYTA, EYTABubblesThatMapToNoASTALAVISTABubbles, \
		   EYTABubblesThatMapToExactlyOneASTALAVISTABubble, EYTABubblesThatMapToMoreThanOneASTALAVISTABubbles, EYTA_AS_bubbles


def parseGTF(trGTFFilename):
	print >> sys.stderr, "Parsing %s..."%trGTFFilename
	gene2Transcripts, transcript2Gene = {}, {}

	with open(trGTFFilename) as trGTFFile:
		for line in trGTFFile:
			lineSplit = line.split()
			if lineSplit[2]=="transcript":
					geneId = lineSplit[9][1:-2]
					transcriptId = lineSplit[13][1:-2]
					
					if geneId not in gene2Transcripts:
						gene2Transcripts[geneId]=[]
					gene2Transcripts[geneId].append(transcriptId)
					transcript2Gene[transcriptId]=geneId
	print >> sys.stderr, "Done!"
	return gene2Transcripts, transcript2Gene


def novelPathIsInTranscriptomeMinusShallowTranscriptome(transcriptIdsOfTheGene, transcriptomeMinusShallowTranscriptomeDict, novelPath, logFile):
	for transcriptId in transcriptIdsOfTheGene:
		if transcriptId in transcriptomeMinusShallowTranscriptomeDict:
			transcript = transcriptomeMinusShallowTranscriptomeDict[transcriptId]
			#check first if we have an exact matching
			if novelPath.getSeq() in transcript.getSeq() or novelPath.getSeq().reverse_complement() in transcript.getSeq():
				return True

			#check if we have an inexact match with >X% of matches
			novelSeqLength = len(novelPath.getSeq())
			FWAlignment = pairwise2.align.localms(novelPath.getSeq(), transcript.getSeq(), 1, 0, -10000, -10000, one_alignment_only=True)[0]
			RCAlignment = pairwise2.align.localms(novelPath.getSeq().reverse_complement(), transcript.getSeq(), 1, 0, -10000, -10000, one_alignment_only=True)[0]
			threshold=0.95
			if FWAlignment[2]>=threshold*novelSeqLength:
				logFile.write("Inexact match:\n")
				logFile.write(format_alignment(*FWAlignment)+"\n")
				return True
			if RCAlignment[2]>=threshold*novelSeqLength:
				logFile.write("Inexact match:\n")
				logFile.write(format_alignment(*RCAlignment)+"\n")
				return True

	return False


def novelPathIsInAnywhereInTheTranscriptome(fullTranscriptome, knownPathGene, novelPath, bubble, transcript2Gene):
	possibleGenes=[]
	for transcript in fullTranscriptome:
		if novelPath.getSeq() in transcript.getSeq() or novelPath.getSeq().reverse_complement() in transcript.getSeq():
				possibleGenes.append(transcript2Gene[transcript.getHeader()])
	return len(possibleGenes)>0

#saveFile(ASTALAVISTABubblesToBeFoundFile, "ASTALAVISTABubblesToBeFoundFile")
def saveFile(bubbles, filename):
	with open(filename, "w") as file:
		for bubble in bubbles:
			file.write(str(bubble)+"\n")
	print "Check file %s"%filename

import argparse
if __name__=="__main__":
	#parse the args
	parser = argparse.ArgumentParser(description='Checks EYTA performance')
	parser.add_argument('--short_reads', dest="sr", required=True, help='Short reads generated with wgsim')
	parser.add_argument('--transcriptomeFasta', dest="trFasta", required=True, help='Full transcriptome in Fasta Format')
	parser.add_argument('--transcriptomeGTF', dest="trGTF", required=True, help='Full transcriptome in GTF Format')
	parser.add_argument('--shallow_tr', dest="shallowTr", required=True, help='Shallow transcriptome generated with generate_shallow_transcriptome.py')
	parser.add_argument('--ASTALAVISTA_AS_bubbles', dest="ASTALAVISTA_AS_bubbles", required=True, help='ASTALAVISTA bubbles as FASTA files')
	parser.add_argument('--EYTA_AS_bubbles', dest="EYTA_AS_bubbles", required=True, help='*EYTA.alternative_paths.detailed.alternative_splicing')
	args = parser.parse_args()

	#get the shallow transcriptome
	shallowTranscriptome = getTranscriptome(args.shallowTr)

	#get transcriptomeMinusShallowTranscriptomeDict
	fullTranscriptome = getTranscriptome(args.trFasta)
	transcriptomeMinusShallowTranscriptomeDict = {transcript.getHeader(): transcript for transcript in fullTranscriptome if transcript not in shallowTranscriptome}

	#get gene2Transcripts and transcript2Gene dicts
	gene2Transcripts, transcript2Gene = parseGTF(args.trGTF)

	#get the transcriptsExpressedInShortReads
	transcriptsExpressedInShortReads = getTranscriptsExpressedInShortReads(args.sr)


	#get the ASTALAVISTA bubbles that have to be found by EYTA
	ASTALAVISTABubblesToBeFound = getASTALAVISTABubblesToBeFound(transcriptsExpressedInShortReads, shallowTranscriptome, args.ASTALAVISTA_AS_bubbles)
	print "%d ASTALAVISTA bubbles to be found by EYTA..."%len(ASTALAVISTABubblesToBeFound)

	#core function that get the precision and recall
	ASTALAVISTABubblesFoundByEYTA, ASTALAVISTABubblesNotFoundByEYTA, EYTABubblesThatMapToNoASTALAVISTABubbles, \
	EYTABubblesThatMapToExactlyOneASTALAVISTABubble, EYTABubblesThatMapToMoreThanOneASTALAVISTABubbles, EYTA_AS_bubbles = \
	getASTALAVISTABubblesFoundByEYTA(ASTALAVISTABubblesToBeFound, args.EYTA_AS_bubbles)


	#however, the precision is not good probably due to some ASTALAVISTA bubbles, maybe our own method bugs, or the fact that bubbles in DBGs are more general than bubbles in splicing graphs
	#here, we try to understand better the precision:
	#We go through all the bubbles in EYTABubblesThatMapToNoASTALAVISTABubbles, and check if the novel path is in any transcript of the gene of
	#of the known path. The transcript considered should not be in the shallow transcriptome
	EYTABubblesThatMapToATranscriptOfTheSameGene=[]
	EYTABubblesNeedingExplanation=[]
	with open("EYTABubblesThatMapToATranscriptOfTheSameGene.log", "w") as logFile:
		for bubbleIndex, bubble in enumerate(EYTABubblesThatMapToNoASTALAVISTABubbles):
			print >> sys.stderr, "Checking EYTABubblesThatMapToATranscriptOfTheSameGene - %.2f"%(float(bubbleIndex)/len(EYTABubblesThatMapToNoASTALAVISTABubbles)*100)
			if bubble.getUpperPath().isNovel() and \
				novelPathIsInTranscriptomeMinusShallowTranscriptome(gene2Transcripts[transcript2Gene[bubble.getLowerPath().getTranscript()]], transcriptomeMinusShallowTranscriptomeDict, bubble.getUpperPath(), logFile):
				EYTABubblesThatMapToATranscriptOfTheSameGene.append(bubble)
			elif bubble.getLowerPath().isNovel() and \
				novelPathIsInTranscriptomeMinusShallowTranscriptome(gene2Transcripts[transcript2Gene[bubble.getUpperPath().getTranscript()]], transcriptomeMinusShallowTranscriptomeDict, bubble.getLowerPath(), logFile):
				EYTABubblesThatMapToATranscriptOfTheSameGene.append(bubble)
			else:
				EYTABubblesNeedingExplanation.append(bubble)


	print >> sys.stderr, "Checking EYTABubblesChimeric..."
	EYTABubblesChimeric=[]
	EYTABubblesThatIReallyCantExplain=[]
	for i, bubble in enumerate(EYTABubblesNeedingExplanation):
		if bubble.getUpperPath().isNovel() and \
			novelPathIsInAnywhereInTheTranscriptome(fullTranscriptome, transcript2Gene[bubble.getLowerPath().getTranscript()], bubble.getUpperPath(), bubble, transcript2Gene):
			EYTABubblesChimeric.append(bubble)
		elif bubble.getLowerPath().isNovel() and \
			novelPathIsInAnywhereInTheTranscriptome(fullTranscriptome, transcript2Gene[bubble.getUpperPath().getTranscript()], bubble.getLowerPath(), bubble, transcript2Gene):
			EYTABubblesChimeric.append(bubble)
		else:
			EYTABubblesThatIReallyCantExplain.append(bubble)



	print "%d ASTALAVISTA bubbles found by EYTA..." % len(ASTALAVISTABubblesFoundByEYTA)
	print "%d ASTALAVISTA bubbles NOT found by EYTA..." % len(ASTALAVISTABubblesNotFoundByEYTA)
	print "%d EYTA bubbles mapping to exactly 1 ASTALAVISTA bubble..." % len(EYTABubblesThatMapToExactlyOneASTALAVISTABubble)
	print "%d EYTA bubbles mapping to >1 ASTALAVISTA bubble..." % len(EYTABubblesThatMapToMoreThanOneASTALAVISTABubbles)
	print "%d EYTA bubbles that do not map to any ASTALAVISTA bubble..." % len(EYTABubblesThatMapToNoASTALAVISTABubbles)
	print "From the above, %d EYTA bubbles map to another transcript of the SAME GENE of the shallow transcript (probably ASTALAVISTA error)..." % len(EYTABubblesThatMapToATranscriptOfTheSameGene)
	print "From the above, %d EYTA bubbles are chimeric (map to ANOTHER GENE)..." % len(EYTABubblesChimeric)
	print "From the above, %d EYTA bubbles that I cant explain..." % len(EYTABubblesThatIReallyCantExplain)


	print "Precision: %f"%(float(len(EYTABubblesThatMapToExactlyOneASTALAVISTABubble)+len(EYTABubblesThatMapToMoreThanOneASTALAVISTABubbles)+len(EYTABubblesThatMapToATranscriptOfTheSameGene)) / \
							 (len(EYTA_AS_bubbles)))
	print "Recall: %f" % (float(len(ASTALAVISTABubblesFoundByEYTA)) / len(ASTALAVISTABubblesToBeFound))


	saveFile(ASTALAVISTABubblesNotFoundByEYTA, "ASTALAVISTABubblesNotFoundByEYTA")
	saveFile(EYTABubblesChimeric, "EYTABubblesChimeric")
	saveFile(EYTABubblesThatIReallyCantExplain, "EYTABubblesThatIReallyCantExplain")
	saveFile(ASTALAVISTABubblesFoundByEYTA, "ASTALAVISTABubblesFoundByEYTA")
	saveFile(ASTALAVISTABubblesToBeFound, "ASTALAVISTABubblesToBeFound")
