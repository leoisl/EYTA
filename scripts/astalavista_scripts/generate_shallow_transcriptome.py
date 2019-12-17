#print the fasta with the shallow transcriptome
from Bio import SeqIO
def getTranscript2Sequence(transcritomeFASTAFilename, minLength):
	"""
	Returns transcripts2Sequence, sequences must have size >= minLength
	"""
	transcript2Sequence={}

	for record in SeqIO.parse(transcritomeFASTAFilename, "fasta"):
		if len(record.seq)>=minLength:
			transcript2Sequence[record.id]=record.seq
		else:
			print "skip sequence '%s' as it is shorter than %d!"%(record.id, minLength)

	return transcript2Sequence


def getGene2TranscriptsFromGTF(gtfFilename, transcript2Sequence):
	"""
	returns a dict where key is gene ID and value is a list of transcript IDs
	the transcripts should be the valid ones, i.e. in transcript2Sequence
	"""
	gene2Transcripts={}
	with open(gtfFilename) as gtfFile:
		for line in gtfFile:
			if line[0]=="#":
				continue
			lineSplit = line.split()
			feature = lineSplit[2]

			if feature=="transcript":
				geneId = lineSplit[lineSplit.index("gene_id")+1][1:-2] #[1:-2] is to remove the begin and the end stuff that we don't care, i.e. gene_id "ENSG00000279493";
				transcriptId = lineSplit[lineSplit.index("transcript_id")+1][1:-2] #[1:-2] is to remove the begin and the end stuff that we don't care, i.e. transcript_id "ENST00000456328";

				#check if the transcript is valid
				if transcriptId in transcript2Sequence:
					#add it
					if geneId not in gene2Transcripts:
						gene2Transcripts[geneId]=[]
					gene2Transcripts[geneId].append(transcriptId)
	return gene2Transcripts

import random
def simulateShallowTranscriptome(gene2Transcripts, proportion):
	gene2ShallowTranscripts={}
	for gene, transcripts in gene2Transcripts.items():
		nbOfTranscripts = max(1, int(proportion*len(transcripts))) #nb of transcript in the shallow transcriptome is at least 1
		gene2ShallowTranscripts[gene] = random.sample(transcripts, nbOfTranscripts)
	return gene2ShallowTranscripts

#print the fasta with the shallow transcriptome
def printFastaShallowTranscriptome(gene2ShallowTranscripts, transcript2Sequence, outputFilename):
	with open(outputFilename, "w") as outputFile:
		for gene, transcripts in gene2ShallowTranscripts.items():
			for transcript in transcripts:
				outputFile.write(">%s\n%s\n"%(transcript, transcript2Sequence[transcript]))

import argparse
if __name__=="__main__":
	#parse the args
	parser = argparse.ArgumentParser(description='Generate shallow transcriptome')
	parser.add_argument('--gtf', dest="transcritomeGTF", required=True, help='The GTF file describing the transcriptome')
	parser.add_argument('--gtfFasta', dest="transcritomeFASTA", required=True, help='The transcripts as fasta (generated with gffread)')
	parser.add_argument('--proportion', dest="proportion", required=True, type=float, help='proportion of transcripts to be simulated from a gene')
	parser.add_argument('--min_length_transcript', dest="minLength", required=True, type=int, help='Only transcripts with length larger than this will be considered. This is to cope with small transcripts - if you simulated reads with 150 bp, you cant simulate transcripts smaller than this, so this should be set to 150 (otherwise, the shallow transcriptome will have a transcript not sequenced by short reads)')
	parser.add_argument('--output', dest="outputFilename", help='output file', default="shallow_transcriptome.fa")
	args = parser.parse_args()

	#get transcript2Sequence
	transcript2Sequence = getTranscript2Sequence(args.transcritomeFASTA, args.minLength)

	#get gene2Transcripts
	gene2Transcripts = getGene2TranscriptsFromGTF(args.transcritomeGTF, transcript2Sequence)

	#simulate the shallow transcriptome
	gene2ShallowTranscripts = simulateShallowTranscriptome(gene2Transcripts, args.proportion)

	#print the fasta with the shallow transcriptome
	printFastaShallowTranscriptome(gene2ShallowTranscripts, transcript2Sequence, args.outputFilename)