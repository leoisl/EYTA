nonParalogousGenesIDs=[]
with open("Homo_sapiens.GRCh38.90.non_paralogous_genes.txt") as fin:
	for i, line in enumerate(fin):
		if i==0:
			continue
		nonParalogousGenesIDs.append(line.strip())

with open("Homo_sapiens.GRCh38.90.gtf") as fin:
	for line in fin:
		line=line.strip()
		if line[0]=="#":
			print line
		else:
			lineSplit = line.split()
			geneId = lineSplit[9][1:-2]
			if geneId in nonParalogousGenesIDs:
				print line