import os

gtfFilename="Homo_sapiens.GRCh38.90.gtf"
outputFolder="%s_divided_by_gene_id"%gtfFilename

#make the output folder
os.makedirs(outputFolder)

with open(gtfFilename) as gtfFile:
	for line in gtfFile:
		if line[0]=="#":
			continue

		lineSplit = line.strip().split()
		geneId = lineSplit[9][1:-2]
		with open("%s/%s.gtf"%(outputFolder, geneId), "a") as geneGtfFile:
			geneGtfFile.write(line)

