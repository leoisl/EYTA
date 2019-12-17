#fix the headers because when converting to BAM, we get a bug when the header is too large
import sys

with open (sys.argv[1]) as fin, open (sys.argv[1]+".fixed_header.fa", "w") as fout:
		for i, line in enumerate(fin):
			line = line.strip()
			if line.startswith(">Bubble"):
				lineSplit = line.split("_")
				fout.write("%s_%s_%d\n"%(lineSplit[0], lineSplit[1], i/2))
			else:
				fout.write("%s\n"%line)