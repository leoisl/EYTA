from intervaltree import Interval, IntervalTree
tree = IntervalTree()

with open ("Homo_sapiens.GRCh38.90.gtf") as fin:
	for line in fin:
		lineSplit =  line.strip().split()
		if lineSplit[2]=="gene":
			coord1 = int(lineSplit[3])
			coord2 = int(lineSplit[4])
			begin = min(coord1, coord2)
			end = max(coord1, coord2)+1
			print "begin: %d, end: %d, gene:%s"%(begin, end, lineSplit[13])
			if len(tree[begin:end])==0: #i.e. this gene is not overlapping with any other gene
				pass #ok
			else:
				print "OVERLAPPING!!!"

