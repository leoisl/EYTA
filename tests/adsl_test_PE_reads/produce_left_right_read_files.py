#python produce_left_right_read_files.py left_reads_ids /pandata/benoit/MCF7_sipp/trimmed_data/fastq/siCTL_N1_GGCTAC_R1_trim_right_cutadapt_match_trim_to_100_match.fastq.gz /pandata/benoit/MCF7_sipp/trimmed_data/fastq/siCTL_N1_GGCTAC_R2_trim_right_cutadapt_match_trim_to_100_match.fastq.gz
import sys
import gzip
import itertools

with open(sys.argv[1]) as left_reads_ids_file, gzip.open(sys.argv[2]) as left_reads, gzip.open(sys.argv[3]) as right_reads, open("left_reads.fq", "w") as left_reads_out, open("right_reads.fq", "w") as right_reads_out:
	left_reads_ids=[line.rstrip() for line in left_reads_ids_file]
	doIPrint=0
	i=-1
	for left_line, right_line in itertools.izip(left_reads, right_reads):
		i+=1
		#just print
		if doIPrint>0:
			left_reads_out.write(left_line)
			right_reads_out.write(right_line)
			doIPrint-=1
			continue

		#check header
		left_line=left_line.rstrip()
		right_line=right_line.rstrip()
		if i%4==0:
			#print "Checking read:", left_line
			read_id = left_line[1:].split()[0]
			if read_id in left_reads_ids:
				print "Added:", read_id
				left_reads_out.write(left_line+"\n")
				right_reads_out.write(right_line+"\n")
				doIPrint=3
			#else:
			#	print "Pass!"




