"""
this works like a grep, but bubble-wise, not sequence wise
args:
python script.py <input file> <strings to match ... - several>
"""
import sys

with open(sys.argv[1]) as inputFile:
	lines=[line.strip() for line in inputFile]
bubbles = [ tuple(lines[i:i+4]) for i in xrange(0, len(lines), 4)]

patterns = sys.argv[2:]
for bubble in bubbles:
	present=True
	for pattern in patterns:
		if pattern not in bubble[0] and pattern not in bubble[2]:
			present=False

	if present:
		print "%s\n%s\n%s\n%s"%bubble
