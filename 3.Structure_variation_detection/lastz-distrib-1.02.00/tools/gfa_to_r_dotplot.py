#!/usr/bin/env python
"""
Program to convert a GFA file to a file suitable for dot-plot in R
------------------------------------------------------------------

:Author: Bob Harris (rsharris@bx.psu.edu)

"""

import sys,re


segPatt  = "^a +(?P<pos1>[0-9]+)(?P<strand1>[+-])/(?P<pos2>[0-9]+)(?P<strand2>[+-]) +(?P<length>[0-9]+) +(?P<score>[-0-9.]+)"
segRe    = re.compile(segPatt)
seg2Patt = "^a +(?P<pos1>[0-9]+)/(?P<pos2>[0-9]+) +(?P<length>[0-9]+) +(?P<score>[-0-9.]+)"
seg2Re   = re.compile(seg2Patt)


def usage(s=None):
	message = """
gfa_to_r_dotplot < gfa_file > r_data_file
"""

	if (s == None): sys.exit (message)
	else:           sys.exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	if (len(sys.argv) != 1):
		usage("give me no arguments")

	# process the stats

	first = True

	seq1Start = None

	idNum = 0
	for line in sys.stdin:
		line = line.rstrip()

		if (line.startswith("s ")):
			(info1,info2) = parse_s_record(line)
			(filename1,start1,end1,strand1,contig1) = info1
			(filename2,start2,end2,strand2,contig2) = info2

			if (seq1Start == None):
				assert (strand1+strand2 == "++"), "first s-record has to be ++"
				(seq1Name,seq1Start,seq1End,seq1Contig) \
				    = (filename1,start1,end1,contig1)
				(seq2Name,seq2Start,seq2End,seq2Contig) \
				    = (filename2,start2,end2,contig2)
				seq1Len = seq1End+1 - seq1Start
				seq2Len = seq2End+1 - seq2Start
			else:
				assert (start1  == seq1Start)   \
				   and (end1    == seq1End)     \
				   and (contig1 == seq1Contig)  \
				   and (start2  == seq2Start)   \
				   and (end2    == seq2End)     \
				   and (contig2 == seq2Contig), \
				      "s-records have to have same sequence lengths"

		elif (line.startswith("A ")):
			assert (seq1Start != None), "missing s-record"
			idNum += 1

		elif (line.startswith("a ")):
			assert (seq1Start != None), "missing s-record"
			(b1,s1,b2,s2,length,score) = parse_a_record(line)
			if (s1 == "+"): (b1,e1) = (b1,b1+length-1)
			else:           (b1,e1) = (seq1Len+2-b1-length,seq1Len+1-b1)
			if (s2 == "+"): (b2,e2) = (b2,b2+length-1)
			else:           (b2,e2) = (seq2Len+2-b2-length,seq2Len+1-b2)
			if (first): first = False
			else:       print "NA NA %s %s NA" % (score,s1+s2)
			print "%s %s %s %s %s" % (seq1Start+b1,seq2Start+b2,score,s1+s2,idNum)
			print "%s %s %s %s %s" % (seq1Start+e1,seq2Start+e2,score,s1+s2,idNum)

	print >>sys.stderr, "dd = read.table(\"datafilename\",header=F)"
	print >>sys.stderr, "plot(dd[,1],dd[,2],type=\"l\","
	print >>sys.stderr, "     xlim=c(%s,%s),ylim=c(%s,%s)," \
	                  % (seq1Start,seq1End,seq2Start,seq2End)
	print >>sys.stderr, "     xlab=%s," % seq1Name
	print >>sys.stderr, "     ylab=%s)" % seq2Name


def parse_s_record(line):
	try:
		fields = line.split()
		if (len(fields) != 11): raise
		info1 = parse_s_info(fields[1:6])
		info2 = parse_s_info(fields[6:11])
	except ValueError:
		assert (False), "bad s-record: %s" % line
	return (info1,info2)


def parse_s_info(fields):
	filename = fields[0]
	start    = int(fields[1])
	end      = int(fields[2])
	contig   = int(fields[4])
	if   (fields[3] == "0"): strand = "+"
	elif (fields[3] == "1"): strand = "-"
	else:                    raise ValueError
	return (filename,start,end,strand,contig)


def parse_a_record(line):
	m = segRe.match(line)
	if (m == None): m = seg2Re.match(line)
	assert (m != None), "bad a-record: %s" % line

	pos1   = int(m.group("pos1"))
	pos2   = int(m.group("pos2"))
	length = int(m.group("length"))

	try:
		strand1 = m.group("strand1")
		strand2 = m.group("strand2")
	except:
		strand1 = strand2 = "+"

	try:
		score = int(m.group("score"))
	except:
		score = float(m.group("score"))

	return (pos1,strand1,pos2,strand2,length,score)


if __name__ == "__main__": main()
