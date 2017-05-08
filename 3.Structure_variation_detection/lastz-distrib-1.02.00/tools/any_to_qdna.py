#!/usr/bin/env python
"""
Convert any file to a LASTZ quantum dna file, just by appending qdna headers
----------------------------------------------------------------------------

Qdna file format is shown below (omitting "named properties, which we don't
use).  We simply create all the headers and copy the file as the "data
sequence".
	
	offset 0x00: C4 B4 71 97   big endian magic number (97 71 B4 C4 => little endian)
	offset 0x04: 00 00 02 00   version 2.0 (fourth byte is sub version)
	offset 0x08: 00 00 00 14   header length (in bytes, including this field)
	offset 0x0C: xx xx xx xx   S, offset (from file start) to data sequence
	offset 0x10: xx xx xx xx   N, offset to name, 0 indicates no name
	offset 0x14: xx xx xx xx   length of data sequence (counted in 'items')
	offset 0x18: 00 00 00 00   (offset to named properties, not used)
	offset    N: ...           name (zero-terminated string)
	offset    S: ...           data sequence

Alternatively, the option --simple creates an old-style qdna file, which
consists of a different magic number prepended to the data sequence.

:Author: Bob Harris (rsharris@bx.psu.edu)
"""

import sys
import os


def main():
	global debug

	qdnaOldMagic = 0xF656659EL	# big endian magic number for older qdna files
	qdnaMagic    = 0xC4B47197L	# big endian magic number for qdna files
	qdnaVersion  = 0x00000200L 

	##########
	# parse the command line
	##########

	debug  = []
	name   = None
	strip  = False
	simple = False

	# pick off options

	args = sys.argv[1:]
	while (len(args) > 0):
		arg = args.pop(0)
		val = None
		fields = arg.split("=",1)
		if (len(fields) == 2):
			arg = fields[0]
			val = fields[1]
			if (val == ""):
				assert (False), "missing a value in %s=" % arg

		if (arg == "--name") and (val != None):
			name = val
		elif (arg in ["--striplinebreaks","--strip"]) and (val == None):
			strip = True
		elif (arg in ["--simple","--old"]) and (val == None):
			simple = True
		elif (arg == "--debug") and (val != None):
			debug.append(val)
		elif (arg.startswith("--")) or (val != None):
			assert (False), "unknown argument: %s" % arg
		else:
			assert (False), "unknown argument: %s" % arg

	assert (not simple) or (name == None), \
	      "simple qdna file cannot carry a sequence name"

	##########
	# read the fasta file
	##########

	seq = []
	for line in sys.stdin:
		if (strip): line = line.rstrip()
		seq += [line]
	seq = "".join(seq)

	##########
	# write the qdna file
	##########

	if (not simple):
		headerLen = 20
		if (name == None):
			nameOffset = 0
			seqOffset  = headerLen + 8;
		else:
			nameOffset = headerLen + 8;
			seqOffset  = nameOffset + len(name) + 1

	# prepend magic number

	if (simple):
		write_4(sys.stdout,qdnaOldMagic)
	else:
		write_4(sys.stdout,qdnaMagic)

	# write the rest of the header

	if (not simple):
		write_4(sys.stdout,qdnaVersion)
		write_4(sys.stdout,headerLen)
		write_4(sys.stdout,seqOffset)
		write_4(sys.stdout,nameOffset)
		write_4(sys.stdout,len(seq))
		write_4(sys.stdout,0)

		if (name != None):
			sys.stdout.write(name)
			sys.stdout.write(chr(0))

	# write the sequence

	sys.stdout.write(seq)


def write_4(f,val):
	f.write (chr((val >> 24) & 0xFF))
	f.write (chr((val >> 16) & 0xFF))
	f.write (chr((val >>  8) & 0xFF))
	f.write (chr( val        & 0xFF))


if __name__ == "__main__": main()
