#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Created on Oct 4, 2019
by Natasha Picciani
usage: addTag.py assembly tag outdir

[assembly] : path to transcriptome assembly
[tag] : tag to add to each header
[outdir] : directory for output file

'''

import sys
from pathlib import PurePosixPath

# Add a name tag to transcript headers
def addTag(transcriptome):
	with open(transcriptome) as t:
		with open(outfile, 'w') as out:
			for line in t:
				line=line.strip('\n')
				if line.startswith('>'):
					out.write(line + '_' + tag + '\n')
				else:
					out.write(line + '\n')
				

if len(sys.argv)<2:
	sys.stderr.write(__doc__)
	sys.stderr.write("Please provide arguments")


else:
    
	transcriptome= sys.argv[1] #short headers are recommended for the starting fastafile; e.g ">comp0_c0_seq1"
	tag= sys.argv[2]
	outdir=sys.argv[3]
	
    
	outname = PurePosixPath(transcriptome).stem
	outfile = outdir + '/' + outname + '.tagged'
	addTag(transcriptome)
