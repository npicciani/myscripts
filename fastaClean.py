#!/usr/bin/env python
"""
fastaClean.py - version 1.0

developed by Steve Haddock at http://practicalcomputing.org/node/127
modified by Natasha Picciani
 
Usage:
  fastaClean.py infile.fasta > outfile.fasta
"""

import re
import sys
  
SearchStr= r'>.+bsT_(comp.+seq\d+),.+'
Replacement = r'>\1'


if len(sys.argv)<2:
	sys.stderr.write(__doc__)
	sys.stderr.write("Current search string: '{}'\n".format(SearchStr))

else:
	fastaFileName= sys.argv[1]
	fastaFile= open(fastaFileName, 'rU')
	p = re.compile(SearchStr)


	for line in fastaFile:
		line=line.strip('\n')
		if line.startswith('>'):
			newline = p.sub(Replacement,line)
			print (newline)
		else:
			print (line)

fastaFile.close()