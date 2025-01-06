#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Created on Aug 21, 2020
by Natasha Picciani
Python 3.7
usage: getSequences_onefasta.py [transcriptList] [transcriptome] [outDir]

'''

import sys
import re
from pathlib import PurePosixPath
from Bio import SeqIO
from Bio import AlignIO

if len(sys.argv)<2:
	sys.stderr.write(__doc__)
	sys.stderr.write("Please provide a list of transcripts and a transcriptome".format(SearchStr))

else:
	transcriptList= sys.argv[1]
	transcriptomeFile=sys.argv[2]
	outputDir=sys.argv[3]

	transcriptomeName = PurePosixPath(transcriptomeFile).stem

	querylist = []
	with open(transcriptList,'r') as file:
		for line in file:
			line=line.strip('\n')
			querylist.extend(line.split(' '))

	transcriptome = open(transcriptomeFile, 'r')
	SequencesFile = outputDir + "/" + transcriptomeName + ".selected.fasta"

	Sequences =[]
	for record in SeqIO.parse(transcriptome,'fasta'):
		for item in querylist:
			if item.strip('\n') == record.id:
				Sequences.append(record)
			
	SeqIO.write(Sequences,SequencesFile,'fasta')
	transcriptome.close()
