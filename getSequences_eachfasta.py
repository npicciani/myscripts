#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Created on Oct 1, 2019
by Natasha Picciani

'''

import re
import glob2
import os
from Bio import SeqIO
from pathlib import PurePosixPath


def makeFolder(folderPath):

	"""Create a folder.

	Keyword arguments:
	folderPath -- path to new directory

	"""
	try:
		os.mkdir(folderPath)
	except OSError:
		print("Failed making directory" + folderPath)

def getSequences(transcriptome, transcripts, outdir):

	"""Pull sequences from a transcriptome based on a list of transcripts. 
	Return a separate fasta file for each transcript.

	Keyword arguments:
	transcriptome -- path to transcriptome assembly from which sequences will be retrieved
	transcripts -- path to list of transcripts to be retrieved
	outdir -- path to directory for output fasta files

	"""

	transcriptList=[]

	with open(transcripts) as infile:
		for line in infile:
			line=line.strip('\n')
			transcriptList.append(line)

	sequences =[]
	with open(transcriptome) as file:	
		for record in SeqIO.parse(file,'fasta'):
			for item in transcriptList:
				if item.strip('\n') == record.id:
					sequences.append(record)
	
	for item in sequences:
		outfile= outdir + "/" + item.id + ".fasta"
		SeqIO.write(item, outfile, 'fasta')


# User defined paths
rootdir='/home/picciani/local/datasets/hydra/cell_markers_aep/'
transcriptome='/home/picciani/local/datasets/hydra/hydra_aep_transcriptome/aepLRv2.fa'
filenames = glob2.glob(rootdir + "/*.transcripts") # list files with transcript list in root directory

for file in filenames:
	folderPath = rootdir + "/" + PurePosixPath(file).stem
	makeFolder(folderPath) # make folders for each transcript list file
	getSequences(transcriptome,file,folderPath) # generate fasta files from each list of transcripts
 
