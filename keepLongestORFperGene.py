#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created by Natasha Picciani on Oct-14-2019 using Python 3.7
# Contact natasha.picciani@gmail.com with questions/comments

"""
Read a FASTA file with open reading frames (ORFs), and keep the longest
ORF for each gene. Return a fasta file with longest ORFs per gene.

"""

import re
import argparse
from pathlib import PurePosixPath
from Bio import SeqIO


def getGeneID(header, type):
	
	"""
	
	Search the header of a fasta file and return the gene identifier of a sequence."

	Arguments:
	header -- the header string (the record.id value of a sequence when using SeqIO package).
	type -- the type of gene identifier; type_1 = "compXX_cXX_seqXX", type_2="TRINITY_DNXXXXX_cX_gX_iX",
			type_3="SegXX.XX.XXXX".

	"""

	if type == "type_1":
		searchStr='(comp\d+_c\d+)_seq\d+'
	if type == "type_2":
		searchStr='(TRINITY_DN\d+_c\d+_g\d+)_i\d+'
	if type == "type_3":
		searchStr='(Seg\d+)\..+'
	parts=header.split('\t')
	p=re.match(searchStr, parts[0])
	geneID=p.group(1)
	return geneID

def keepLongest(fastaFile, identifier_type, output_directory):

	"""	
	
	Read a fasta file with ORFs (open reading frames) and select the longest ORF per each gene. Return a FASTA file.

	Arguments:
	fastaFile -- FASTA file with ORF sequences
	identifier_type -- the type of gene identifier; type_1 = "compXX_cXX_seqXX", type_2="TRINITY_DNXXXXX_cX_gX_iX",
			type_3="SegXX.XX.XXXX".
	output_directory -- path to directory for placing output file

	"""
	
	seqs={}
	currentGeneID=''
	longestORFperGene = output_directory + "/" + PurePosixPath(fastaFile).stem + "_longestORFperGene.fasta"
# 	unique = 0 #uncomment for testing
# 	duplicates = 0 #uncomment for testing
	
	with open(longestORFperGene,"w") as outfile:

		for record in SeqIO.parse(fastaFile,"fasta"): #find gene names
			geneID = getGeneID(record.id, identifier_type)
			if currentGeneID == '': #start the loop with an empty GeneID string
				currentGeneID = geneID
				seqs[record.description]=[geneID,record.seq] #add first sequence to dictionary "seqs", which has geneID as a list value
			else:
				if geneID != currentGeneID: #if the previous geneID is different from that of the current record
# 					unique += 1 #uncomment for testing
					seqs[record.description]=[geneID,record.seq] #add unique sequences to dictionary "seqs"
					currentGeneID = geneID #re-set current geneID     
				else: #if the geneID is the same as that of the previous record				
# 					duplicates += 1 #uncomment for testing
					for key, value in seqs.items(): #search for the previous record in the dictionary "seqs"
						if value[0] == currentGeneID: #and save its sequence length and name
							dupeLength=(len(value[1]))
							dupeID=key
					if dupeLength<len(record.seq): #if length of previous record is smaller that of than current record
						del seqs[dupeID] #get rid of previous record
						seqs[record.description]=[geneID,record.seq] #add current record to dictionary
						currentGeneID = geneID
 
		for name, seq in seqs.items(): #print sequences to an output file
			print (">" + name, file=outfile)
			print (seq[1], file=outfile)  


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Read a FASTA file with open reading frames(ORFs), and keep the longest ORF per gene')

	parser.add_argument('-f', metavar = 'FASTA_file', type = str, required = True,
					   help = "FASTA file with ORFs")
	parser.add_argument('-o', metavar = 'output_directory', type = str, required = True,
					   default = ".", help = "directory for placing output file")
	parser.add_argument('-identifier', metavar = 'identifier_type', type = str, choices=['type_1', 'type_2', 'type_3'], required = False,
						default = "type_1", help = "type of gene identifier in the ORF FASTA file")
	parser.add_argument('-version', action='version', version='1.0')

	args = parser.parse_args()

	
	outdir = args.o
	fastaFile = args.f
	idType = args.identifier
			  
	# Select longest ORFs per gene
	keepLongest(fastaFile, idType, output_directory=outdir)


