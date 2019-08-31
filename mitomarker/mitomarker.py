#!/usr/bin/python
# -*- coding: utf-8 -*-
# Created by Natasha Picciani on Aug-23-2019 using Python 3.7
# Contact natasha.picciani@gmail with questions/comments

"""Use result files from MitoFinder to tag (append MT_) the mitochondrial
genes in a GTF file."""

import glob2
import re
from pathlib import PurePosixPath
import argparse

def concatenateFiles(type, target_directory, output_directory):
	"""Concatenate the files of a type in a target directory and return the path to concatenated file.

	Keyword arguments:
	type -- the file type to be concatenated
	target_directory -- path to directory with the targeted files
	output_directory -- path to directory where the concatenated file should be placed

	Note: Mitofinder produces several gff files with detailed info for each mitochondrial gene found
	but does not produce a final concatenated file with the info for all of them

	"""

	filenames = glob2.glob(target_directory + "/*." + type) # list of all files of a type in the directory
	outfile = output_directory + "/concatenated." + type
	with open(outfile, 'w') as f:
		for file in filenames:
			with open(file) as infile:
				f.write(infile.read())
	return outfile

def mitochondrialGenes(concatenated_gff):
	"""List mitochondrial genes in the concatenated gff file.

	Keyword arguments:
	concatenated_gff -- path to concatenated gff file
	"""
	mitochondrialList=[]
	with open(concatenated_gff, 'r') as infile:
		for line in infile:
			line=line.strip('\n')
			elementList=line.split('\t')
			mitochondrialList.append(elementList[0])
	return mitochondrialList

def mitochondrialGTF(concatenated_gff, output_directory):
	"""Convert GFF file with information for mitochondrial genes to GTF and return the path to GTF file.

	Keyword arguments:
	concatenated_gff -- path to concatenated GFF file
	output_directory -- path to directory where GTF file should be placed

	Important: transcript name must follow 'compXX_cXX_seqXX' 

	Notes: Make a separate GTF file with tagged mitochondrial genes (MT_) so that later 
	it is merged it with the post-processed original GTF file (after removal of mitochondrial genes).

	"""

	output = output_directory + "/concatenated.gtf"
	with open(concatenated_gff, 'r') as infile:
		with open(output, 'w') as outfile:
			for line in infile:
				line=line.strip('\n')
				elementList=line.split('\t')
				p=re.search(r'((comp\d+_c\d+)_seq\d+)', elementList[0])
				geneID=p.group(2)
				transcriptID=p.group(1)
				outfile.write(elementList[0]+'\tmitofinder\tgene\t'+elementList[3]+'\t'+elementList[4]+'\t'+elementList[5]+'\t'+elementList[6]+'\t'+elementList[7]+'\tgene_id "'+geneID+'"; transcript_id "'+transcriptID+'"; gene_name "'+'MT_'+elementList[0]+'_'+elementList[8]+'";'+'\n')    
	return output

def cleanGTF(originalGTF, mitochondrialList, output_directory):
	""" Remove lines with mitochondrial genes from original GTF file and returns path to the new clean GTF file.

	Keyword arguments:
	originalGTF -- path to original GTF file
	mitochondrialList -- list of mitochondrial genes to search for
	output_directory -- path to directory where clean GTF file should be placed

	Note: A single transcript can contain several mitochondrial gene products as per annotation from Mitofinder. 
	Instead of editing those lines, best to just remove them and add the mitochondrial gtf block we generated above.

	"""
	originalGTF_basename = PurePosixPath(originalGTF).stem
	output = output_directory + "/" + originalGTF_basename + ".clean.gtf"

	with open(originalGTF, 'r') as infile:
		with open (output, 'w') as outfile:
			for line in infile:
				if not line.startswith(tuple(mitochondrialList)):
					outfile.write(line)
	return output

def concatenateGTF(cleanGTF, mitochondrialGTF, output_directory):
	"""Concatenate clean and mitochondrial GTF files and return the path to concatenated file.

	Keyword arguments:
	cleanGTF -- path to clean GTF file
	mitochondrialGTF -- path to mitochondrial GTF file
	output_directory -- path to directory where concatenated gtf file should be placed

	"""
	originalGTF_basename = PurePosixPath(originalGTF).stem
	finalGTF = output_directory + "/" + originalGTF_basename + ".final.gtf"

	with open(cleanGTF,'r') as clean:
		with open(mitochondrialGTF,'r') as mitochondrial:
			with open(finalGTF,'w') as final:
				for line in clean:
					final.write(line)
				for line in mitochondrial:
					final.write(line)
	return finalGTF

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Label the mitochondrial genes found by MitoFinder in a reference GTF file')

	parser.add_argument('-o', metavar = 'output_directory', type = str, required = True,
					   default = ".", help = "directory for placing output files")
	parser.add_argument('-m', metavar = 'mitofinder_folder', type = str, required = True, 
						help = "folder with mitofinder final results")
	parser.add_argument('-gtf', metavar = 'GTF_file', type = str, required = True,
					   help = "original GTF file with mitochondrial genes")
	parser.add_argument('-version', action='version', version='%(prog)s 0.1.0')

	args = parser.parse_args()

	outDir = args.o
	mitofolder = args.m
	originalGTF = args.gtf
			  
	# Concatenate GFF files produced by MitoFinder
	concatenated_gff = concatenateFiles(type="gff",target_directory=mitofolder, output_directory=outDir)
	# List mitochondrial genes in gff file.
	mitoList = mitochondrialGenes(concatenated_gff)
	# Produce a small GTF file with the mitochondrial genes only based on MitoFinder output
	mitoGTF = mitochondrialGTF(concatenated_gff, output_directory=outDir)
	# Remove lines with mitochondrial genes from original GTF file 
	cleanGTF_file = cleanGTF(originalGTF, mitoList, output_directory=outDir)
	# Concatenate the mitochondrial GTF with clean GTF 
	concatGTF_file = concatenateGTF(cleanGTF_file, mitoGTF, output_directory=outDir)