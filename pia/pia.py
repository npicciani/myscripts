#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Created on May 6, 2019
by Natasha Picciani, modified PIA script (see Speiser et al. 2014)
usage: pia.py fastaFile outputDir
version 1.0


'''

import sys
import re
import subprocess
from ete3 import Tree
from pathlib import PurePosixPath
from Bio import SeqIO
from Bio import AlignIO

# Predicting the Open Reading Frames 
def make_ORFs(fastaFile):
	longORFs_log = subprocess.run([transdecoder_path, \
                    "-t", fastaFile, "-m", "70", "--output_dir", outputDir],\
                   check=True, capture_output=True)
	return longORFs_log

# Building a blast database
def make_blastdb(longest_orfs):
	blastdb_log = subprocess.run([makeblastdb_path,"-in",longest_orfs,\
                              "-dbtype","prot","-parse_seqids","-out",blastdbName], \
                             check=True, capture_output=True)
	return blastdb_log

# Blasting the database with an opsin baitsfile
def blastp(blastdbName):
	opsinCandidates_log = subprocess.run([blastp_path,"-db", blastdbName,\
                                      "-query", baits,"-out",opsinsBlastpFile,\
                                      "-evalue","1e-10","-outfmt","7"], check=True, capture_output=True)
	return opsinCandidates_log

# Adding candidate opsin sequences to a previous alignment
def mafft(opsinsFile):
	alignedOpsins = subprocess.check_output([mafft_path,"--add",opsinsFile,\
											"--reorder",originalAlignment], universal_newlines=True)
	return alignedOpsins

# Placing opsin sequences in a pre-calculated opsin tree
def raxml(opsinsPhylipFile):
	opsinTree = subprocess.run([raxml_path,"-f", "v",\
                                      "-s", opsinsPhylipFile,"-t",treefile,"-w",outputDir,\
                                      "-m","PROTGAMMAWAG","-T","8","-n", "Opsins"], check=True, capture_output=True)
	return opsinTree

# User-defined paths

transdecoder_path = "/home/picciani/local/bin/TransDecoder.LongOrfs" #path to Transdecoder.LongOrfs 
makeblastdb_path = "/home/picciani/local/bin/ncbi-blast-2.8.1+/bin/makeblastdb" #path to makeblastdb
blastp_path = "/home/picciani/local/bin/ncbi-blast-2.8.1+/bin/blastp" #path to blastp
mafft_path = "/home/picciani/local/bin/anaconda2/bin/mafft" #path to mafft
raxml_path = "/home/picciani/local/bin/anaconda2/bin/raxmlHPC-PTHREADS-SSE3" #path to raxmlHPC-PTHREADS-SSE3

baits= "/home/picciani/local/datasets/pia/baits1217.fasta" #baitfile to search for a gene family
originalAlignment = "/home/picciani/local/datasets/pia/picciani/all_0512_gt1_rs50_s65.aln" #alignment used to build the pre-calculated gene tree
treefile = "/home/picciani/local/datasets/pia/picciani/rep17.treefile" #pre-calculated gene tree

##############

if len(sys.argv)<2:
	sys.stderr.write(__doc__)
	sys.stderr.write("Current fasta file '{}'\n".format(fastaName))


else:
    
	fastaFile= sys.argv[1] #short headers are recommended for the starting fastafile; e.g ">comp0_c0_seq1"
	outputDir= sys.argv[2]
    
	fastaName = PurePosixPath(fastaFile).stem


# Using Transdecoder to predict the Open Reading Frames
	
	print ("Predicting open reading frames in " + fastaName)        
	longORFs_log = make_ORFs(fastaFile)
	longORFsFile = outputDir + "/longest_orfs.pep"
    
# Correcting the headers of the longest orfs file

	LineNumber=0

	newLongORFsFile = outputDir + "/longest_orfs.pep.shortheaders"
	newLongORFs =open(newLongORFsFile, 'w')

	with open (longORFsFile) as longorfs:
	
		for line in longorfs:
			line=line.strip('\n')
		
			if line.startswith('>'):
				line=line.split(' ')
				newline=line[0]
				newLongORFs.write(newline + '\n')
		
			else:
				newLongORFs.write(line + '\n')
	
		LineNumber += 1

# Building a blast database
	
	print ("Building blast database for " + fastaName)
	blastdbName = outputDir + "/" + fastaName + ".blastdb"
	blastdb_log = make_blastdb(newLongORFsFile)

# Blasting the database with an opsin baitsfile
	
	print ("Blasting the database with the opsin baitsfile")
	opsinsBlastpFile = outputDir + "/" + fastaName + ".opsins_blast_results"
	blastp_log = blastp(blastdbName)

# Exporting the second column from the blastp results file as a blastp hits list

	LineNumber = 0

	inFile = open(opsinsBlastpFile, 'r')
	opsinListInitial= outputDir + "/" + fastaName + ".candidate_opsins_seq_list_initial"
	outFile = open(opsinListInitial, "w")

	for line in inFile:
	
		if line[0] == "#":
			LineNumber= LineNumber + 1
			continue
		
		line=line.strip('\n')
	
		ElementList=line.split('\t')
		blasthits=ElementList[1]
	
		print (blasthits, file=outFile)
			 
		LineNumber= LineNumber + 1

	inFile.close()
	outFile.close()

# Removing duplicate lines in the initial blastp opsin list without sorting

	linesSeen = set()
	opsinList = outputDir + "/" + fastaName + ".candidate_opsins_seq_list"
	outfile = open(opsinList, 'w')
	for line in open(opsinListInitial, 'r'):
		if line not in linesSeen:
			outfile.write(line)
			linesSeen.add(line)
	outfile.close()

# Pulling the opsin sequences from the original ORF peptide file based on the blastp hits list

	print ("Pulling opsin sequences from ORF peptide file")
	my_list = []
	with open (opsinList) as opsinlist:
		
		for line in opsinlist:
			my_list.extend(line.split(' '))

	transcriptome = open(newLongORFsFile, 'r')
	opsinsFile = outputDir + "/" + fastaName + ".candidate_opsins.fasta"

	selected_transcripts = []
	for record in SeqIO.parse(transcriptome,'fasta'):
		for item in my_list:
			if item.strip() == record.id:
				selected_transcripts.append(record)
			
	SeqIO.write(selected_transcripts, opsinsFile,"fasta")
	transcriptome.close()
	
# Adding the opsin sequences to a previous alignment	

	print ("Adding opsin sequences to a previous alignment using Mafft")
	alignedOpsins = mafft(opsinsFile)
	
	alignedOpsinsFile = outputDir + "/" + fastaName + ".candidate_opsins.aln"
	with open (alignedOpsinsFile, 'w') as output:
		output.write (alignedOpsins)

# Converting the aligned opsins file from fasta to phylip format

	print ("Converting fasta to phylip format")
	opsinsPhylipFile= outputDir + "/" + fastaName + ".candidate_opsins.phy"
	opsinsPhylip = open(opsinsPhylipFile, 'w')

	opsins =open(alignedOpsinsFile,'r')
	lines = 0
	for line in opsins:
		if line.startswith('>'):
				lines += 1
			
	alignment = AlignIO.read(open(alignedOpsinsFile), "fasta")
	length = alignment.get_alignment_length()

	opsinsPhylip.write("%s\t%i"%(lines,length)+'\n')

	for record in alignment:
		opsinsPhylip.write(record.id + '\t' + str(record.seq)+ '\n')

	opsins.close()
	opsinsPhylip.close()
	
# Placing candidate opsin sequences in a pre-calculated opsin tree
	
	print ("Placing candidate opsin sequences in pre-calculated opsin tree using RaxML")
	opsintree = raxml(opsinsPhylipFile)
 
# Removing characters added by RAxML
	
	print("Removing characters added by RAxML")
	raxmlTree = outputDir + "/RAxML_labelledTree.Opsins"
	
	OutFilePath = outputDir + "/RAxML_labelledTree_Corrected.Opsins.tre"
	OutFile = open(OutFilePath, 'w')
	
	pattern = re.compile(r'\[I\d+\]')
	replacement = r''
	
	with open(raxmlTree, 'r') as tree:
		for line in tree:
			newtree = pattern.sub(replacement,line)
			OutFile.write(newtree)
	
	OutFile.close()
	
# Rooting the tree, then detaching the subtree with opsins only (no other GPCR)
 	
	print ("Rooting the tree then detaching the opsin clade")
	treeCorrected = outputDir + "/RAxML_labelledTree_Corrected.Opsins.tre"

	tree = Tree(treeCorrected, format=5)
	tree.set_outgroup(tree&"196_Hpoly_JEL142_contig01071")

	opsinAnc = tree.get_common_ancestor("Patiria_miniata_ops5_chaopsin", "Cavia_porcellus_H0W479")
	opsinClade = opsinAnc.detach()
	opsinClade.write(format=5, outfile= outputDir + "/RAxML_labelledTree_OpsinClade.tre")
	
# Listing the leaf names that start with QUERY in the detached opsin tree

	print ("Making a list with the queries that fell in the opsin clade")
	opsinClade = outputDir + "/RAxML_labelledTree_OpsinClade.tre"
	tree = Tree(opsinClade, format=5)
	
	opsinQueryListPath = outputDir + "/"+ fastaName + ".opsinQueryList.txt"
	outFile = open(opsinQueryListPath, 'w')
	
	for leaf in tree:
		if leaf.name.startswith('QUERY'):
			queries = leaf.name
			sequences = re.search(r'QUERY___(.+)', queries)
			print (sequences.group(1), file = outFile)
			
	outFile.close()

# Retrieving sequences in the query list from the transcriptome

	print ("Retrieving the queries that fell in the opsin clade from transcriptome")
	
	querylist = []
	with open(opsinQueryListPath,'r') as file:
		for line in file:
			line=line.strip('\n')
			querylist.extend(line.split(' '))
	
	transcriptome = open(newLongORFsFile, 'r')
	opsinSequencesFile = outputDir + "/" + fastaName + ".opsins.fasta"
	
	opsinSequences =[]
	for record in SeqIO.parse(transcriptome,'fasta'):
		for item in querylist:
			if item.strip('\n') == record.id:
				opsinSequences.append(record)
				
	SeqIO.write(opsinSequences,opsinSequencesFile,'fasta')
	transcriptome.close()
