#!/usr/bin/env python
# Python 3.7
# -*- coding: utf-8 -*-

'''
created on June 3, 2019
by Natasha Picciani 
usage: dataMining_v3.py taxaList taxaRestriction outputDir

This script searches the NCBI Nucleotide database, pull sequences for a list of taxa, 
and cluster them using MCL

Changes in version 2:
- added an argument to restrict the search to a certain group of NCBI taxa (e.g. Metazoa, Arthropoda)

Changes in version 3:
- Depeng added a sleep time to prevent HTTP errors

'''

import re
import subprocess
from pathlib import PurePosixPath
from Bio import Entrez
from Bio import SeqIO
import blast_to_mcl_py3_ed as bm
import write_fasta_files_from_mcl_py3 as wt
import os,errno
import sys
from urllib.error import HTTPError

def mkdir_p(path): #created by Andrew Swafford
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise
    return(path)


def nameFix(seqNew,seqOld,outputPath):
	seqShort = SeqIO.to_dict(SeqIO.parse(seqNew,'fasta'))
	seqOriginal = SeqIO.to_dict(SeqIO.parse(seqOld,'fasta')) #no white spaces in the headers

	with open(seqOld,'r') as mask:
		seqO = mask.read().splitlines()
	list = []
	for name in seqShort.keys():
		for k in seqO:
			if name in k:
				list.extend(k.split('\n'))            
	matchList = []
	for line in list:
		line = re.sub(r">(.+)",r"\1",line)
		matchList.append(line)

	with open(outputPath, "w") as output:
			for seq_id in (line.strip() for line in matchList):
				if seq_id:
					SeqIO.write(seqOriginal[seq_id],output,'fasta')

# Path to programs
mclPath = "/home/depengli/anaconda3/bin/mcl"
blastnPath = "/home/depengli/anaconda3/bin/blastn"
makeBlastDBPath = "/home/depengli/anaconda3/bin/makeblastdb"

# Email settings for NCBI and user inputs
Entrez.email = "n_picciani@hotmail.com" #user defined
taxaListPath = sys.argv[1] #list of taxa for the search
taxaRestrict= sys.argv[2] #major group to restrict the search
outDir = sys.argv[3] #output directory

# Listing the IDs of sequences available on the NCBI nucleotide database from a list of taxa
initialAcc=[]
with open (taxaListPath, 'r') as taxaList:
    for line in taxaList:
        line=line.strip('\n')
        searchTerm = line + "[organism]" + " AND " + taxaRestrict + "[organism]"
        print (searchTerm)
        handle = Entrez.esearch(db="nucleotide", retmax=100000, term=searchTerm, idtype="acc")
        record = Entrez.read(handle)
        handle.close()
        initialAcc.extend(record["IdList"])
        
# Removing any possible duplicates in the ID list
accNumbersFilePath = outDir + "/accessionNumbers.txt"
finalAcc = []
with open(accNumbersFilePath, 'w') as outfile:
    for item in initialAcc:
        if item not in finalAcc:
            finalAcc.append(item)
    for item in finalAcc:
        print("%s" %item, file=outfile)
        
print("Exported list of accession numbers")

# Pulling fasta files from list of accession numbers
fastaRecordsPath = outDir + "/fastaRecords.fasta"
fastaRecords=[]
for item in finalAcc:
#try/except block to prevent HTTPError /depengli
    try:
        handle = Entrez.efetch(db="nucleotide", id=item, rettype="fasta", retmode="text") #automatically uses an HTTP POST if there are over 200 identifiers
    except HTTPError:
        time.sleep(5)
        handle = Entrez.efetch(db="nucleotide", id=item, rettype="fasta", retmode="text")
    fastaRecords.extend(handle)
    handle.close()

with open(fastaRecordsPath, 'w') as f:
    for item in fastaRecords:
        f.write("%s" %item)

print("Exported fasta file")


# Replacing white spaces in the fasta headers
pattern = re.compile(r' ')
replacement = r'_'
fastaFileUnderscore = outDir + "/fastaRecords_underscores.fasta"
outFile = open(fastaFileUnderscore,'w')

with open(fastaRecordsPath, 'r') as file:
    for line in file:
        if line.startswith('>'):
            newline = pattern.sub(replacement,line)
            outFile.write(newline)
        else:
            outFile.write(line)
            
outFile.close()

# Shortening the sequence names
pattern = re.compile(r'(>\w+\.\d+)\s.+')
replacement = r'\1'
fastaFile = outDir + "/fastaRecords_shortheaders.fasta"
outFile = open(fastaFile,'w')

with open(fastaRecordsPath, 'r') as file:
    for line in file:
        if line.startswith('>'):
            newline = pattern.sub(replacement,line)
            outFile.write(newline)
        else:
            outFile.write(line)
            
outFile.close()

# Building a blast database 
fastaName = PurePosixPath(fastaFile).stem

print ("Building blast database for " + fastaName)

blastdbName = outDir + "/" + fastaName + ".blastdb"
blastdb_log = subprocess.run([makeBlastDBPath,"-in",fastaFile,\
                              "-dbtype","nucl","-parse_seqids","-out",blastdbName], \
                             check=True, capture_output=True)

# Blasting the database with a nucleotide query
rawblast = outDir + "/" + fastaName + ".all.rawblast"
blastn_log = subprocess.run([blastnPath,"-db", blastdbName,\
                                    "-query", fastaFile,"-evalue","10","-max_target_seqs","1000","-out",rawblast,\
                                    "-outfmt","6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore"], \
                                    check=True, capture_output=True)
print("Blasted all against all")

# Clustering with MCL
fracList = [0.3,0.4,0.5] #list of hit fractions
infList = [1,1.4,2] #list of inflation values

for frac in fracList:
	fracOutDir = outDir + "/hit_fraction_%s"%frac 
	mkdir_p(fracOutDir) #making a directory for hit-fraction files
	bm.blast_to_mcl(rawblast, frac,out_dir = fracOutDir) #convert blast format to MCL input format
	mclFracFile = fracOutDir + "/fastaRecords_shortheaders.all.rawblast.hit-frac%s.minusLogEvalue"%frac
	for inf in infList:
		inflation = str(inf)
		infDir = fracOutDir + "/" + inflation
		mkdir_p(infDir)
		infOutFile = infDir+"/"+os.path.basename(mclFracFile)+".e5_I%s"%inflation
		mcl_log = subprocess.run([mclPath,mclFracFile,"--abc","-te",\
                              "5","-tf","gq(5)","-I",inflation,"-o",infOutFile], \
                             check=True, capture_output=True)
		wt.mcl_to_fasta(fastaFile,infOutFile,5,infDir)
		for filename in os.listdir(infDir): #converting the short headers of fasta files back to their original names		
			if filename.endswith(".fa"):
				base = os.path.basename(filename).split(".")
				outputPath = os.path.normpath(os.path.join(infDir,"repaired_" + filename))
				filenamePath = infDir + "/" + filename
				nameFix(filenamePath,fastaFileUnderscore,outputPath)


print("Clustered sequences using MCL with hit fractions %.1f, %.1f and %.1f as well as inflation values %.1f, %.1f and %.1f" % (fracList[0],fracList[1],fracList[2],infList[0],infList[1],infList[2]))
