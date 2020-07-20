"""
Published in Yang, Y., and Smith, S.A. (2014). Orthology inference in nonmodel organ-
isms using transcriptomes and low-coverage genomes: improving accu-
racy and matrix occupancy for phylogenomics. Mol. Biol. Evol. 31, 3081–3092.


Read the concatenated fasta file either with or without ends cut
write individual fasta files for each cluster
"""

import sys,os
from Bio import SeqIO

def mcl_to_fasta(all_fasta,mcl_outfile,minimal_taxa,outdir):
	print ("Reading mcl output file")
	clusterDICT = {} #key is seqID, value is clusterID
	count = 0
	if outdir[-1] != "/": outdir += "/"
	with open(mcl_outfile,"rU") as infile:
		for line in infile:
			if len(line) < 3: continue #ignore empty lines
			spls = line.strip().split('\t')
			if len(set(i.split("@")[0] for i in spls)) >= minimal_taxa:
				count += 1
				clusterID = str(count)
				for seqID in spls:
					clusterDICT[seqID] = clusterID
	print (count,"clusters with at least",minimal_taxa,"taxa read")
					
	print ("Reading the fasta file")
	handle = open(all_fasta,"r")
	for record in SeqIO.parse(handle,"fasta"):
		seqid,seq = str(record.id),str(record.seq)
		try:
			clusterID = clusterDICT[seqid]
			with open(outdir+"cluster"+clusterID+".fa","a") as outfile:
				outfile.write(">"+seqid+"\n"+seq+"\n")
		except:
			pass # Those seqs that did not go in a cluster with enough taxa
			# will not be in clusterDICT
	handle.close()
	

if __name__ =="__main__":
	if len(sys.argv) != 5:
		print ("usage: write_fasta_files_from_mcl.py all_fasta mcl_outfile minimal_taxa outDIR")
		sys.exit()
	
	mcl_to_fasta(all_fasta=sys.argv[1],mcl_outfile=sys.argv[2],\
				 minimal_taxa=int(sys.argv[3]),outdir=sys.argv[4])
	
