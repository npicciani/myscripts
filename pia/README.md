# PIA (Phylogenetically Informed Annotation)

This is a python version of the PIA pipeline published by Speiser et al. (2014) modified to retrieve cnidarian opsins from a
transcriptome assembly. 
 
Usage: pia.py fastafile outdir 
 
Arguments:

fastafile -- path to the transcriptome fasta file  
outdir -- path to the output directory	 
 
 
## User defined paths

First set the path to a few things in the User-defined Paths block (lines 57-65):

transdecoder_path (download from https://github.com/TransDecoder/TransDecoder)  
makeblastdb (download NCBI blast+ executables version 2.8)  
blastp (download NCBI blast+ executables version 2.8)  
mafft v7.407 (download from https://mafft.cbrc.jp/alignment/software)  
raxmlHPC-PTHREADS-SSE3 version 8.2.12 (download from https://github.com/stamatak/standard-RAxML)

This script uses the three files (baitfile, sequence alignment, and phylogenetic tree; also available at 
this repository) for retrieving candidate opsin sequences and then placing those sequences in a previously 
published phylogenetic tree (in this case, from Picciani et al. 2018). 

Once you download the baitfile, alignment and tree, change the path to those files on lines 63-65 before running the script.


## Output files

This script will return all outputs from transdecoder, makeblastdb, blastp, mafft and RaXML.   
The final opsin file will be written with the extension ".opsins.fasta"


## Running this script to mine sequences from a different gene family

This script is especifically designed to work for mining opsin genes from transcriptomes. But with a different set of baitsfile, 
alignment and phylogenetic tree, it can be used to search for genes from other families. Notice that when doing that, 
you also need to modify the placement step by setting a new outgroup and defining the boundaries for the clade in your 
phylogenetic tree that corresponds to your gene family of interest as well as the outgroup to be used for rooting the tree.

You can change those as follows:

On line 241:

	tree.set_outgroup(tree&"196_Hpoly_JEL142_contig01071")

becomes:

	tree.set_outgroup(tree&"New_outgroup_name") #include the name of the new outgroup
 
 
On line 243:

	opsinAnc = tree.get_common_ancestor("Patiria_miniata_ops5_chaopsin", "Cavia_porcellus_H0W479")

becomes:

	opsinAnc = tree.get_common_ancestor("taxa_1", "taxa_2") #include the name of two terminals whose most recent common ancestor is that of the gene family of interest
