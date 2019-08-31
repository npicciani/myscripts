# MitoMarker

Use the results from MitoFinder (https://github.com/RemiAllio/MitoFinder) to add a tag (MT_) to mitochondrial
genes in the gene name field of a reference GTF file.
 
Usage: mitomarker.py [-h] -o output_directory -m mitofinder_folder -gtf GTF_file [-version]

Label the mitochondrial genes found by MitoFinder in a reference GTF file  

|        Option       |   Description    |
| :------------------ | -:- |
|-h, --help           | show this help message and exit.|
|-o output_directory  | directory for placing output files.|
|-m mitofinder_folder | folder with mitofinder final results.|
|-gtf GTF_file        | original GTF file with mitochondrial genes.|
|-version             | show program's version number and exit.|

 
 
## Notes

The reference GTF file to be tagged by this script was originally based on Trinity transcriptome assemblies of nonmodel species.
Therefore, transcript names follow the format "compXX_cXX_seqXX". If your transcript/gene/protein name follows a different format,
you can edit lines 69-71 to search for a different expression pattern and set the search groups that correspond to geneID and
transcriptID.