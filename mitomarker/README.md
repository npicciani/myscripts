# MitoMarker

Use the results from MitoFinder (https://github.com/RemiAllio/MitoFinder) to add a tag (MT_) to mitochondrial
genes in the gene name field of a reference GTF file.
 
Usage: mitomarker.py [-h] -o output_directory -m mitofinder_folder -gtf GTF_file -identifier identifier_type [-version]

Label the mitochondrial genes found by MitoFinder in a reference GTF file.  


|        Option       							|   Description                              										  |
| :---           						        | :---                                       										  |
|-h, --help     						        | show this help message and exit.           										  |
|-o output_directory 						    | directory for placing output files.        										  |
|-m mitofinder_folder 							| folder with mitofinder final results.      									      |
|-gtf GTF_file       						    | original GTF file with mitochondrial genes.  										  |
|-identifier identifier_type					| type of gene/transcript identifier in the transcriptome assembly (type_1 or type_2) |
|-version             | show program's version number and exit.    |

 
 
## Note on Gene/Transcript Identifier	

The reference GTF file to be tagged by this script was originally based on transcriptome assemblies built with old Trinity versions.
Therefore, the default gene/transcript identifier (type_1) uses transcript names that follow the format "compXX_cXX_seqXX".  
If the assembly was made using a newer version of Trinity, the alternative gene/transcript identifier (type_2) will follow the format "TRINITY_DNXXXXX_cX_gX_iX".
You can set any of these two types by passing the argument -identifier when calling MitoMarker.