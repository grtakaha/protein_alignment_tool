Collection of scripts used to visualize protein sequences.

NOTE: The following do not have command-line support for different BLAST and Clustal Omega options. They are currently set with (presumably) default parameters. Settings can be found within the scripts for those who are curious.

NOTE: As per NCBI guidelines, please do not use search_proteins.py for more than 100 searches in a 24 hour period. Limits are in place in an attempt to comply with guidelines listed here:\
	&emsp;&emsp;https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html#rest\
Currently working on local alternative that will use the NCBI BLAST+ CLI instead.

NOTE: Please follow usage guidelines found here for both search_proteins.py and alignment.py:
	&emsp;&emsp;https://ebi-biows.gitdocs.ebi.ac.uk/documentation/webservices/

NOTE: These scripts make use of EMBL-EBI and/or NCBI resources. References for tools and databases used here include:\

EMBL-EBI's Job Dispatcher:\
Madeira, Fábio, Matt Pearce, Adrian R. N. Tivey, Prasad Basutkar, Joon Lee, Ossama Edbali, Nandana Madhusoodanan, Anton Kolesnikov, and Rodrigo Lopez. \
“Search and Sequence Analysis Tools Services from EMBL-EBI in 2022.” \
Nucleic Acids Research 50, no. W1 (July 5, 2022): W276–79. https://doi.org/10.1093/nar/gkac240.

UniProt:\
The UniProt Consortium. \
“UniProt: The Universal Protein Knowledgebase in 2023.” \
Nucleic Acids Research 51, no. D1 (January 6, 2023): D523–31. https://doi.org/10.1093/nar/gkac1052.

NCBI:\
Sayers, Eric W, Evan E Bolton, J Rodney Brister, Kathi Canese, Jessica Chan, Donald C Comeau, Ryan Connor, et al. \
“Database Resources of the National Center for Biotechnology Information.” \
Nucleic Acids Research 50, no. D1 (December 1, 2021): D20–26. https://doi.org/10.1093/nar/gkab1112.

Clustal Omega:\
Sievers, Fabian, Andreas Wilm, David Dineen, Toby J Gibson, Kevin Karplus, Weizhong Li, Rodrigo Lopez, et al. \
“Fast, Scalable Generation of High‐quality Protein Multiple Sequence Alignments Using Clustal Omega.” \
Molecular Systems Biology 7, no. 1 (January 2011): 539. https://doi.org/10.1038/msb.2011.75.

Sievers, Fabian, and Desmond G. Higgins. \
“Clustal Omega for Making Accurate Alignments of Many Protein Sequences.” \
Protein Science: A Publication of the Protein Society 27, no. 1 (January 2018): 135–45. https://doi.org/10.1002/pro.3290.

BLAST+:\
Camacho, Christiam, George Coulouris, Vahram Avagyan, Ning Ma, Jason Papadopoulos, Kevin Bealer, and Thomas L. Madden. \
“BLAST+: Architecture and Applications.” \
BMC Bioinformatics 10, no. 1 (December 2009): 1–9. https://doi.org/10.1186/1471-2105-10-421.

https://blast.ncbi.nlm.nih.gov/doc/blast-help/references.html#references

# REQUIREMENTS
General:\
	&emsp;&emsp;Internet connection (when running blast.py, retrieve_annotations.py, and/or alignment.py)\
	&emsp;&emsp;Python 3.7+

Python Libraries:\
	&emsp;&emsp;pandas\
	&emsp;&emsp;requests

# search_proteins.py

Takes one or more protein sequences (FASTA format) as input and BLASTs them against UniProt databases (uniprotkb_refprotswissprot).

NOTE: --stype dna is currently not supported in any form. May enter an infinite loop. Please do not use --stype dna until updated.

INPUT: FASTA-formatted file with at least one sequence

OUTPUT: A set of directories - one for each sequence in the original input file - that contain the following:\
	&emsp;&emsp;the BLAST results for that sequence (the query) against UniProt databases in both outfmt6 ([QUERY].tsv) and readable form ([QUERY].out)\
	&emsp;&emsp;individual FASTA files with UniProt sequences for each BLAST hit\
	&emsp;&emsp;one FASTA file containing all protein sequences, including the query sequence (all.fasta)

If used in a "multi" run, downstream commands will be run on each resulting collection of outputs.

For example: blast.py will yield multiple all.fasta (one for each query), which can be sent to both retrieve_annotations.py and alignment.py. This is why "blast annotate align" is a valid input for the --order optional argument.

Example usage from main_tool.py:\
  &ensp;python main_tool.py [-i INFILE] [-o OUT_DIRECTORY] blast [-h] [-s STYPE] [-e EMAIL] [-nr NUM_RES]

Example usage:\
  &ensp;python blast.py [-h] [-i INFILE] [-o OUT_DIRECTORY] [-s STYPE] [-e EMAIL] [-nr NUM_RES]

optional arguments:\
  &ensp;-h, --help&emsp;&emsp;&emsp;show this help message and exit\
  &ensp;-i INFILE, --infile INFILE\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Full path of input file.\
  &ensp;-o OUT_DIRECTORY, --out_directory OUT_DIRECTORY\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Full path of output directory. Must end with "/".\
  &ensp;-s STYPE, --stype STYPE\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Sequence type ("protein" or "dna"). Use "dna" if aligning RNA sequences too. If run using multi, use only protein sequences and "protein" in the --stype optional argument.\
  &ensp;-e EMAIL, --email EMAIL\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Personal email. Used to submit BLAST and Clustal Omega jobs.\
  &ensp;-nr NUM_RES, --num_res NUM_RES\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Number of results.


# retrieve_annotations.py

Takes one or more UniProt protein sequences (FASTA format) as input and retrieves annotations for those sequences. 

NOTE: Please only use UniProt sequences here until this feature is updated. May run infinitely if it cannot find a given entry.

INPUT: A FASTA-formatted, CLUSTAL_NUM-formatted, or CLUSTAL-formatted file with at least one protein sequence. ALL protein sequences must be from UniProt (specifically, they must have a UniProt accession at the beginning of their names). At least three sequences are necessary if used in a multi run that includes alignment.py (align).

OUTPUT: A collection of files that includes the following:\
	&emsp;&emsp;individual annotation files (.ann), one for each unique sequence in the input file\
	&emsp;&emsp;one combined annotation file that includes all annotations for this collection of sequences (all.ann)

NOTE: all.ann can be used as input for clustal_to_svg.py. See **ANNOTATION FORMAT** for help formatting annotations by hand.

Example usage from main_tool.py:\
  &ensp;python main_tool.py [-i INFILE] [-o OUT_DIRECTORY] annotate [-h]

Example usage:\
  &ensp;python retrieve_annotations.py [-h] [-i INFILE] [-o OUT_DIRECTORY]

optional arguments:\
  &ensp;-h, --help &emsp;&emsp;&emsp;show this help message and exit\
  &ensp;-i INFILE, --infile INFILE\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Full path of input file.\
  &ensp;-o OUT_DIRECTORY, --out_directory OUT_DIRECTORY\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Full path of output directory. Must end with "/".


# alignment.py

Takes at least three protein sequences as input and aligns them using Clustal Omega.

INPUT: A FASTA-formatted file with at least three sequences

OUTPUT: An alignment (.clustal_num) and associated percent identity matrix (.pim) of the given FASTA file

Example usage from main_tool.py:\
  &ensp;python main_tool.py [-i INFILE] [-o OUT_DIRECTORY] align [-h] [-s STYPE] [-e EMAIL] [-t TITLE]

Example usage:\
  &ensp;python alignment.py [-h] [-i INFILE] [-o OUT_DIRECTORY] [-s STYPE] [-e EMAIL] [-t TITLE]

optional arguments:\
  &ensp;-h, --help &emsp;&emsp;&emsp;show this help message and exit\
  &ensp;-i INFILE, --infile INFILE\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Full path of input file.\
  &ensp;-o OUT_DIRECTORY, --out_directory OUT_DIRECTORY\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Full path of output directory. Must end with "/".\
  &ensp;-s STYPE, --stype STYPE\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Sequence type ("protein" or "dna"). Use "dna" if aligning RNA sequences too. If run using multi, use only protein sequences and "protein" in the --stype optional argument.\
  &ensp;-e EMAIL, --email EMAIL\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Personal email. Used to submit BLAST and Clustal Omega jobs.\
  &ensp;-t TITLE, --title TITLE\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Title for the output alignment (.clustal_num) and percent identity matrix (.pim). Example: alignment1 -> alignment1.clustal_num, alignment1.pim


# clustal_to_svg.py

Reformats a .clustal_num or .clustal alignment into an editable Inkscape SVG. Currently annotates conserved residues (automatic, not optional) and active site residues (requires an input annotation file, optional).

INPUT: A CLUSTAL or CLUSTAL_NUM file

OUTPUT: A sequential set of SVGs (.svg), numbered 0, 1, 2, etc., with formatted alignments and associated conserved residues and/or annotations.

Example usage from main_tool.py:\
  &ensp;python main_tool.py [-i INFILE] [-o OUT_DIRECTORY] svg [-h] [-c CODES] [-n NUMS] [-u UNIPROT_FORMAT]
                                                         &nbsp;&ensp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;[-a ANNOTATIONS]

Example usage:\
  &ensp;python clustal_to_svg.py [-h] [-i INFILE] [-o OUT_DIRECTORY] [-c CODES] [-n NUMS] [-u UNIPROT_FORMAT]
                           &ensp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;[-a ANNOTATIONS]

optional arguments:\
  &ensp;-h, --help &emsp;&emsp;&emsp;show this help message and exit\
  &ensp;-i INFILE, --infile INFILE\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Full path of input file.\
  &ensp;-o OUT_DIRECTORY, --out_directory OUT_DIRECTORY\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Full path of output directory. Must end with "/".\
  &ensp;-c CODES, --codes CODES\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Default FALSE. If TRUE, will add Clustal Omega conservation codes to the bottom of each aligned block.\
  &ensp;-n NUMS, --nums NUMS\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Default FALSE. If TRUE, will add total residue numbers to the right side of every line.\
  &ensp;-u UNIPROT_FORMAT, --uniprot_format UNIPROT_FORMAT\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Default FALSE. If TRUE, will truncate accessions according to UniProt formatting. Example: sp|P00784|PAPA1_CARPA -> PAPA1_CARPA\
  &ensp;-a ANNOTATIONS, --annotations ANNOTATIONS\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Full path to annotation file. Currently only supports active site annotations. Others will be ignored. If run using multi, annotations can either be provided separately, or acquired from UniProt by including "annotate" in the --order optional argument.


# main_tool.py

Runs one, two, or all of the above in the order given.

INPUT: Depends on which tool(s) are being executed. Should be an acceptable input of the first tool being executed. Inputs for runs that start with annotate.py may additionally be limited by what can be passed to downstream tools.  Example: multi --order annotate svg -> must use .clustal or .clustal_num file as input (.fasta file cannot be used by clustal_to_svg.py)

OUTPUT: Depends on which tool(s) are being executed. Each tool will have its own output if it is included in a run. All outputs will be split into separate directories in a multi run that includes blast.py.

Example usage:\
  &ensp;python main_tool.py [-i INFILE] [-o OUT_DIRECTORY] {blast,annotate,align,svg,multi} [-h] [-ord ORDER [ORDER...]]\
 &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;[-s STYPE] [-e EMAIL]\
 &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;[-nr NUM_RES] [-t TITLE]\
 &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;[-c CODES] [-n NUMS]\
 &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;[-u UNIPROT_FORMAT]\
 &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;[-a ANNOTATIONS]

Centralized Tool Manager

positional arguments:\
  &ensp;{blast,annotate,align,svg,multi}\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Tool to execute

optional arguments:\
  &ensp;-h, --help &emsp;&emsp;&emsp;show this help message and exit\
  &ensp;-i INFILE, --infile INFILE\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Full path of input file.\
  &ensp;-o OUT_DIRECTORY, --out_directory OUT_DIRECTORY\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Full path of output directory. Must end with "/".\
  &ensp;-ord ORDER [ORDER ...], --order ORDER [ORDER ...]\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Order of tools to run if "multi" is used as a positional argument. There are currently limited ways to run multi (inputs and outputs will vary depending on start and end):\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;blast annotate align svg\
			&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;blast annotate align\
			&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;blast annotate\
			&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;blast align annotate svg\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;blast align annotate\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;blast align svg\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;blast align\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;blast\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;annotate align svg\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;annotate align\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;annotate svg\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;annotate\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;align annotate svg\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;align annotate\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;align svg\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;align\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;svg\
  &ensp;-s STYPE, --stype STYPE\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Sequence type ("protein" or "dna"). Use "dna" if aligning RNA sequences too. If run using multi, use only protein sequences and "protein" in the --stype optional argument.\
  &ensp;-e EMAIL, --email EMAIL\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Personal email. Used to submit BLAST and Clustal Omega jobs.\
  &ensp;-nr NUM_RES, --num_res NUM_RES\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Number of results.\
  &ensp;-t TITLE, --title TITLE\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Title for the output alignment (.clustal_num) and percent identity matrix (.pim). Example: alignment1 -> alignment1.clustal_num, alignment1.pim\
  &ensp;-c CODES, --codes CODES\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Default FALSE. If TRUE, will add Clustal Omega conservation codes to the bottom of each aligned block.\
  &ensp;-n NUMS, --nums NUMS\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Default FALSE. If TRUE, will add total residue numbers to the right side of every line.\
  &ensp;-u UNIPROT_FORMAT, --uniprot_format UNIPROT_FORMAT\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Default FALSE. If TRUE, will truncate accessions according to UniProt formatting. Example: sp|P00784|PAPA1_CARPA -> PAPA1_CARPA\
  &ensp;-a ANNOTATIONS, --annotations ANNOTATIONS\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Full path to annotation file. Currently only supports active site annotations. Others will be ignored. If run using multi, annotations can either be provided separately, or acquired from UniProt by including "annotate" in the --order optional argument.


# ANNOTATION FORMAT

As of 2024.03.19, only annotations of type "Active site" will be used. This will be updated in the future.

Annotation files that include the following columns and VALUES (tab-delimited) can be used as inputs for clustal_to_svg.py:

		prot	whole_prot	type	location.start.value	location.end.value
	ARBITRARY_INDEX	UNIPROT_FORMAT_ACC	FULL_ACCESSION	ANNOTATION_TYPE	START	END

A real example might look like:

		prot	whole_prot	type	location.start.value	location.end.value
	0	PAPA1_CARPA	sp|P00784|PAPA1_CARPA	Active site	158	158

A truncated, but real example of a valid annotation file can be found in annotation_example.ann.
