Collection of scripts used to visualize protein sequences.
 
NOTE: The following do not have command-line support for different BLAST and Clustal Omega options. They are currently set with (presumably) default parameters. Settings can be found within the scripts for those who are curious.

# blast.py

blast.py - Tool that takes one or more protein sequences (FASTA format) as input and BLASTs them against UniProt databases (RefProt and SwissProt).

NOTE: --stype dna is currently not supported in any form. May enter an infinite loop. Please do not use --stype dna until updated.

## INPUT: FASTA-formatted file with at least one sequence.

## OUTPUT: a set of directories - one for each sequence in the original input file - that contain the following:\
	&emsp;&emsp;the BLAST results for that sequence (the query) against UniProt databases in both outfmt6 ([QUERY].tsv) and readable form ([QUERY].out)\
	&emsp;&emsp;individual FASTA files with UniProt sequences for each BLAST hit\
	&emsp;&emsp;one FASTA file containing all protein sequences, including the query sequence (all.fasta)

If used in a "multi" run, downstream commands will be run on each resulting collection of outputs.

For example: blast.py will yield multiple all.fasta (one for each query), which can be sent to both retrieve_annotations.py and alignment.py\
	&emsp;&emsp;This is why "blast annotate align" is a valid input for the --order optional argument.

Example usage:\
	&emsp;&emsp;python blast.py [-h] [-i INFILE] [-o OUT_DIRECTORY] [-s STYPE] [-e EMAIL] [-nr NUM_RES]

Example usage from main_tool.py:\
	&emsp;&emsp;python main_tool.py [-i INFILE] [-o OUT_DIRECTORY] blast [-h] [-s STYPE] [-e EMAIL] [-nr NUM_RES]

optional arguments:\
  &ensp;-h, --help            &emsp;&emsp;&emsp;show this help message and exit\
  &ensp;-i INFILE, --infile INFILE\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Full path of input file.\
  &ensp;-o OUT_DIRECTORY, --out_directory OUT_DIRECTORY\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Full path of output directory. Must end with "/".\
  &ensp;-s STYPE, --stype STYPE\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Sequence type ("protein" or "dna"). Use "dna" if aligning RNA sequences too. If run using multi, use only protein sequences and "protein" in the --stype optional argument.\
  &ensp;-e EMAIL, --email EMAIL\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Personal email. Used to submit BLAST and Clustal Omega jobs.\
  &ensp;-nr NUM_RES, --num_res NUM_RES\
                        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Number of results.


# retrieve_annotations.py

retrieve_annotations.py - Tool that takes one or more UniProt protein sequences (FASTA format) as input - UniProt is important; have not added support for non-UniProt inputs; may run infinitely if it cannot find an entry - and retrieves annotations for those sequences. Please only use UniProt sequences here until this feature is updated.
Infile must be a FASTA-formatted, CLUSTAL_NUM-formatted, or CLUSTAL-formatted file with at least one protein sequence. ALL protein sequences must be from UniProt (specifically, they must have a UniProt accession at the beginning of their names). At least three sequences are necessary if used in a multi run that includes alignment.py (align).
Outputs include:
	a collection of individual annotation files (.ann), one for each sequence in the input FASTA file
	one combined annotation file that includes all annotations for this collection of sequences (all.ann)
all.ann can be used as input for clustal_to_svg.py. See ***** ANNOTATION FORMAT ***** for help on formatting annotations by hand.

Example usage: 
python retrieve_annotations.py [-h] [-i INFILE] [-o OUT_DIRECTORY]

Example usage from main_tool.py: 
python main_tool.py [-i INFILE] [-o OUT_DIRECTORY] annotate [-h]

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Full path of input file.
  -o OUT_DIRECTORY, --out_directory OUT_DIRECTORY
                        Full path of output directory. Must end with "/".
***** END OF RETRIEVE_ANNOTATIONS.PY INSTRUCTIONS *****

***** ALIGNMENT.PY INSTRUCTIONS *****
alignment.py - Tool that takes at least three protein sequences as input and aligns them using Clustal Omega.
Infile must be a FASTA-formatted file with at least three sequences.
Outputs include:
	an alignment (.clustal_num) and associated percent identity matrix (.pim) of the given FASTA file
	
Example usage: 
python alignment.py [-h] [-i INFILE] [-o OUT_DIRECTORY] [-s STYPE] [-e EMAIL] [-t TITLE]

Example usage from main_tool.py: 
python main_tool.py [-i INFILE] [-o OUT_DIRECTORY] align [-h] [-s STYPE] [-e EMAIL] [-t TITLE]

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Full path of input file.
  -o OUT_DIRECTORY, --out_directory OUT_DIRECTORY
                        Full path of output directory. Must end with "/".
  -s STYPE, --stype STYPE
                        Sequence type ("protein" or "dna"). Use "dna" if aligning RNA sequences too. If run using multi, use only protein sequences and "protein" in the --stype optional argument.
  -e EMAIL, --email EMAIL
                        Personal email. Used to submit BLAST and Clustal Omega jobs.
  -t TITLE, --title TITLE
                        Title for the output alignment (.clustal_num) and percent identity matrix (.pim). Example: alignment1 -> alignment1.clustal_num, alignment1.pim
***** END OF ALIGNMENT.PY INSTRUCTIONS *****

***** CLUSTAL_TO_SVG.PY INSTRUCTIONS *****
clustal_to_svg.py - Tool that reformats a .clustal_num or .clustal alignment into an editable Inkscape SVG. Currently annotates conserved residues (automatic, not optional) and active site residues (requires an input annotation file, optional).
Infile must be a CLUSTAL or CLUSTAL_NUM file.
Output is a sequential set of SVGs (.svg), numbered 0, 1, 2, etc., with formatted alignments and associted conserved residues and/or annotations.

Example usage from main_tool.py: 
python clustal_to_svg.py [-h] [-i INFILE] [-o OUT_DIRECTORY] [-c CODES] [-n NUMS] [-u UNIPROT_FORMAT]
                                                   [-a ANNOTATIONS]
Example usage from main_tool.py: 
python main_tool.py [-i INFILE] [-o OUT_DIRECTORY] svg [-h] [-c CODES] [-n NUMS] [-u UNIPROT_FORMAT]
                                                       [-a ANNOTATIONS]

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Full path of input file.
  -o OUT_DIRECTORY, --out_directory OUT_DIRECTORY
                        Full path of output directory. Must end with "/".
  -c CODES, --codes CODES
                        Default FALSE. If TRUE, will add Clustal Omega conservation codes to the bottom of each aligned block.
  -n NUMS, --nums NUMS
                        Default FALSE. If TRUE, will add total residue numbers to the right side of every line.
  -u UNIPROT_FORMAT, --uniprot_format UNIPROT_FORMAT
                        Default FALSE. If TRUE, will truncate accessions according to UniProt formatting. Example: sp|P00784|PAPA1_CARPA -> PAPA1_CARPA
  -a ANNOTATIONS, --annotations ANNOTATIONS
                        Full path to annotation file. Currently only supports active site annotations. Others will be ignored. If run using multi, annotations can either be provided separately, or acquired from UniProt by including "annotate" in the --order optional argument.
***** END OF CLUSTAL_TO_SVG.PY INSTRUCTIONS *****

***** MAIN_TOOL.PY INSTRUCTIONS *****
main_tool.py - Runs one, two, or all of the above in the order given
Infile is dependent on which tool(s) are being executed. Infile should match that of the first tool being executed. Example: --order blast align svg -> use the infile format for blast.py

Example usage: 
python main_tool.py [-i INFILE] [-o OUT_DIRECTORY] {blast,annotate,align,svg,multi} [-h] [-ord ORDER [ORDER ...]] 
                                                                                    [-s STYPE] [-e EMAIL]
                                                                                    [-nr NUM_RES] [-t TITLE] [-c CODES] [-n NUMS]
                                                                                    [-u UNIPROT_FORMAT] [-a ANNOTATIONS]

Centralized Tool Manager

positional arguments:
  {blast,annotate,align,svg,multi}
                        Tool to execute

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Full path of input file.
  -o OUT_DIRECTORY, --out_directory OUT_DIRECTORY
                        Full path of output directory. Must end with "/".
  -ord ORDER [ORDER ...], --order ORDER [ORDER ...]
                        Order of tools to run if "multi" is used as a positional argument. There are currently limited ways to run multi (inputs and outputs will vary depending on start and end):
                        blast annotate align svg
			blast align annotate svg
			blast annotate align
			blast align annotate
                        blast annotate
                        blast align svg
                        blast align
                        blast
                        annotate align svg
                        annotate align
                        annotate svg
                        annotate
                        align svg
                        align
                        svg
  -s STYPE, --stype STYPE
                        Sequence type ("protein" or "dna"). Use "dna" if aligning RNA sequences too. If run using multi, use only protein sequences and "protein" in the --stype optional argument.
  -e EMAIL, --email EMAIL
                        Personal email. Used to submit BLAST and Clustal Omega jobs.
  -nr NUM_RES, --num_res NUM_RES
                        Number of results.
  -t TITLE, --title TITLE
                        Title for the output alignment (.clustal_num) and percent identity matrix (.pim). Example: alignment1 -> alignment1.clustal_num, alignment1.pim
  -c CODES, --codes CODES
                        Default FALSE. If TRUE, will add Clustal Omega conservation codes to the bottom of each aligned block.
  -n NUMS, --nums NUMS
                        Default FALSE. If TRUE, will add total residue numbers to the right side of every line.
  -u UNIPROT_FORMAT, --uniprot_format UNIPROT_FORMAT
                        Default FALSE. If TRUE, will truncate accessions according to UniProt formatting. Example: sp|P00784|PAPA1_CARPA -> PAPA1_CARPA
  -a ANNOTATIONS, --annotations ANNOTATIONS
                        Full path to annotation file. Currently only supports active site annotations. Others will be ignored. If run using multi, annotations can either be provided separately, or acquired from UniProt by including "annotate" in the --order optional argument.
***** END OF MAIN_TOOL.PY INSTRUCTIONS *****

***** ANNOTATION FORMAT *****
As of 2024.03.19, only annotations of type "Active site" will be used. This will be updated in the future.
Annotation files that include the following columns and VALUES (tab-delimited) can be used as inputs for clustal_to_svg.py:
	prot	whole_prot	type	location.start.value	location.end.value
ARBITRARY_INDEX	UNIPROT_FORMAT_ACC	FULL_ACCESSION	ANNOTATION_TYPE	START	END

Example of a valid annotation file:

	prot	whole_prot	type	description	evidences	location.start.value	location.start.modifier	location.end.value	location.end.modifier	featureId	featureCrossReferences	ligand.name	ligand.id	ligand.note	alternativeSequence.originalSequence	alternativeSequence.alternativeSequences
0	A0A164XSU3_DAUCS	tr|A0A164XSU3|A0A164XSU3_DAUCS	Signal		[{'evidenceCode': 'ECO:0000256', 'source': 'SAM', 'id': 'SignalP'}]	1	EXACT	27	EXACT							
1	A0A164XSU3_DAUCS	tr|A0A164XSU3|A0A164XSU3_DAUCS	Chain	Cysteine protease	[{'evidenceCode': 'ECO:0000256', 'source': 'SAM', 'id': 'SignalP'}]	28	EXACT	354	EXACT	PRO_5018759784						
2	A0A164XSU3_DAUCS	tr|A0A164XSU3|A0A164XSU3_DAUCS	Domain	Cathepsin propeptide inhibitor	[{'evidenceCode': 'ECO:0000259', 'source': 'SMART', 'id': 'SM00848'}]	50	EXACT	106	EXACT							
3	A0A164XSU3_DAUCS	tr|A0A164XSU3|A0A164XSU3_DAUCS	Domain	Peptidase C1A papain C-terminal	[{'evidenceCode': 'ECO:0000259', 'source': 'SMART', 'id': 'SM00645'}]	136	EXACT	351	EXACT							
0	PAPA2_CARPA	sp|P14080|PAPA2_CARPA	Signal		[{'evidenceCode': 'ECO:0000255'}]	1	EXACT	18	EXACT							
1	PAPA2_CARPA	sp|P14080|PAPA2_CARPA	Propeptide	Activation peptide	[{'evidenceCode': 'ECO:0000269', 'source': 'PubMed', 'id': '2106878'}, {'evidenceCode': 'ECO:0000269', 'source': 'PubMed', 'id': '2500950'}]	19	EXACT	134	EXACT	PRO_0000026408						
2	PAPA2_CARPA	sp|P14080|PAPA2_CARPA	Chain	Chymopapain		135	EXACT	352	EXACT	PRO_0000026409						
3	PAPA2_CARPA	sp|P14080|PAPA2_CARPA	Active site		[{'evidenceCode': 'ECO:0000255', 'source': 'PROSITE-ProRule', 'id': 'PRU10088'}]	159	EXACT	159	EXACT							
4	PAPA2_CARPA	sp|P14080|PAPA2_CARPA	Active site		[{'evidenceCode': 'ECO:0000255', 'source': 'PROSITE-ProRule', 'id': 'PRU10089'}]	293	EXACT	293	EXACT							
5	PAPA2_CARPA	sp|P14080|PAPA2_CARPA	Active site		[{'evidenceCode': 'ECO:0000255', 'source': 'PROSITE-ProRule', 'id': 'PRU10090'}]	313	EXACT	313	EXACT							
6	PAPA2_CARPA	sp|P14080|PAPA2_CARPA	Glycosylation	N-linked (GlcNAc...) asparagine	[{'evidenceCode': 'ECO:0000255', 'source': 'PROSITE-ProRule', 'id': 'PRU00498'}]	86	EXACT	86	EXACT							
							
***** END OF ANNOTATION FORMAT *****

Requires internet connection.
