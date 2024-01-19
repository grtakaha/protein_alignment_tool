# protein_alignment_tool
<<<<<<< Updated upstream
NOTE: The following do not have command-line support for different BLAST and Clustal Omega options. They are currently set with (presumably) default parameters. Settings can be found within the scripts for those who are curious.

uniprot_fastas.py - Tool that takes one or more protein sequences (FASTA format) as input and BLASTs them against UniProt databases (RefProt and SwissProt). Saves annotations and FASTA sequences for each hit.

alignment.py - Tool that takes at least three protein sequences as input and aligns them using Clustal Omega.

clustal_to_svg.py - Tool that reformats a .clustal_num or .clustal alignment into an editable Inkscape SVG. Does not yet support annotations, but will soon.

main_tool.py - Runs one, two, or all of the above in the order listed. use "python main_tool.py -h" (or something similar via command-line) for help running this script.

Requires internet connection.
=======
Tool that takes a Clustal or Clustal_Num file as input and outputs an editable (Inkscape) SVG for use in publications. Has some annotation functionality. More to come.

2023.10.18: copied last version of files from local laptop. unsure of functionality
>>>>>>> Stashed changes
