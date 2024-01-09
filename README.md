# protein_alignment_tool
NOTE: The following do not have command-line support for different BLAST and Clustal Omega options. They are 
\n
\nuniprot_fastas.py - Tool that takes one or more protein sequences (FASTA format) as input and BLASTs them against UniProt databases (RefProt and SwissProt). Saves annotations and FASTA sequences for each hit.
\nalignment.py - Tool that takes at least three protein sequences as input and aligns them using Clustal Omega.
\nclustal_to_svg.py - Tool that reformats a .clustal_num or .clustal alignment into an editable Inkscape SVG. Does not yet support annotations, but will soon.
\nmain_tool.py - Runs one, two, or all of the above in the order listed. use "python main_tool.py -h" (or something similar via command-line) for help running this script.
\n
\nRequires internet connection.
