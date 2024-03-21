'''
Helper functions for reading files.

Extra imports:
    pandas
        https://pandas.pydata.org/docs/getting_started/install.html
    protein_classes
        https://github.com/grtakaha/protein_alignment_tool/blob/main/protein_classes.py

Functions:
    find_path(str, str) -> str
    fasta_to_df(str) -> pandas.DataFrame
    read_alignment(str, int) -> protein_classes.Alignment
'''
import os
import sys
import pandas as pd
import protein_classes as pc

# change action to rf, wf, rd, wd instead for read vs. write, file vs. directory
# change text output accordingly.
def find_path(path, action):
    '''
    Returns full path of inputs based on current working directory.

            Parameters:
                    path (str): Full or relative path based on current working directory.
                                "./", "../", etc. are acceptable.
                    action (str): Either "r" (for input files) or "w" (for output directories).

            Returns:
                    abspath (str): Full path of the given input.
    '''

    if path[-1] == "/":
        abspath = os.path.abspath(path) + "/"
    else:
        abspath = os.path.abspath(path)

    if action == "r":
        if not os.path.isfile(abspath):
            print(f"No file found at {abspath}. Exiting.\n", flush=True)
            sys.exit() # check if working
    elif action == "w":
        parents = "/".join(abspath.split("/")[:-1])
        if not os.path.isdir(parents):
            print(f"Output directory does not exist. Creating parent(s): {parents}\n", flush=True)
            os.makedirs(parents)
    return abspath

def fasta_to_df(file):
    '''
    Takes in a FASTA-formatted file and returns a pandas DataFrame of the input sequences.

            Parameters:
                    file (str): Full path to a FASTA-formatted file

            Returns:
                    fasta_df (pandas.DataFrame): pandas DataFrame with short form UniProt ID indexes
                                                 columns=[Accession, IDs, Sequence]
    '''

    with open(file, "r", encoding="utf-8") as fasta:
        accessions = []
        ids = []
        seqs = []
        for line in fasta:
            if line[0] == ">":
                accessions.append(line.strip())
                # get rid of ">" and remove "|" if they exist (for uniprot entries)
                ids.append(line.strip().split(" ")[0][1:].split("|")[-1])
                seqs.append("")
            else:
                seqs[-1] += line.strip()
    return pd.DataFrame({"Accession": accessions, "IDs": ids, "Sequence": seqs}).set_index("IDs")

def read_alignment(infile, max_header=16):
    '''
    Takes in a .clustal or .clustal_num file and returns a protein_classes Alignment object.

            Parameters:
                    infile (str): Full path to a CLUSTAL- or CLUSTAL_NUM-formatted file.
                    max_header (int): Maximum length of displayed protein names
                                      when an Alignment object is written to an SVG.

            Returns:
                    alignment (protein_classes.Alignment): Alignment object created from the
                                                           given Clustal Omega alignment.
    '''

    # max_header is the max length (before 2 spaces are added) of the protein name + spaces.
    # example: max_header=10; protein name + spaces = 12
    # If the given protein name exceeds this value, it will be truncated.

    alignment = pc.Alignment() # empty alignment
    protein_key = {} # does not contain protein objects; just names and aligned seqs
    codes = ""
    header_len = 1000 # set so that it just uses an empty line if there is no header_len set
    with open(infile, "r", encoding="utf-8") as clust:
        clust.readline() # includes "\n" at the end
        for line in clust:
            if line == "\n":
                continue # skip newlines

            line = line.rstrip("\n") # remove newlines
            if line[0] != " ": # find lines with protein seqs
                l_spl = line.split(" ")
                name = l_spl[0]
                # used for codes; depends on codes coming after a block
                header_len = len(name) + line.count(" ")
                if name not in protein_key:
                    protein_key[name] = ""
                protein_key[name] += l_spl[-1].split("\t")[0] # sequence
            else:
                codes += line[header_len:] # consider adding own code checker

    # populate alignment
    alignment.add_codes(codes)
    # will add 2 spaces regardless; max_header is max accession length
    alignment.set_max_header(max_header)
    for name in protein_key:
        alignment.add_protein(pc.Protein(name, protein_key[name]))

    alignment.set_conserved_res() # set conserved residues for the completed alignment
    return alignment
