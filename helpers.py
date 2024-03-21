import pandas as pd
import os
import protein_classes as pc

def find_path(path, action):
    ''' Returns full path based on current working directory.
    Current working directory is where this script was run, not where it was stored.
    os.path.abspath() handles ".", "..", "../..", "./path", etc
    '''
    if path[-1] == "/":
        abspath = os.path.abspath(path) + "/"
    else:
        abspath = os.path.abspath(path)

    if action == "r":
        if not os.path.isfile(abspath):
            print(f"No file found at {abspath}. Exiting.\n", flush=True)
            exit() # check if working
    elif action == "w":
        parents = "/".join(abspath.split("/")[:-1])
        if not os.path.isdir(parents):
            print(f"Output directory does not exist. Creating parent(s): {parents}\n", flush=True)
            os.makedirs(parents)
    return abspath


def fasta_to_df(file):
    with open(file, "r") as fasta:
        accessions = []
        ids = []
        seqs = []
        for line in fasta:
            if line[0] == ">":
                accessions.append(line.strip())
                ids.append(line.strip().split(" ")[0][1:].split("|")[-1]) # get rid of ">" and remove "|" if they exist (for uniprot entries)
                seqs.append("")
            else:
                seqs[-1] += line.strip()
    return pd.DataFrame({"Accession": accessions, "IDs": ids, "Sequence": seqs}).set_index("IDs")

def read_alignment(infile, max_header=16):
    ''' max_header is the max length (before 2 spaces are added) of the protein name + spaces.
    example: max_header=10; protein name + spaces = 12
    If the given protein name exceeds this value, it will be truncated.
    '''
    alignment = pc.Alignment() # empty alignment
    protein_key = {} # does not contain protein objects; just names and aligned seqs
    codes = ""
    header_len = 1000 # set so that it just uses an empty line if there is no header_len set
    with open(infile, "r") as clust:
        header = clust.readline() # includes "\n" at the end
        for line in clust:
            if line == "\n":
                continue # skip newlines
            else:
                line = line.rstrip("\n") # remove newlines
                if line[0] != " ": # find lines with protein seqs
                    l = line.split(" ")
                    name = l[0]
                    header_len = len(name) + line.count(" ") # used for codes; depends on codes coming after a block
                    if name not in protein_key:
                        protein_key[name] = ""
                    protein_key[name] += l[-1].split("\t")[0] # sequence
                else:
                    codes += line[header_len:] # consider adding own code checker
    
    # populate alignment
    alignment.add_codes(codes)
    alignment.set_max_header(max_header) # will add 2 spaces regardless; max_header is max accession length
    for name in protein_key:
        alignment.add_protein(pc.Protein(name, protein_key[name]))
    
    alignment.set_conserved_res() # set conserved residues for the completed alignment
    return alignment
