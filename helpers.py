import pandas as pd
import os

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
