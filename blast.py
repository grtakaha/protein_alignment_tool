import argparse
import requests
import os
import time
import pandas as pd
import json
import api_funcs as af
from helpers import find_path, fasta_to_df

def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--out_directory", default="./")
    parser.add_argument("-s", "--stype", default="protein") # program is tied to stype (dna:blastx, rna:blastx, protein:blastp)
    parser.add_argument("-e", "--email", default="")
    parser.add_argument("-nr", "--num_res", default="10")
    
    return parser.parse_args()
        
def main():

    args = parse_args()

    infile = find_path(args.infile, action="r").replace("\\", "/")
    print(f"Processing sequences from {infile}\n", flush=True)
    infile_df = fasta_to_df(infile)
    
    out_directory = find_path(args.out_directory, action="w").replace("\\", "/")
    print(f"Storing outputs in {out_directory}\n", flush=True)
    
    # add in translation feature later maybe...or just remove dna
    stype = args.stype
    if stype == "dna":
        program = "blastx"
    elif stype == "protein":
        program = "blastp"
    else:
        print("Given stype not found. Please specify \"-stype dna\" for blastx OR \"-stype protein\" for blastp", flush=True)
        print("Defaulting to stype=\"protein\" and program=\"blastx\"\n", flush=True)
        stype = "protein"
        program = "blastx"    
    
    blast_dict = {}
    for protein in infile_df.index.values:
        sequence = infile_df.loc[protein]["Sequence"]
        accession = infile_df.loc[protein]["Accession"] # includes ">"
        
        blast_id = af.blast(args.email, program, stype, protein, sequence, num_res=args.num_res)
        blast_dict[protein] = blast_id
        
        prot_directory = find_path(f"{out_directory}/{protein}/", action="w")     
        with open(f"{prot_directory}/all.fasta", "w") as all_fasta:
            all_fasta.write(f"{accession[0]}QUERY_{accession[1:]}\n{sequence}\n") # just overwriting it if it exists
            
    for protein in blast_dict:
        # run blast_ids all at once and then retrieve all at once
        blast_results = af.get_blast_results(blast_dict[protein]) # consider retrieving ids all at once, then getting results one by one after running
        
        # prot_directory should already be created; no find_path necessary
        prot_directory = f"{out_directory}/{protein}/"
        
        with open(f"{prot_directory}/{protein}.out", "w") as out:
            out.write(blast_results[0])
        with open(f"{prot_directory}/{protein}.tsv", "w") as tsv:
            tsv.write(blast_results[1])         
        
        features_df_list = [] # unused?
        # parse blast_results[1] (tsv form)
        for line in blast_results[1].split("\n")[1:-1]: # skip first line (headers) and last line (empty); don't bother reading into a dataframe
            hit = line.split("\t")[2]
            print(f"Found BLAST hit: {hit}\n", flush=True)
            hit_fasta = af.get_fasta(hit)
            hit_filename = hit_fasta.split("\n")[0].split(" ")[0].split("|")[-1]
            with open(f"{prot_directory}/{hit_filename}.fasta", "w") as fasta:
                fasta.write(hit_fasta)
            with open(f"{prot_directory}/all.fasta", "a") as all_fasta:
                all_fasta.write(hit_fasta)
                
if __name__ == "__main__":
    main()