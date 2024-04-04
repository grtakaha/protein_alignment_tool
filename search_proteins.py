'''
BLASTs FASTA sequences against UniProt databases.

Extra imports:
    api_funcs
        https://github.com/grtakaha/protein_alignment_tool/blob/main/api_funcs.py
    helpers
        https://github.com/grtakaha/protein_alignment_tool/blob/main/helpers.py

Functions:
    parse_args() -> Namespace
    main()

Command-Line Arguments:
    --infile
    --out_directory
    --stype
    --email
    --num_res
'''

import argparse
import api_funcs as af
from helpers import find_path, fasta_to_df

def parse_args():
    '''
    Takes in command-line arguments and returns an argparse Namespace object.

        Returns:
            arguments (Namespace): Namespace with command-line arguments.
    '''

    parser = argparse.ArgumentParser(prog="UniProt BLAST script",
                                     description="BLASTs FASTA sequences against UniProt databases")

    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--out_directory", default="./")
    parser.add_argument("-s", "--stype", default="protein")
    parser.add_argument("-e", "--email", default="")
    parser.add_argument("-nr", "--num_res", default="10")

    return parser.parse_args()

def main():
    '''
    Parses an input FASTA file and saves UniProt BLAST results in separate directories.

        Outputs:
            One directory for each sequence (query), each with the following files:
                Table ([QUERY].tsv) and readable ([QUERY].out) BLAST results for that sequence.
                Individual FASTA files with UniProt sequences for each BLAST hit
                One FASTA file containing all protein sequences, including the query sequence.
    '''

    # TODO: Add in translation feature later maybe...or just remove dna
    args = parse_args()

    infile = find_path(args.infile, action="r").replace("\\", "/")
    print(f"Processing sequences from {infile}\n", flush=True)
    infile_df = fasta_to_df(infile)

    out_directory = find_path(args.out_directory, action="w").replace("\\", "/")
    print(f"Storing outputs in {out_directory}\n", flush=True)

    # TODO: BLAST all at once in a batch query instead of one at a time.
    blast_dict = {}
    for protein in infile_df.index.values:
        sequence = infile_df.loc[protein]["Sequence"]
        accession = infile_df.loc[protein]["Accession"] # includes ">"

        blast_id = af.blast(args.email, args.stype, protein, sequence, num_res=args.num_res)
        blast_dict[protein] = blast_id

        prot_directory = find_path(f"{out_directory}/{protein}/", action="w")

        # Overwrites output all.fasta if it exists.
        with open(f"{prot_directory}/all.fasta", "w", encoding="utf-8") as all_fasta:
            all_fasta.write(f"{accession[0]}QUERY_{accession[1:]}\n{sequence}\n")

    # TODO: BLAST all at once in a batch query instead of one at a time.
    for protein in blast_dict:
        blast_results = af.get_blast_results(blast_dict[protein])

        # prot_directory should already be created; no find_path necessary
        prot_directory = f"{out_directory}/{protein}/"

        with open(f"{prot_directory}/{protein}.out", "w", encoding="utf-8") as out:
            out.write(blast_results[0])
        with open(f"{prot_directory}/{protein}.tsv", "w", encoding="utf-8") as tsv:
            tsv.write(blast_results[1])

        # Parse blast_results[1] (tsv form)
        # Skip first line (headers) and last line (empty)
        for line in blast_results[1].split("\n")[1:-1]:
            hit = line.split("\t")[2]
            print(f"Found BLAST hit: {hit}\n", flush=True)
            hit_fasta = af.get_fasta(hit)
            hit_filename = hit_fasta.split("\n")[0].split(" ")[0].split("|")[-1]
            with open(f"{prot_directory}/{hit_filename}.fasta", "w", encoding="utf-8") as fasta:
                fasta.write(hit_fasta)
            with open(f"{prot_directory}/all.fasta", "a", encoding="utf-8") as all_fasta:
                all_fasta.write(hit_fasta)

        #with open(infile, "r", encoding="utf-8") as infile_temp:
            #all_seqs = infile_temp.read()
        #print(all_seqs)
        #blast_id = af.blast(args.email, args.stype, all_seqs, num_res=args.num_res)
        ## TODO: BLAST all at once in a batch query instead of one at a time.
        #blast_dict = {}
        #for protein in infile_df.index.values:
            #sequence = infile_df.loc[protein]["Sequence"]
            #accession = infile_df.loc[protein]["Accession"] # includes ">"
    
            #blast_dict[protein] = blast_id
    
            #prot_directory = find_path(f"{out_directory}/{protein}/", action="w")
    
            ## Overwrites output all.fasta if it exists.
            #with open(f"{prot_directory}/all.fasta", "w", encoding="utf-8") as all_fasta:
                #all_fasta.write(f"{accession[0]}QUERY_{accession[1:]}\n{sequence}\n")
        #print(blast_id)
        #blast_results = af.get_blast_results(blast_id)
        #print(blast_results[0])
        #print(blast_results[1])
                
if __name__ == "__main__":
    main()
