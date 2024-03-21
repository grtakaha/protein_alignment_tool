import argparse
import requests
import os
import time
import pandas as pd
import json
import api_funcs as af
from helpers import find_path

def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--out_directory", default="./test_runs/test/")
    parser.add_argument("-s", "--stype", default="protein")
    parser.add_argument("-e", "--email", default="")
    parser.add_argument("-t", "--title", default="alignment")

    return parser.parse_args()

def main():
    
    args = parse_args()
    
    infile = find_path(args.infile, action="r").replace("\\", "/")
    print(f"Processing sequences from {infile} \n", flush=True)
    with open(infile, "r") as f:
        seqs = "".join(f.readlines())
    
    out_directory = find_path(args.out_directory, action="w").replace("\\", "/")
    print(f"Storing outputs in {out_directory}\n", flush=True)
    
    stype = args.stype
    if stype not in ["dna", "protein"]:
        print("Given stype not found. Please specify \"-stype dna\" OR \"-stype protein\n", flush=True)
        print("Defaulting to stype=\"protein\"\n", flush=True)
        stype = "protein"  
    
    title = args.title
    
    ID = af.align(args.email, stype, title, seqs)
    alignment = af.get_alignment(ID)
    pim = af.get_pim(ID)
    
    alignment_file = f"{out_directory}/{title}.clustal_num"
    pim_file = f"{out_directory}/{title}.pim"
    
    print(f"Saving {title} alignment to {alignment_file}.", flush=True)
    with open(alignment_file, "w") as f:
        f.write(alignment) # just overwriting it if it exists
    with open(pim_file, "w") as p:
        p.write(pim) # just overwriting if it exists
        
if __name__ == "__main__":
    main()