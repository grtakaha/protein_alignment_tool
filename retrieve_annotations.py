import argparse
import requests
import os
import time
import pandas as pd
import json
import clustal_to_svg as cts
import api_funcs as af
from helpers import find_path, fasta_to_df

def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--infile") # one or more sequences in FASTA format
    parser.add_argument("-o", "--out_directory", default="./test_runs/test/")

    return parser.parse_args()   

def remove_query(df):
    for p in df.index.values:
        acc = df.loc[p, "Accession"]
        if type(acc) == pd.Series:
            for i in range(len(acc)):
                a = acc[i]
                if "QUERY_" in a:
                    df = df[df["Accession"] != a]       
        else:
            if "QUERY_" in acc:
                df = df.drop(df[df["Accession"] == acc].index)
    return df
    
def main():

    args = parse_args()

    infile = find_path(args.infile, action="r").replace("\\", "/")   
    print(f"Processing sequences from {infile}\n", flush=True)
    proteins = []
    if infile.endswith(".fasta"):
        infile_df = fasta_to_df(infile)
        infile_df = remove_query(infile_df)
        for p in infile_df.index.values:
            if type(infile_df.loc[p, "Accession"]) == pd.Series:
                # only use the first accession that IS NOT the query
                acc = infile_df.loc[p, "Accession"][0][1:]
                print(f"Duplicate accessions for ID: {p}.\nProceeding with {acc}.\n", flush=True)
            else:
                acc = infile_df.loc[p, "Accession"][1:]
            proteins.append((p, acc.split(" ")[0])) # cannot handle when something is labeled "QUERY", but has the same uniprot ID
    elif infile.endswith(".clustal") or infile.endswith(".clustal_num"):
        alignment = cts.read_alignment(infile)
        for whole_prot in alignment.proteins:
            proteins.append((whole_prot.split(" ")[0].split("|")[-1], whole_prot))
    else:
        print(f"Infile does not have a \".fasta\", \".clustal\", or \".clustal_num\" extension.\nAssuming FASTA file.\n", flush=True)
        infile_df = fasta_to_df(infile)
        infile_df = remove_query(infile_df)
        for p in infile_df.index.values:
            if type(infile_df.loc[p, "Accession"]) == pd.Series:
                # only use the first accession that IS NOT the query
                acc = infile_df.loc[p, "Accession"][0][1:]
                print(f"Duplicate accessions for ID: {p}.\nProceeding with {acc}.\n", flush=True)
            else:
                acc = infile_df.loc[p, "Accession"][1:]
            proteins.append((p, acc.split(" ")[0])) # cannot handle when something is labeled "QUERY", but has the same uniprot ID                
    
    out_directory = find_path(args.out_directory, action="w").replace("\\", "/")
    print(f"Storing outputs in {out_directory}\n", flush=True)
    
    features_df_list = []
    for p in proteins:
        prot = p[0]
        whole_prot = p[1]
        metadata = af.get_metadata(prot)

        features = metadata.get("features")
        if features is not None:
            features_df = pd.json_normalize(features)
        else:
            print(f"Could not retrieve features for {prot}...continuing with empty DataFrame\n", flush=True)
            features_df = pd.DataFrame(columns=["type","description","evidences","location.start.value","location.start.modifier","location.end.value","location.end.modifier"])

        features_df.to_csv(f"{out_directory}/{prot}.ann", sep="\t")
        
        features_df.insert(0, "prot", prot)
        features_df.insert(1, "whole_prot", whole_prot)
        features_df_list.append(features_df)
        
    all_features_df = pd.concat(features_df_list, axis=0, join='outer')    
    all_features_df.to_csv(f"{out_directory}/all.ann", sep="\t")        
        
if __name__ == "__main__":
    main()