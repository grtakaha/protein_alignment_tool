import argparse
import requests
import os
import time
import pandas as pd
import json
import clustal_to_svg as cts

class api_tools():
    
    def get_metadata(ID):
        # for now, just returns entire json...consider extracting annotations only
        url_metadata = f"https://rest.uniprot.org/uniprotkb/{ID}"
        
        # rerun the request until it returns 200
        # it's possible this becomes an endless loop...deal with that later
        # Set the headers to accept JSON
        headers = {"Accept": "application/json"}      
        while True:
            try:
                response_metadata = requests.get(url_metadata, headers)
                if response_metadata.status_code != 200:
                    print(f"Annotations status code: {response_metadata.status_code}. Retrying...\n", flush=True)
                else:
                    break
            except Exception as e:
                print(f"Exception occurred when trying to retrieve {ID} annotations:\n{e}\nRetrying...\n", flush=True)            
            time.sleep(5) # wait 5 seconds in-between tries
        print(f"Annotations retrieved for {ID}.\n", flush=True)
        
        return response_metadata.json()   
    
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

    infile = find_path(args.infile, action="r")    
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
            #infile_df = fasta_to_df(infile)
            #infile_df = remove_query(infile_df)
            #for p in infile_df.index.values:
                #if type(infile_df.loc[p, "Accession"]) == pd.Series:
                    ## only use the first accession that IS NOT the query
                    #acc = infile_df.loc[p, "Accession"][0][1:]
                    #print(f"Duplicate accessions for ID: {p}.\nProceeding with {acc}.\n", flush=True)
                #else:
                    #acc = infile_df.loc[p, "Accession"][1:]
                #proteins.append((p, acc.split(" ")[0])) # cannot handle when something is labeled "QUERY", but has the same uniprot ID  
    else:
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
    
    out_directory = find_path(args.out_directory, action="w")
    print(f"Storing outputs in {out_directory}\n", flush=True)
    
    features_df_list = []
    for p in proteins:
        prot = p[0]
        whole_prot = p[1]
        metadata = api_tools.get_metadata(prot)

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