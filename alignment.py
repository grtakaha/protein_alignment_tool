import argparse
import requests
import os
import time
import pandas as pd
import json

class api_tools():
    # acc and seq do no have newlines
    def align(email, stype, title, seqs):
        # for now, acc is "accession", no ">"
        print(f"Aligning {title}...\n", flush=True)
        # Define the data to be sent in the form
        data = {
            'email': email,
            'outfmt': 'clustal_num',
            'order': 'aligned',
            'title': title,
            'sequence': seqs, # will be in FASTA format: f'>{acc}\n{seq}\n>{acc}\n{seq}\n'
            'stype': stype
        }
        
        # Define the URL
        url = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run"
        
        # rerun the request until it returns 200
        # it's possible this becomes an endless loop...deal with that later
        while True:
            try:
                response = requests.post(url, data=data)
                if response.status_code != 200:
                    print(f"Alignment submission status code: {response.status_code}. Retrying...\n", flush=True)
                else:
                    break
            except Exception as e:
                print(f"Alignment submission caused exception: {e}. Retrying...\n", flush=True)
            time.sleep(5) # wait 5 second in-between tries
            
        print(f"Alignment for {title} submitted with ID: {response.text}.\n", flush=True)
        # Print the response
        return response.text
    
    def get_alignment(ID):
        url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{ID}/aln-clustal_num"
        
        # rerun the request until it returns 200
        # it's possible this becomes an endless loop...deal with that later
        while True:
            try:
                response = requests.get(url)
                if response.status_code != 200:
                    print(f"Alignment retrieval status code: {response.status_code}. Retrying...\n", flush=True)
                else:
                    break
            except Exception as e:
                print(f"Exception occurred when trying to retrieve {ID}:\n{e}\nRetrying...\n", flush=True)
            time.sleep(5) # wait 5 seconds in-between tries
        print(f"Alignment found for {ID}.\n", flush=True)
        
        return response.text
    
    def get_pim(ID):
        url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{ID}/pim"
        
        # rerun the request until it returns 200
        # it's possible this becomes an endless loop...deal with that later
        while True:
            try:
                response = requests.get(url)
                if response.status_code != 200:
                    print(f"PIM retrieval status code: {response.status_code}. Retrying...\n", flush=True)
                else:
                    break
            except Exception as e:
                print(f"Exception occurred when trying to retrieve PIM for {ID}:\n{e}\nRetrying...\n", flush=True)
            time.sleep(5) # wait 5 seconds in-between tries
        print(f"PIM found for {ID}.\n", flush=True)
        
        return response.text
    
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
    
    infile = find_path(args.infile, action="r")
    print(f"Processing sequences from {infile} \n", flush=True)
    with open(infile, "r") as f:
        seqs = "".join(f.readlines())
    
    out_directory = find_path(args.out_directory, action="w")
    print(f"Storing outputs in {out_directory}\n", flush=True)
    
    stype = args.stype
    if stype not in ["dna", "protein"]:
        print("Given stype not found. Please specify \"-stype dna\" OR \"-stype protein\n", flush=True)
        print("Defaulting to stype=\"protein\"\n", flush=True)
        stype = "protein"  
    
    title = args.title
    
    ID = api_tools.align(args.email, stype, title, seqs)
    alignment = api_tools.get_alignment(ID)
    pim = api_tools.get_pim(ID)
    
    alignment_file = f"{out_directory}/{title}.clustal_num"
    pim_file = f"{out_directory}/{title}.pim"
    
    print(f"Saving {title} alignment to {alignment_file}.", flush=True)
    with open(alignment_file, "w") as f:
        f.write(alignment) # just overwriting it if it exists
    with open(pim_file, "w") as p:
        p.write(pim) # just overwriting if it exists
        
if __name__ == "__main__":
    main()