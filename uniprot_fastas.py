import argparse
import requests
import os
import time
import pandas as pd
import json

class api_tools():
    # acc and seq do no have newlines
    def blast(email, program, stype, acc, seq, num_res="250"):
        # for now, acc is "accession", no ">"
        print(f"BLASTing {acc}...\n", flush=True)
        # Define the data to be sent in the form
        data = {
            'email': email,
            'program': program,
            'matrix': 'BLOSUM62',
            'alignments': num_res, # for now, make these the same
            'scores': num_res,
            'exp': '10',
            'filter': 'F',
            'gapalign': 'true',
            'compstats': 'F',
            'align': '0',
            'stype': stype,
            'sequence': f'>{acc}\n{seq}',
            'database': 'uniprotkb_refprotswissprot'
        }
        
        # Define the URL
        url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run"
        
        # rerun the request until it returns 200
        # it's possible this becomes an endless loop...deal with that later
        while True:
            try:
                response = requests.post(url, data=data)
                if response.status_code != 200:
                    print(f"BLAST query status code: {response.status_code}. Retrying...\n", flush=True)
                else:
                    print(f"{acc} successfully BLASTed.\n", flush=True)
                    break
            except Exception as e:
                print(f"BLAST query caused exception: {e}. Retrying...\n", flush=True)
            time.sleep(5) # wait 5 second in-between tries
        
        # Print the response
        return response.text
    
    def get_blast_results(ID):
        url_out = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{ID}/out"
        url_tsv = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{ID}/tsv"
        
        # rerun the request until it returns 200
        # it's possible this becomes an endless loop...deal with that later
        while True:
            try:
                response_out = requests.get(url_out)
                if response_out.status_code != 200:
                    print(f"Default BLAST output status code: {response_out.status_code}. Retrying...\n", flush=True)
                else:
                    print(f"Default BLAST output found for {ID}.\n", flush=True)
                    break
            except Exception as e:
                print(f"Default BLAST caused exception: {e}. Retrying...\n", flush=True)
            time.sleep(5) # wait 5 seconds in-between tries
           
        # rerun the request until it returns 200
        # it's possible this becomes an endless loop...deal with that later
        while True:
            try:
                response_tsv = requests.get(url_tsv)
                if response_tsv.status_code != 200:
                    print(f"TSV BLAST output status code: {response_tsv.status_code}. Retrying...\n", flush=True)
                else:
                    print(f"TSV BLAST output found for {ID}.\n", flush=True)
                    break
            except Exception as e:
                print(f"TSV BLAST caused exception: {e}. Retrying...\n", flush=True)
            time.sleep(5) # wait 5 seconds in-between tries
            
        return (response_out.text, response_tsv.text) # tuple of strings
    
    def get_fasta(ID):
        url_fasta = f"https://rest.uniprot.org/uniprotkb/{ID}.fasta"
        
        # rerun the request until it returns 200
        # it's possible this becomes an endless loop...deal with that later
        while True:
            try:
                response = requests.get(url_fasta)
                if response.status_code != 200:
                    print(f"FASTA status code: {response.status_code}. Retrying...\n", flush=True)
                else:
                    break
            except Exception as e:
                print(f"Exception occurred when trying to retrieve {ID} FASTA:\n{e}\nRetrying...\n", flush=True)
            time.sleep(5) # wait 5 seconds in-between tries
        print(f"FASTA found for {ID}.\n", flush=True)
        
        return response.text
    
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

def arguments():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-infile")
    parser.add_argument("-out_directory", default="./test_runs/test/")
    parser.add_argument("-stype", default="protein") # program is tied to stype (dna:blastx, rna:blastx, protein:blastp)
    parser.add_argument("-email")
    parser.add_argument("-num_res", default="10")
    
    args = parser.parse_args()
    
    options = {}
    options["infile"] = args.infile
    options["out_directory"] = args.out_directory
    options["stype"] = args.stype
    options["email"] = args.email
    options["num_res"] = args.num_res
    
    return options   
    
def main(options):

    infile = find_path(options["infile"], action="r")    
    print(f"Processing sequences from {infile}\n", flush=True)
    infile_df = fasta_to_df(infile)
    
    out_directory = find_path(options["out_directory"], action="w")
    print(f"Storing outputs in {out_directory}\n", flush=True)
    
    stype = options["stype"]
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
        
        blast_id = api_tools.blast(options["email"], program, stype, protein, sequence, num_res=options["num_res"])
        blast_dict[protein] = blast_id
        
        prot_directory = find_path(f"{out_directory}/{protein}/", action="w")     
        with open(f"{prot_directory}/all.fasta", "w") as all_fasta:
            all_fasta.write(f"{accession[0]}QUERY_{accession[1:]}\n{sequence}\n") # just overwriting it if it exists

    for protein in blast_dict:
        # run blast_ids all at once and then retrieve all at once
        blast_results = api_tools.get_blast_results(blast_dict[protein]) # consider retrieving ids all at once, then getting results one by one after running
        
        # prot_directory should already be created; no find_path necessary
        prot_directory = f"{out_directory}/{protein}/"
        
        with open(f"{prot_directory}/{protein}.out", "w") as out:
            out.write(blast_results[0])
        with open(f"{prot_directory}/{protein}.tsv", "w") as tsv:
            tsv.write(blast_results[1])         
        
        # parse blast_results[1] (tsv form)
        for line in blast_results[1].split("\n")[1:-1]: # skip first line (headers) and last line (empty); don't bother reading into a dataframe
            hit = line.split("\t")[2]
            print(f"Found BLAST hit: {hit}\n", flush=True)
            hit_fasta = api_tools.get_fasta(hit)
            hit_filename = hit_fasta.split("\n")[0].split(" ")[0].split("|")[-1]
            with open(f"{prot_directory}/{hit_filename}.fasta", "w") as fasta:
                fasta.write(hit_fasta)
            hit_metadata = api_tools.get_metadata(hit)
            with open(f"{prot_directory}/{hit_filename}.json", "w") as json:
                json.write(str(hit_metadata))
            with open(f"{prot_directory}/all.fasta", "a") as all_fasta:
                all_fasta.write(hit_fasta)

            features = hit_metadata.get("features")
            if features is not None:
                features_df = pd.json_normalize(features)
            else:
                print(f"Could not retrieve features for {hit}...continuing with empty DataFrame\n", flush=True)
                features_df = pd.DataFrame(columns=["type","description","evidences","location.start.value","location.start.modifier","location.end.value","location.end.modifier"])
    
            features_df.to_csv(f"{prot_directory}/{hit_filename}.ann", sep="\t")
        
if __name__ == "__main__":
    main(arguments())