# set of api functions
import requests
import time
import json

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
