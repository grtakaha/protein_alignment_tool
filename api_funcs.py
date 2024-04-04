"""
Helper functions for sending various API requests.

Extra imports:
    requests
        https://pypi.org/project/requests/

Functions:
    blast(str, str, str, str, str, str) -> str
    get_blast_results(str) -> (str, str)
    get_fasta(str) -> str
    align(str, str, str, str) -> str
    get_alignment(str) -> str
    get_pim(str) -> str
    get_metadata(str) -> dict
"""

##### CONDENSE "ID" KIND OF REQUESTS INTO ONE FUNCTION #####
# TODO: Limit the number of requests that can be submitted at once.
# TODO: Get this to work with CLI tools instead; run locally.

import time
import sys
import requests

##### SET UP SO THAT IT WILL BLAST MULTIPLE QUERIES AT ONCE AND PARSE #####
# acc and seq do no have newlines
def blast(email, stype, acc, seq, num_res="10"):
    """
    Takes in BLAST parameters and returns a BLAST job ID against UniProt databases.

        Parameters:
            email (str):   Personal email used for BLAST runs.
            program (str): Type of BLAST to run.
                           "blastx" for nucleotide query
                           "blastp" for protein query
            stype (str):   Sequence type.
                           "protein" for protein.
                           "dna" for nucleotide.
            acc (str):     A sequence accession, processed in UniProt style.
                           >sp|P00784|PAPA1_CARPA -> PAPA1_CARPA
            seq (str):     As of 2024.03.22, a protein sequence.
            num_res (str): Number of BLAST results to return.
                           This is a string integer.

        Returns:
            response (str): A BLAST submission ID.
    """

    # Ensure that queries are spaced out.
    time.sleep(10)

    # Variable stype is tied to program.
    if stype == "dna":
        program = "blastx"
    elif stype == "protein":
        program = "blastp"
    else:
        print("Given stype not found. Please specify " +
              "\"-stype dna\" for blastx OR \"-stype protein\" for blastp",
              flush=True)
        print("Defaulting to stype=\"protein\" and program=\"blastx\"\n", flush=True)
        stype = "protein"
        program = "blastp"

    # Variable acc is "accession", without ">"
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
    current_request = "BLAST query"
    while True:
        try:
            response = requests.post(url, data=data)
            if response.status_code != 200:
                print(f"{current_request} status code: " +
                      f"{response.status_code}. Retrying...\n", flush=True)
            else:
                print(f"{acc} successfully BLASTed.\n", flush=True)
                break

        except requests.exceptions.ConnectionError as errc:
            print(f"{current_request} caused a connection error: " +
                  f"{errc}. Retrying...\n", flush=True)
        except requests.exceptions.RequestException as err:
            print(f"{current_request} caused exception: {err}. Exiting...\n",
                  flush=True)
            sys.exit()

        time.sleep(10) # wait 5 second in-between tries

    # Print the response
    return response.text

def get_blast_results(bid):
    """
    Takes in a BLAST bid and returns a tuple of BLAST results.

        Parameters:
            bid (str): A BLAST submission bid.

        Returns:
            response (tup): A tuple of human readable and tsv BLAST results.
    """
    # Ensure that queries are spaced out.   
    time.sleep(10)
    
    url_out = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{bid}/out"
    url_tsv = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{bid}/tsv"

    # rerun the request until it returns 200
    # it's possible this becomes an endless loop...deal with that later
    current_request = "Default BLAST retrieval"
    while True:
        try:
            response_out = requests.get(url_out)
            if response_out.status_code != 200:
                print(f"{current_request} output status code: " +
                      f"{response_out.status_code}. Retrying...\n",
                      flush=True)
            else:
                print(f"{current_request} output found for {bid}.\n", flush=True)
                break

        except requests.exceptions.ConnectionError as errc:
            print(f"{current_request} caused a connection error: " +
                  f"{errc}. Retrying...\n", flush=True)
        except requests.exceptions.RequestException as err:
            print(f"{current_request} caused exception: {err}. Exiting...\n",
                  flush=True)
            sys.exit()

        time.sleep(60) # wait 5 seconds in-between tries

    time.sleep(10)
    # rerun the request until it returns 200
    # it's possible this becomes an endless loop...deal with that later
    current_request = "TSV BLAST retrieval"
    while True:
        try:
            response_tsv = requests.get(url_tsv)
            if response_tsv.status_code != 200:
                print(f"{current_request} output status code: " +
                      f"{response_tsv.status_code}. Retrying...\n",
                      flush=True)
            else:
                print(f"{current_request} output found for {bid}.\n", flush=True)
                break

        except requests.exceptions.ConnectionError as errc:
            print(f"{current_request} caused a connection error: " +
                  f"{errc}. Retrying...\n", flush=True)
        except requests.exceptions.RequestException as err:
            print(f"{current_request} caused exception: {err}. Exiting...\n",
                  flush=True)
            sys.exit()

        time.sleep(60) # wait 5 seconds in-between tries

    return (response_out.text, response_tsv.text) # tuple of strings

def get_fasta(fid):
    """
    Takes in a UniProt fid and returns that fid's FASTA-formatted sequence.

        Parameters:
            fid (str): A UniProt fid.

        Returns:
            response (str): A FASTA-formatted UniProt sequence.
    """

    url_fasta = f"https://rest.uniprot.org/uniprotkb/{fid}.fasta"

    # rerun the request until it returns 200
    # it's possible this becomes an endless loop...deal with that later
    current_request = "FASTA retrieval"
    while True:
        try:
            response = requests.get(url_fasta)
            if response.status_code != 200:
                print(f"{current_request} status code: " +
                      f"{response.status_code}. Retrying...\n", flush=True)
            else:
                break

        except requests.exceptions.ConnectionError as errc:
            print(f"{current_request} caused a connection error: " +
                  f"{errc}. Retrying...\n", flush=True)
        except requests.exceptions.RequestException as err:
            print(f"{current_request} caused exception: {err}. Exiting...\n",
                  flush=True)
            sys.exit()

        time.sleep(10) # wait 5 seconds in-between tries
    print(f"FASTA found for {fid}.\n", flush=True)

    return response.text

def align(email, stype, title, seqs):
    """
    Takes in Clustal Omega parameters and returns a Clustal job ID.

        Parameters:
            email (str): Personal email used for Clustal runs.
            stype (str): Sequence type.
                         "protein" for protein.
                         "dna" for nucleotide.
            title (str): A title for the output alignment.
            seqs (str):  A set of sequences to be aligned.

        Returns:
            response (str): A Clsutal Omega submission ID.
    """

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
    current_request = "Alignment submission"
    while True:
        try:
            response = requests.post(url, data=data)
            if response.status_code != 200:
                print(f"{current_request} status code: " +
                      f"{response.status_code}. Retrying...\n", flush=True)
            else:
                break

        except requests.exceptions.ConnectionError as errc:
            print(f"{current_request} caused a connection error: " +
                  f"{errc}. Retrying...\n", flush=True)
        except requests.exceptions.RequestException as err:
            print(f"{current_request} caused exception: {err}. Exiting...\n",
                  flush=True)
            sys.exit()

        time.sleep(10) # wait 5 second in-between tries

    print(f"Alignment for {title} submitted with ID: {response.text}.\n", flush=True)
    # Print the response
    return response.text

def get_alignment(aid):
    """
    Takes in a Clustal Omega aid and returns a Clustal alignment.

        Parameters:
            aid (str): A Clustal Omega submission aid.

        Returns:
            response (str): A Clustal alignment.
    """

    url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{aid}/aln-clustal_num"

    # rerun the request until it returns 200
    # it's possible this becomes an endless loop...deal with that later
    current_request = "Alignment retrieval"
    while True:
        try:
            response = requests.get(url)
            if response.status_code != 200:
                print(f"{current_request} status code: " +
                      f"{response.status_code}. Retrying...\n", flush=True)
            else:
                break

        except requests.exceptions.ConnectionError as errc:
            print(f"{current_request} caused a connection error: " +
                  f"{errc}. Retrying...\n", flush=True)
        except requests.exceptions.RequestException as err:
            print(f"{current_request} caused exception: {err}. Exiting...\n",
                  flush=True)
            sys.exit()

        time.sleep(10) # wait 5 seconds in-between tries
    print(f"Alignment found for {aid}.\n", flush=True)

    return response.text

def get_pim(pid):
    """
    Takes in a Clustal Omega pid and returns a percent pidentity matrix (PIM).

        Parameters:
            pid (str): A Clustal Omega submission pid.

        Returns:
            response (str): A PIM based on the given pid's alignment.
    """

    url = f"https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{pid}/pim"

    # rerun the request until it returns 200
    # it's possible this becomes an endless loop...deal with that later
    current_request = "PIM retrieval"
    while True:
        try:
            response = requests.get(url)
            if response.status_code != 200:
                print(f"{current_request} status code: " +
                      f"{response.status_code}. Retrying...\n", flush=True)
            else:
                break

        except requests.exceptions.ConnectionError as errc:
            print(f"{current_request} caused a connection error: " +
                  f"{errc}. Retrying...\n", flush=True)
        except requests.exceptions.RequestException as err:
            print(f"{current_request} caused exception: {err}. Exiting...\n",
                  flush=True)
            sys.exit()

        time.sleep(10) # wait 5 seconds in-between tries
    print(f"PIM found for {pid}.\n", flush=True)

    return response.text

def get_metadata(mid):
    """
    Takes in a UniProt mid and returns that mid's JSON metadata.

        Parameters:
            mid (str): A UniProt mid.

        Returns:
            response (str): A JSON of the given mid's metadata.
    """

    # for now, just returns entire json...consmider extracting annotations only
    url_metadata = f"https://rest.uniprot.org/uniprotkb/{mid}"

    # rerun the request until it returns 200
    # it's possible this becomes an endless loop...deal with that later
    # Set the headers to accept JSON
    headers = {"Accept": "application/json"}
    current_request = "Metadata retrieval"
    while True:
        try:
            response_metadata = requests.get(url_metadata, headers)
            if response_metadata.status_code != 200:
                print(f"{current_request} status code: " +
                      f"{response_metadata.status_code}. Retrying...\n", flush=True)
            else:
                break

        except requests.exceptions.ConnectionError as errc:
            print(f"{current_request} caused a connection error: " +
                  f"{errc}. Retrying...\n", flush=True)
        except requests.exceptions.RequestException as err:
            print(f"{current_request} caused exception: {err}. Exiting...\n",
                  flush=True)
            sys.exit()

            time.sleep(10) # wait 5 seconds in-between tries
    print(f"Metadata retrieved for {mid}.\n", flush=True)

    return response_metadata.json()
