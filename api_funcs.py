"""
Helper functions for retrieving information.

Extra imports:
    requests
        https://pypi.org/project/requests/
    helpers
        https://github.com/grtakaha/protein_alignment_tool/blob/main/helpers.py

Functions:
    get_ftp(str, str, str) -> str or None
    verify_sprot()
    get_fasta(str) -> str
    blast(str, str, str, str)
    align(str, str, str, str)
    get_metadata(str) -> dict
"""

# time.sleep() pauses are placed in an attempt to comply with usage guidelines:
# https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html#rest
# https://ebi-biows.gitdocs.ebi.ac.uk/documentation/#fair-use-policy

# TODO: Consider condensing "ID"-type requests into one function.
# TODO: Fix potential endless loop if no exception is raised and status != 200.
# TODO: Make it so that non-UniProt entries don't run infinitely.

import os
import time
import sys
import subprocess
import shutil
import urllib
from contextlib import closing
import gzip
import requests
from helpers import find_path

def get_ftp(filename, out_dir, action="save"):
    """
    Retrieves a UniProt file with the given filename.

        Parameters:
            filename (str): A filename present in UniProt's
                            FTP server (ftp://ftp.uniprot.org), stored at:
                            /pub/databases/uniprot/current_release/knowledgebase/complete/
                            Ex. reldate.txt, README, uniprot_sprot.fasta.gz
            out_dir (str):  Directory for storage of retrieved file.
            action (str):   "save" or "read"
                            save - Saves retrieved file as [OUT_DIR]/[FILENAME].
                            read - Returns retrieved file as text.

        Returns:
            content (str): String of content from the given FTP request.

        Outputs:
            One directory (SwissProt) with the following files:
                README - Swiss-Prot README.
                reldate.txt - UniProt database release dates.
                uniprot_sprot.fasta.gz - Compressed Swiss-Prot database.
                uniprot_sprot.fasta - Decompressed Swiss-Prot database.
                uniprot_sprot.pdb, .phr, .pin, .pot, .psq, .ptf, .pto
                    - Files generated by makeblastdb.
    """

    # Base URL.
    url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/"

    # Rerun the request until it returns 200.
    current_request = f"Retrieving {filename}"
    # Content will remain None if action != "read".
    content = None

    while True:
        try:
            # This does not check for status code.
            # See if Python3 urllib has a way to do this.
            with closing(urllib.request.urlopen(f"{url}/{filename}")) as response:
                print(f"{filename} found.\n", flush=True)
                if action == "save":
                    # Content is just saved at the specified location.
                    with open(f"{out_dir}/{filename}", 'wb') as ftp_file:
                        shutil.copyfileobj(response, ftp_file)
                elif action == "read":
                    # Content is decoded and read into content variable.
                    # Assume utf-8; response.headers.get_content_charset() returned None.
                    content = response.read().decode("utf-8")
                break
        except urllib.error.URLError as urlerr:
            print(f"{current_request} caused a URL error: " +
                  f"{urlerr}. Exiting...\n", flush=True)
            sys.exit()

    # Will return either None or the content of the ftp request.
    return content

def verify_sprot():
    """
    Checks for the existence of and/or downloads Swiss-Prot files.

        Outputs:
            One directory (SwissProt) with the following files:
                README - Swiss-Prot README.
                reldate.txt - UniProt database release dates.
                uniprot_sprot.fasta.gz - Compressed Swiss-Prot database.
                uniprot_sprot.fasta - Decompressed Swiss-Prot database.
                uniprot_sprot.pdb, .phr, .pin, .pot, .psq, .ptf, .pto
                    - Files generated by makeblastdb.
    """

    sprot_path = find_path(f"{os.path.abspath(os.path.dirname(__file__))}/SwissProt/",
                           "w", "d")
    sprot_exists = True
    current_reldate = get_ftp("reldate.txt", sprot_path, action="read")

    # Checks for .fasta, rather than .fasta.gz.
    # .fasta will be used for makeblastdb.
    if os.path.isfile(f"{sprot_path}/uniprot_sprot.fasta"):
        print("Swiss-Prot database found. Checking release...", flush=True)

        # Currently makes sure the entire file matches (not just Swiss-Prot).
        if os.path.isfile(f"{sprot_path}/reldate.txt"):
            with open(f"{sprot_path}/reldate.txt", "r", encoding="utf-8") as old_reldate:
                if old_reldate.read() == current_reldate:
                    print("Release dates are current.", flush=True)
                else:
                    print("Release dates are not current.\n",
                          flush=True)
                    sprot_exists = False
        else:
            print("No release dates found.\n",
                  flush=True)
            sprot_exists = False
    else:
        print("No Swiss-Prot database found.\n",
              flush=True)
        sprot_exists = False

    blastdb_exists = True
    # Download and unpack Swiss-Prot files.
    if not sprot_exists:
        print("Downloading latest Swiss-Prot release.\n", flush=True)

        get_ftp("reldate.txt", sprot_path, action="save")
        get_ftp("README", sprot_path, action="save")
        get_ftp("uniprot_sprot.fasta.gz", sprot_path, action="save")

        # Immediately unpack uniprot_sprot.fasta.gz.
        with gzip.open(f"{sprot_path}/uniprot_sprot.fasta.gz", 'rb') as sprot_gz:
            with open(f"{sprot_path}/uniprot_sprot.fasta", "wb") as sprot_fasta:
                sprot_fasta.write(sprot_gz.read())

        print(f"Swiss-Prot files saved to: {sprot_path}", flush=True)

        blastdb_exists = False

    else:
        print("Checking for Swiss-Prot BLAST database...", flush=True)

        # Function check_output returns bytes, not str.
        output = subprocess.check_output(["blastdbcmd",
                        "-db", f"{sprot_path}/uniprot_sprot.fasta",
                        "-info"]).decode("utf-8")
        if output.startswith("BLAST Database error"):
            print("No BLAST database found.\n", flush=True)
            blastdb_exists = False
        else:
            print("BLAST database found for existing " +
                  "Swiss-Prot files. Continuing...\n",
                  flush=True)

    # Make BLAST database out of Swiss-Prot FASTA.
    if not blastdb_exists:
        print("Creating BLAST database...\n", flush=True)
        subprocess.run(["makeblastdb",
                        "-in", f"{sprot_path}/uniprot_sprot.fasta",
                        "-dbtype", "prot",
                        "-title", "Swiss-Prot"], check=True)

    print("BLAST and Swiss-Prot files successfully saved.\n", flush=True)

def blast(infile, stype, out_prefix, num_res="5"):
    """
    BLASTs a given FASTA-formatted query against a local Swiss-Prot database.

        Parameters:
            infile (str):     FASTA query file.
            stype (str):      Sequence type.
                              "protein" for protein.
                              "dna" for nucleotide.
            out_prefix (str): Prefix for storage of output TSV.
                              Ex. [OUT_PREFIX].tsv
            num_res (str):    String integer of BLAST results to return.

        Outputs:
            BLAST results in outfmt 6 (TSV) form.

    """

    sprot_path = find_path(f"{os.path.abspath(os.path.dirname(__file__))}/SwissProt/",
                           "w", "d")
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

    subprocess.run([program,
                    "-db", f"{sprot_path}/uniprot_sprot.fasta",
                    "-query", infile,
                    "-outfmt", "6",
                    "-out", f"{out_prefix}.tsv",
                    "-num_alignments", num_res], check=True)

def get_fasta(fid):
    """
    Takes in a UniProt ID and returns that ID's FASTA-formatted sequence.

        Parameters:
            fid (str): A UniProt ID.

        Returns:
            response (str): A FASTA-formatted UniProt sequence.
    """

    url_fasta = f"https://rest.uniprot.org/uniprotkb/{fid}.fasta"

    # Rerun the request until it returns 200.
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

        # If query fails, try again after 10 seconds.
        time.sleep(10)

    print(f"FASTA found for {fid}.\n", flush=True)

    return response.text

def align(infile, stype, out_directory, title):
    """
    Aligns given FASTA file using command-line Clustal Omega.

        Parameters:
            infile (str):        FASTA file with at least three sequences.
            stype (str):         Sequence type.
                                 "protein" for protein.
                                 "dna" for DNA.
            out_directory (str): Directory for storage of Clustal outputs.
            title (str):         A title for the output alignment.

        Outputs:
            An alignment (.clustal_num) and percent identity matrix (.pim)
            created from the input FASTA file.

    """

    print(f"Aligning {title}...\n", flush=True)
    stype_conv = {"protein":"Protein", "dna":"DNA", "rna":"RNA"}
    subprocess.run(["clustalo",
                    "--infile", infile,
                    "--outfile", f"{out_directory}/{title}.clustal",
                    "--seqtype", stype_conv[stype],
                    "--distmat-out", f"{out_directory}/{title}.pim",
                    "--percent-id", "--full",
                    "--outfmt", "clu", "--force"], check=True)

def get_metadata(mid):
    """
    Takes in a UniProt ID and returns that ID's JSON metadata.

        Parameters:
            mid (str): A UniProt ID.

        Returns:
            response (str): A JSON of the given ID's metadata.
    """

    # TODO: Consider returning only annotations.
    # Returns entire JSON.
    url_metadata = f"https://rest.uniprot.org/uniprotkb/{mid}"
    print(f"Retrieving metadata for {mid}.", flush=True)
    # Set headers to accept JSON.
    headers = {"Accept": "application/json"}

    # For this request specfically, terminate if the first request is not 200.
    # This is to avoid sending requests for non-UniProt accessions.
    current_request = "Metadata retrieval"
    while True:
        try:
            response_metadata = requests.get(url_metadata, headers)
            if response_metadata.status_code != 200:
                print(f"{current_request} status code: " +
                      f"{response_metadata.status_code}.\n" +
                      f"Failed to retrieve metadata for {mid}. Continuing...\n",
                      flush=True)
                return None # Handle Nonetype in script that calls this.
            else:
                break

        except requests.exceptions.ConnectionError as errc:
            print(f"{current_request} caused a connection error: " +
                  f"{errc}. Retrying...\n", flush=True)
        except requests.exceptions.RequestException as err:
            print(f"{current_request} caused exception: {err}. Exiting...\n",
                  flush=True)
            sys.exit()

        # If query fails, try again after 10 seconds.
        time.sleep(10)

    print(f"Metadata retrieved for {mid}.\n", flush=True)

    return response_metadata.json()
