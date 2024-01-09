import argparse
import requests
import os
import time
import pandas as pd
import json
import uniprot_fastas as uf
import alignment
import clustal_to_svg as cts

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

def arguments():
    parser = argparse.ArgumentParser()
    
    # fix to maintain newlines...later...low priority
    parser.add_argument("-module", nargs="?", default="all",
                        help='''Module(s) to run.
                        (Note that different modules require different inputs.)
                        uniprot: Run uniprot_fastas.py. Input: protein .fasta file (minimum 1 sequence).
                        align: Run alignment.py. Input: protein .fasta file (minimum 3 sequences).
                        clust_to_svg: Run clustal_to_svg.py. Input: .clustal_num or .clustal alignment.
                        uniprot_align: Run uniprot_fastas.py, then alignment.py. Input: protein .fasta file (minimum 1 sequence).
                        align_svg: Run alignment.py, then clustal_to_svg.py. Input: protein .fasta file (minimum 3 sequences).
                        all: Run all of the above, in order. Input: .fasta file (minimum 1 sequence).''')
    parser.add_argument("-infile", help="Input file. Different modules require different inputs.")
    parser.add_argument("-out_directory", default=None) # if None, set to infile directory later(?)
    parser.add_argument("-mod_options", help=
                        '''Module-specific options, in nested Python dictionary format, necessary to run each module.
                        If a file is passed to -control_file then this option will not be used.
                        Nested Python dictionary format: {module1:{option1:o1,option2:o2},module2:{option1:o1,option2:o2,option3:o3}}
                        The quotes in the following examples are important!
                        Ex1. -mod_options "{uniprot:{stype:protein,email:user@user.com,num_res:10},align:{stype:protein,email:user@user.com,title:project1}}"
                        Ex2. -mod_options "{clust_to_svg:{codes:TRUE,nums:FALSE,uniprot_format:TRUE}}"''')
    parser.add_argument("-control_file", default=None, help=
                        "Full file path to a control file with text in the format of -mod_options.")
    
    args = parser.parse_args()
    
    options = {}
    options["module"] = args.module
    options["infile"] = find_path(args.infile, action="r")
    if not args.out_directory:
        options["out_directory"] = "/".join(options["infile"].split("/")[:-1])
    else:
        options["out_directory"] = find_path(args.out_directory, action="w")
    options["out_directory"] = find_path(args.out_directory, action="w")
    options["mod_options"] = args.mod_options
    options["control_file"] = args.control_file
    
    return options

def main(options):
    # if running on its own, will use command-line arguments; otherwise, pass arguments for this function
    # consider renaming
    
    infile = options["infile"]
    print(f"Processing sequences from {infile}\n", flush=True)
    with open(infile, "r") as f:
        seqs = "".join(f.readlines())
    
    # in the out_directory, for cleanliness' sake and ease of access without passing outputs from one module to another, save in a different directory here
    out_directory = options["out_directory"] + "/output/"
    print(f"Storing outputs in {out_directory}\n", flush=True)
    
    if options["control_file"]:
        control_file = find_path(options["control_file"], action="r")
        with open(control_file, "r") as cf:
            # for now be in the same format as command-line options; consider making more human-readable
            opt_str = cf.readline().replace("\'", "\"")
    else:
        opt_str = options["mod_options"]
    
    # makes readable for json.loads
    opt_str = opt_str.replace(" ", "").replace("{", "{\"").replace(":", "\":\"").replace("}", "\"}").replace(",", "\",\"").replace(":\"{", ":{").replace("}\"}", "}}").replace("}\",", "},")
    mod_options = json.loads(opt_str.replace("\"\"", "\""))
    
    # set universal options (may change later)
    universal_options = {"infile":infile, "out_directory":out_directory}
    
    module = options["module"]
    file_list = [(infile, out_directory)] # infile list that can be set for each protein after uniprot_fastas.py
    if module in ["uniprot", "uniprot_align", "all"]:
        # make sure to reset infile here for uniprot_align and all
        opts = {**mod_options["uniprot"], **universal_options}
        print(f"Running uniprot_fastas.py with the following options: {opts}", flush=True)
        uf.main(opts)
        
        # read all resulting folders in current out_directory; and save full file paths as a list; these will be the new infiles for alignment
        # set equal to list of infiles for next thing....maybe...need to fix
        # os.walk() is a generator
        # next() runs it just once, retrieving the subdirectories (for each protein) created by uniprot_fastas.py
        # The output is a list: ["DIRECTORY", [SUBDIR1, SUBDIR2, ..., SUBDIRN], [FILE1, FILE2, ..., FILEN]]
        # Note that this doesn't give full subdirectory paths; just names
        file_list = [] # reset if uniprot_fastas.py is run
        for d in next(os.walk(out_directory))[1]:
            file_list.append((f"{out_directory}/{d}/all.fasta", f"{out_directory}/{d}/")) # the .fasta file of seqs will be stored in this manner
    if module in ["align", "uniprot_align", "align_svg", "all"]:
        # make sure to reset infile here for align_svg and all
        for infile, out_directory in file_list:
            universal_options["infile"] = infile
            universal_options["out_directory"] = out_directory
            opts = {**mod_options["align"], **universal_options}
            print(f"Running alignment.py with the following options: {opts}", flush=True)
            alignment.main(opts)
            
        for i, ftup in enumerate(file_list):
            prot_dir = "/".join(ftup[0].split("/")[:-1])
            title = opts["title"]
            file_list[i] = (f"{prot_dir}/{title}.clustal_num", ftup[1])
    
    if module in ["clust_to_svg", "align_svg", "all"]:
        for infile, out_directory in file_list:
            universal_options["infile"] = infile
            universal_options["out_directory"] = out_directory + "/svg_alignments/" # hardcoded...consider allowing options later
            opts = {**mod_options["clust_to_svg"], **universal_options}
            print(f"Running clustal_to_svg.py with the following options: {opts}", flush=True)
            cts.main(opts)
            
if __name__ == "__main__":
    main(arguments())