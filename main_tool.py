import argparse
import requests
import os
import time
import pandas as pd
import json
import alignment
import clustal_to_svg as cts
import subprocess

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

#def arguments():
    #parser = argparse.ArgumentParser()
    
    ## fix to maintain newlines...later...low priority
    #parser.add_argument("-module", nargs="?", default="all",
                        #help='''Module(s) to run.
                        #(Note that different modules require different inputs.)
                        #uniprot: Run uniprot_fastas.py. Input: protein .fasta file (minimum 1 sequence).
                        #align: Run alignment.py. Input: protein .fasta file (minimum 3 sequences).
                        #clust_to_svg: Run clustal_to_svg.py. Input: .clustal_num or .clustal alignment.
                        #uniprot_align: Run uniprot_fastas.py, then alignment.py. Input: protein .fasta file (minimum 1 sequence).
                        #align_svg: Run alignment.py, then clustal_to_svg.py. Input: protein .fasta file (minimum 3 sequences).
                        #all: Run all of the above, in order. Input: .fasta file (minimum 1 sequence).''')
    #parser.add_argument("-infile", help="Input file. Different modules require different inputs.")
    #parser.add_argument("-out_directory", default=None) # if None, set to infile directory later(?)
    #parser.add_argument("-mod_options", help=
                        #'''Module-specific options, in nested Python dictionary format, necessary to run each module.
                        #If a file is passed to -control_file then this option will not be used.
                        #Nested Python dictionary format: {module1:{option1:o1,option2:o2},module2:{option1:o1,option2:o2,option3:o3}}
                        #The quotes in the following examples are important!
                        #Ex1. -mod_options "{uniprot:{stype:protein,email:user@user.com,num_res:10},align:{stype:protein,email:user@user.com,title:project1}}"
                        #Ex2. -mod_options "{clust_to_svg:{codes:TRUE,nums:FALSE,uniprot_format:TRUE}}"''')
    #parser.add_argument("-control_file", default=None, help=
                        #"Full file path to a control file with text in the format of -mod_options.")
    
    #args = parser.parse_args()
    
    #options = {}
    #options["module"] = args.module
    #options["infile"] = find_path(args.infile, action="r")
    #if not args.out_directory:
        #options["out_directory"] = "/".join(options["infile"].split("/")[:-1])
    #else:
        #options["out_directory"] = find_path(args.out_directory, action="w")
    #options["out_directory"] = find_path(args.out_directory, action="w")
    #options["mod_options"] = args.mod_options
    #options["control_file"] = args.control_file
    
    #return options


# main_tool.py

def parse_args():
    parser = argparse.ArgumentParser(description='Centralized Tool Manager')
    # I don't completely understand, but -i and -o must come before tool_name
    # -i ./infile.fasta -o ./out_dir annotate
    
    parser.add_argument("-i", "--infile", type=str, help="Full path of input file")
    parser.add_argument("-o", "--out_directory", type=str, help="Full path of output directory")

    subparsers = parser.add_subparsers(dest='tool_name', help='Tool to execute') # must be after -i and -o (positional)
    
    # each tool should have its own argument parsing
    ##### CONSIDER REMOBING EMAIL REQUIREMENT ##### DOES SETTING A PLACEHOLDER WORK?
    blast = subparsers.add_parser("blast")
    blast.add_argument("-s", "--stype", default="protein") # program is tied to stype (dna:blastx, rna:blastx, protein:blastp)
    blast.add_argument("-e", "--email", default="") # see if this will run without an email
    blast.add_argument("-nr", "--num_res", default="10")
    
    annotation = subparsers.add_parser("annotate") # instead of having combination calls...consider having an argument that asks if other actions should be performed
    # only uses default -i and -o arguments for retrieve_annotations
    
    alignment = subparsers.add_parser("align")
    alignment.add_argument("-s", "--stype", default="protein")
    alignment.add_argument("-e", "--email", default="")
    alignment.add_argument("-t", "--title", default="alignment")
    
    clust_to_svg = subparsers.add_parser("svg")
    clust_to_svg.add_argument("-c", "--codes", default="FALSE")
    clust_to_svg.add_argument("-n", "--nums", default="FALSE")
    clust_to_svg.add_argument("-u", "--uniprot_format", default="FALSE")
    clust_to_svg.add_argument("-a", "--annotations", default="") # will annotate if this is provided at all
    
    
    # if multitool is set, then it should have all tools from start to finish, with arguments for each tool listed in the form:
    # --multitool tool1 arg1=thing1,arg2=thing2,arg3=thing3 tool2 arg1=thing1,arg2=thing2,arg3=thing3
    # parse the above in this script; NO SPACES
    # only call this if the tool specified is "multi", in the form:
    # python main_tool.py multi -i infile -o out_dir -ra arg1=thing1,arg2=thing2 -a arg1=thing1,arg2=thing2
        
    multi = subparsers.add_parser("multi")
    # consider requesting a multi command file
    # note that it will only run in one direction, not the order of input commands
    #multi.add_argument("-ann", "--annotate", nargs="?")
    #multi.add_argument("-uf", "--uniprot_fastas", nargs="?")
    #multi.add_argument("-a", "--align", nargs="?")
    #multi.add_argument("-cts", "--clust_to_svg", nargs="?")
    multi.add_argument("-ord", "--order", nargs="+") # will return a list of 1 or more tools to use with these arguments
    multi.add_argument("-s", "--stype", default="protein") # program is tied to stype (dna:blastx, rna:blastx, protein:blastp)
    multi.add_argument("-e", "--email", default="")
    multi.add_argument("-nr", "--num_res", default="10")
    multi.add_argument("-t", "--title", default="alignment")
    multi.add_argument("-c", "--codes", default="FALSE")
    multi.add_argument("-n", "--nums", default="FALSE")
    multi.add_argument("-u", "--uniprot_format", default="FALSE")
    multi.add_argument("-a", "--annotations", default="") # will annotate if this is provided at all #### cannot be nonetype
 
    #print(parser.parse_args(["annotate", "-i", "test", "-o", "test2"]))
    return parser.parse_args()

##### TODO: MAKE IT SO THAT FILES CAN BE FOUND ANYWHERE #####
##### TODO: MAKE SURE NON-UNIPROT QUERIES DON'T BREAK ANNOTATION #####
##### TODO: MAKE IT SO THAT ANNOTATIONS AREN'T REQUIRED ##### DONE
def execute_tool(args):
    # return output files? or have a separate function to return output files for a particular tool
    # Execute the selected tool
    if args.tool_name == "annotate": # annotation command main_tool.py ann ...
        subprocess.run(["python", f"{os.path.abspath(os.path.dirname(__file__))}/retrieve_annotations.py",
                        "--infile", args.infile,
                        "--out_directory", args.out_directory])
        ##### SOMETHING WRONG HERE #####
        args.annotations = find_outputs(args)[0] # will change args

    elif args.tool_name == "blast": # replace with split uniprot_fastas into blast and split blast...
        subprocess.run(["python", f"{os.path.abspath(os.path.dirname(__file__))}/blast.py", 
                        "--infile", args.infile, 
                        "--out_directory", args.out_directory,
                        "--stype", args.stype,
                        "--email", args.email,
                        "--num_res", args.num_res])      

    elif args.tool_name == "align":
        subprocess.run(["python", f"{os.path.abspath(os.path.dirname(__file__))}/alignment.py", 
                        "--infile", args.infile, 
                        "--out_directory", args.out_directory,
                        "--stype", args.stype,
                        "--email", args.email,
                        "--title", args.title])
        args.infile = find_outputs(args)[0] # .clustal_num file
        
    elif args.tool_name == "svg":
        #print(args)
        #print(args.annotations)
        subprocess.run(["python", f"{os.path.abspath(os.path.dirname(__file__))}/clustal_to_svg.py", 
                        "--infile", args.infile, 
                        "--out_directory", args.out_directory,
                        "--codes", args.codes,
                        "--nums", args.nums,
                        "--uniprot_format", args.uniprot_format,
                        "--annotations", args.annotations])
    
    # for now, use as if all variables are distinct...
    elif args.tool_name == "multi":
        ##### I THINK WHAT'S HAPPENING IS THAT THE FIRST LOOP - THE ONE THAT HAS blast, annotate, align, svg -
        ##### IS CONTINUING AFTER BLAST WITH THE CURRENT ARGS; NEED TO MAKE SURE TO BREAK OUT OF OUTER LOOP IF BLAST IS CALLED...
        for i, tool in enumerate(args.order): # list of one or more: blast annotate align svg
            # eventually, turn uf into blast... you can get this to do uf things with just blast capabilities
            print(f"Running {tool}.", flush=True)
            args.tool_name = tool
            print(f"POS1: {args}")
            execute_tool(args)
            print(f"POS2: {args}")
            if tool == "blast" and i+1 < len(args.order): # will only continue if there is something after blast
                # run everything from here on out using multi and different inputs
                infiles = find_outputs(args)
                print(f"INFILE LIST: {infiles}")
                ##### MAY ERROR OUT BECAUSE +1 MIGHT NOT BE IN RANGE #####
                args.order = args.order[i+1:] # start after blast call
                for f in infiles:
                    print(f"ORDER: {args.order}")
                    print(f"INFILE: {f}")
                    args.tool_name = "multi"
                    args.infile = f
                    args.out_directory = "/".join(f.split("/")[:-1]) + "/"
                    execute_tool(args)
                break # BREAKS OUT OF OUTER LOOP IF BLAST IS FOUND...PREVENTS EXTRA EXECUTION OF ANNOTATE, ALIGN, AND SVG AFTER BLAST
                
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

# instead of using this, consider returning outputs, then printing them seaprately only if this is run from main/main_tool.py
def find_outputs(args):
    new_inputs = []
    if args.tool_name == "annotate":
        # it has another output for individual annotations
        new_inputs.append(f"{args.out_directory}/all.ann") # tsv for annotations input later
    elif args.tool_name == "blast":
        proteins = fasta_to_df(find_path(args.infile, "r"))
        for p in proteins.index.values:
            new_inputs.append(f"{args.out_directory}/{p}/all.fasta") # fastas of blast hits
    elif args.tool_name == "align":
        new_inputs.append(f"{args.out_directory}/{args.title}.clustal_num")
    elif args.tool_name == "svg":
        pass # not necessary to return anything because svg is the end of the line
    return new_inputs

def main():
    args = parse_args()
    print(args)
    
    execute_tool(args)
    # add execution only after args works
    
    
    

# for now, just be ok with only being able to blast/align one thing at a time
# figure out clustered things later


#def main_archived(options):
    ## if running on its own, will use command-line arguments; otherwise, pass arguments for this function
    ## consider renaming
    
    #infile = options["infile"]
    #print(f"Processing sequences from {infile}\n", flush=True)
    #with open(infile, "r") as f:
        #seqs = "".join(f.readlines())
    
    ## in the out_directory, for cleanliness' sake and ease of access without passing outputs from one module to another, save in a different directory here
    #out_directory = options["out_directory"] + "/output/"
    #print(f"Storing outputs in {out_directory}\n", flush=True)
    
    #if options["control_file"]:
        #control_file = find_path(options["control_file"], action="r")
        #with open(control_file, "r") as cf:
            ## for now be in the same format as command-line options; consider making more human-readable
            #opt_str = cf.readline().replace("\'", "\"")
    #else:
        #opt_str = options["mod_options"]
    
    ## makes readable for json.loads
    #opt_str = opt_str.replace(" ", "").replace("{", "{\"").replace(":", "\":\"").replace("}", "\"}").replace(",", "\",\"").replace(":\"{", ":{").replace("}\"}", "}}").replace("}\",", "},")
    #mod_options = json.loads(opt_str.replace("\"\"", "\""))
    
    ## set universal options (may change later)
    #universal_options = {"infile":infile, "out_directory":out_directory}
    
    #module = options["module"]
    #file_list = [(infile, out_directory)] # infile list that can be set for each protein after uniprot_fastas.py
    #if module in ["uniprot", "uniprot_align", "all"]:
        ## make sure to reset infile here for uniprot_align and all
        #opts = {**mod_options["uniprot"], **universal_options}
        #print(f"Running uniprot_fastas.py with the following options: {opts}", flush=True)
        #uf.main(opts)
        
        ## read all resulting folders in current out_directory; and save full file paths as a list; these will be the new infiles for alignment
        ## set equal to list of infiles for next thing....maybe...need to fix
        ## os.walk() is a generator
        ## next() runs it just once, retrieving the subdirectories (for each protein) created by uniprot_fastas.py
        ## The output is a list: ["DIRECTORY", [SUBDIR1, SUBDIR2, ..., SUBDIRN], [FILE1, FILE2, ..., FILEN]]
        ## Note that this doesn't give full subdirectory paths; just names
        #file_list = [] # reset if uniprot_fastas.py is run
        #for d in next(os.walk(out_directory))[1]:
            #file_list.append((f"{out_directory}/{d}/all.fasta", f"{out_directory}/{d}/")) # the .fasta file of seqs will be stored in this manner
    #if module in ["align", "uniprot_align", "align_svg", "all"]:
        ## make sure to reset infile here for align_svg and all
        #for infile, out_directory in file_list:
            #universal_options["infile"] = infile
            #universal_options["out_directory"] = out_directory
            #opts = {**mod_options["align"], **universal_options}
            #print(f"Running alignment.py with the following options: {opts}", flush=True)
            #alignment.main(opts)
            
        #for i, ftup in enumerate(file_list):
            ##prot_dir = "/".join(ftup[1].split("/")[:-1]) # use the previous_out_directory
            #title = opts["title"]
            #file_list[i] = (f"{out_directory}/{title}.clustal_num", ftup[1])
    
    #if module in ["clust_to_svg", "align_svg", "all"]:
        #for infile, out_directory in file_list:
            ## set annotation files if exist/if they were generated in this run
            #if mod_options.get("clust_to_svg").get("annotations"):
                #mod_options["clust_to_svg"]["annotations"] = find_path(mod_options["clust_to_svg"]["annotations"], "r")
            #elif module in ["align_svg", "all"]:
                #prot_dir = "/".join(infile.split("/")[:-1])
                #mod_options["clust_to_svg"]["annotations"] = find_path(f"{prot_dir}/all.ann", "r")
            #else:
                #print("Proceeding without annotations.\n", flush=True)
                
            #universal_options["infile"] = infile
            #universal_options["out_directory"] = out_directory + "/svg_alignments/" # hardcoded...consider allowing options later
            #opts = {**mod_options["clust_to_svg"], **universal_options}
            #print(f"Running clustal_to_svg.py with the following options: {opts}", flush=True)
            #cts.main(opts)
            
            
### split all of this into simpler parts; make extra script to run things all at once maybe
### add retrieve_annotations.py
if __name__ == "__main__":
    main()