import argparse
import requests
import os
import time
import pandas as pd
import json
import alignment
import clustal_to_svg as cts
import subprocess
from helpers import find_path, fasta_to_df

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
    multi.add_argument("-ord", "--order", nargs="+") # will return a list of 1 or more tools to use with these arguments
    multi.add_argument("-s", "--stype", default="protein") # program is tied to stype (dna:blastx, rna:blastx, protein:blastp)
    multi.add_argument("-e", "--email", default="")
    multi.add_argument("-nr", "--num_res", default="10")
    multi.add_argument("-t", "--title", default="alignment")
    multi.add_argument("-c", "--codes", default="FALSE")
    multi.add_argument("-n", "--nums", default="FALSE")
    multi.add_argument("-u", "--uniprot_format", default="FALSE")
    multi.add_argument("-a", "--annotations", default="") # will annotate if this is provided at all #### cannot be nonetype
 
    return parser.parse_args()

##### TODO: MAKE IT SO THAT FILES CAN BE FOUND ANYWHERE ##### DONE
##### TODO: MAKE SURE NON-UNIPROT QUERIES DON'T BREAK ANNOTATION #####
##### TODO: MAKE IT SO THAT ANNOTATIONS AREN'T REQUIRED ##### DONE
def execute_tool(args):
    # return output files? or have a separate function to return output files for a particular tool
    # Execute the selected tool
    print(f"Running {args.tool_name}...\n", flush=True)
    if args.tool_name == "annotate": # annotation command main_tool.py ann ...
        subprocess.run(["python", f"{os.path.abspath(os.path.dirname(__file__))}/retrieve_annotations.py",
                        "--infile", args.infile,
                        "--out_directory", args.out_directory])
        args.annotations = find_outputs(args)[0] # will change args

    elif args.tool_name == "blast": # replace with split uniprot_fastas into blast and split blast...
        subprocess.run(["python", f"{os.path.abspath(os.path.dirname(__file__))}/search_proteins.py", 
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
            print(f"In multi run, executing {tool}...\n", flush=True)
            args.tool_name = tool
            #print(f"POS1: {args}")
            execute_tool(args)
            #print(f"POS2: {args}")
            if tool == "blast" and i+1 < len(args.order): # will only continue if there is something after blast
                # run everything from here on out using multi and different inputs
                infiles = find_outputs(args)
                print(f"blast tool executed in a multi run.")
                print(f"Remaining calls ({args.order}) will be executed on {infiles}.\n")
                ##### MAY ERROR OUT BECAUSE +1 MIGHT NOT BE IN RANGE #####
                args.order = args.order[i+1:] # start after blast call
                for f in infiles:
                    #print(f"ORDER: {args.order}")
                    #print(f"INFILE: {f}")
                    args.tool_name = "multi"
                    args.infile = f
                    args.out_directory = "/".join(f.split("/")[:-1]) + "/"
                    execute_tool(args)
                break # BREAKS OUT OF OUTER LOOP IF BLAST IS FOUND...PREVENTS EXTRA EXECUTION OF ANNOTATE, ALIGN, AND SVG AFTER BLAST
                
# instead of using this, consider returning outputs, then printing them separately only if this is run from main/main_tool.py
def find_outputs(args):
    new_inputs = []
    if args.tool_name == "annotate":
        # it has another output for individual annotations
        new_inputs.append(f"{args.out_directory}/all.ann".replace("//", "/")) # tsv for annotations input later
    elif args.tool_name == "blast":
        proteins = fasta_to_df(find_path(args.infile, "r"))
        for p in proteins.index.values:
            new_inputs.append(f"{args.out_directory}/{p}/all.fasta".replace("//", "/")) # fastas of blast hits
    elif args.tool_name == "align":
        new_inputs.append(f"{args.out_directory}/{args.title}.clustal_num".replace("//", "/"))
    elif args.tool_name == "svg":
        pass # not necessary to return anything because svg is the end of the line
    return new_inputs

def main():
    args = parse_args()
    #print(args)
    execute_tool(args)
    # add execution only after args works

if __name__ == "__main__":
    main()