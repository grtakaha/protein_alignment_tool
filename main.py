"""
Manager for protein alignment tools.

Extra imports:
    helpers
        https://github.com/grtakaha/protein_alignment_tool/blob/main/helpers.py

Functions:
    parse_args() -> Namespace
    execute_tool(Namespace)
    main()

Command-Line Arguments:
    --infile
    --out_directory
    tool_name
    --order
    --stype
    --email
    --num_res
    --title
    --codes
    --nums
    --uniprot_format
    --annotations
"""

import argparse
import os
import subprocess
from helpers import find_path, fasta_to_df

def parse_args():
    """
    Takes in command-line arguments and returns an argparse Namespace object.

        Returns:
            arguments (Namespace): Namespace with command-line arguments.
    """

    parser = argparse.ArgumentParser(description='Protein Alignment Tool Manager')
    # -i and -o must come before tool_name.
    # Ex. -i ./infile.fasta -o ./out_dir annotate

    parser.add_argument("-i", "--infile", type=str, help="Full path of input file")
    parser.add_argument("-o", "--out_directory", type=str, help="Full path of output directory")

    # must be after -i and -o (positional)
    subparsers = parser.add_subparsers(dest='tool_name', help='Tool to execute')

    # Each tool should have its own argument parsing.
    blast = subparsers.add_parser("blast")
    blast.add_argument("-s", "--stype", default="protein")
    blast.add_argument("-e", "--email", default="")
    blast.add_argument("-nr", "--num_res", default="10")

    annotation = subparsers.add_parser("annotate")

    alignment = subparsers.add_parser("align")
    alignment.add_argument("-s", "--stype", default="protein")
    alignment.add_argument("-e", "--email", default="")
    alignment.add_argument("-t", "--title", default="alignment")

    clust_to_svg = subparsers.add_parser("svg")
    clust_to_svg.add_argument("-c", "--codes", default="FALSE")
    clust_to_svg.add_argument("-n", "--nums", default="FALSE")
    clust_to_svg.add_argument("-u", "--uniprot_format", default="FALSE")
    # Will annotate if this is provided at all; cannot be nonetype.
    clust_to_svg.add_argument("-a", "--annotations", default="")

    multi = subparsers.add_parser("multi")
    # Will provide a list of 1 or more tools to use with these arguments
    multi.add_argument("-ord", "--order", nargs="+")
    multi.add_argument("-s", "--stype", default="protein")
    multi.add_argument("-e", "--email", default="")
    multi.add_argument("-nr", "--num_res", default="10")
    multi.add_argument("-t", "--title", default="alignment")
    multi.add_argument("-c", "--codes", default="FALSE")
    multi.add_argument("-n", "--nums", default="FALSE")
    multi.add_argument("-u", "--uniprot_format", default="FALSE")
    # Will annotate if this is provided at all; cannot be nonetype.
    multi.add_argument("-a", "--annotations", default="")

    return parser.parse_args()

# TODO: Make sure non-UniProt entries don't break annotation.
# TODO: Remove subprocesses in favor of just importing each module.
def execute_tool(args):
    """
    Executes the tool(s) specified in args.

        Parameters:
            args (Namespace): Namespace with command-line arguments.
    """

    # Execute the selected tool.
    print(f"Running {args.tool_name}...\n", flush=True)
    if args.tool_name == "annotate": # annotation command main_tool.py ann ...
        subprocess.run(["python",
                        f"{os.path.abspath(os.path.dirname(__file__))}/retrieve_annotations.py",
                        "--infile", args.infile,
                        "--out_directory", args.out_directory])
        args.annotations = find_outputs(args)[0] # will change args

    elif args.tool_name == "blast":
        subprocess.run(["python",
                        f"{os.path.abspath(os.path.dirname(__file__))}/search_proteins.py",
                        "--infile", args.infile,
                        "--out_directory", args.out_directory,
                        "--stype", args.stype,
                        "--email", args.email,
                        "--num_res", args.num_res])

    elif args.tool_name == "align":
        subprocess.run(["python",
                        f"{os.path.abspath(os.path.dirname(__file__))}/alignment.py",
                        "--infile", args.infile,
                        "--out_directory", args.out_directory,
                        "--stype", args.stype,
                        "--email", args.email,
                        "--title", args.title])
        args.infile = find_outputs(args)[0] # .clustal_num file

    elif args.tool_name == "svg":
        subprocess.run(["python",
                        f"{os.path.abspath(os.path.dirname(__file__))}/clustal_to_svg.py",
                        "--infile", args.infile,
                        "--out_directory", args.out_directory,
                        "--codes", args.codes,
                        "--nums", args.nums,
                        "--uniprot_format", args.uniprot_format,
                        "--annotations", args.annotations])

    elif args.tool_name == "multi":
        # List of one or more: blast annotate align svg.
        for i, tool in enumerate(args.order):
            print(f"In multi run, executing {tool}...\n", flush=True)
            args.tool_name = tool
            execute_tool(args)
            # Will only continue if there is something after blast
            if tool == "blast" and i+1 < len(args.order):
                # Run everything from here on out using multi and different inputs
                infiles = find_outputs(args)
                print("blast tool executed in a multi run.")
                args.order = args.order[i+1:] # Start after blast call.
                print(f"Remaining calls ({args.order}) will be executed on {infiles}.\n")
                for in_f in infiles:
                    args.tool_name = "multi"
                    args.infile = in_f
                    args.out_directory = "/".join(in_f.split("/")[:-1]) + "/"
                    execute_tool(args)
                # Prevents extra execution of annotate, align, or svg after blast.
                break

# TODO: Consider returning outputs, then printing them separately.
def find_outputs(args):
    """
    Returns output file paths for the tool specified in args.

        Parameters:
            args (Namespace): Namespace with command-line arguments.

        Returns:
            new_inputs (list): List of output file paths for the specified tool.
                               Can be used as infiles for the next tool in a multi run.
    """

    new_inputs = []
    if args.tool_name == "annotate":
        # It also outputs individual annotation files.
        # This is a TSV of all annotations for a particular FASTA.
        new_inputs.append(f"{args.out_directory}/all.ann".replace("//", "/"))
    elif args.tool_name == "blast":
        proteins = fasta_to_df(find_path(args.infile, "r"))
        for prot in proteins.index.values:
            # Fastas of blast hits.
            new_inputs.append(f"{args.out_directory}/{prot}/all.fasta".replace("//", "/"))
    elif args.tool_name == "align":
        # TODO: Either make clustalo output .clustal_num, or change this to .clustal.
        new_inputs.append(f"{args.out_directory}/{args.title}.clustal_num".replace("//", "/"))
    elif args.tool_name == "svg":
        pass # Not necessary to return anything because svg is the end of the line.
    return new_inputs

def main():
    """
    Executes the given protein alignment tools (blast, annotate, align, and svg).

        Outputs:
            Outputs depend on the tool(s) being executed.
    """

    args = parse_args()
    execute_tool(args)

if __name__ == "__main__":
    main()
