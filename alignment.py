"""
Retrieves annotations for given UniProt IDs.

Extra imports:
    api_funcs
        https://github.com/grtakaha/protein_alignment_tool/blob/main/api_funcs.py
    helpers
        https://github.com/grtakaha/protein_alignment_tool/blob/main/helpers.py

Functions:
    parse_args() -> Namespace
    main()

Command-Line Arguments:
    --infile
    --out_directory
    --stype
    --email
    --title
"""

import argparse
import api_funcs as af
from helpers import find_path

def parse_args():
    """
    Takes in command-line arguments and returns an argparse Namespace object.

        Returns:
            arguments (Namespace): Namespace with command-line arguments.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--infile", help="Full path of input file.")
    parser.add_argument("-o", "--out_directory", default="./",
                        help="Full path of output directory. Must end with \"/\".")
    parser.add_argument("-s", "--stype", default="protein",
                        help="Sequence type (\"protein\" or \"dna\").")
    parser.add_argument("-e", "--email",
                        help="Personal email. Used to submit BLAST and Clustal Omega jobs.")
    parser.add_argument("-t", "--title", default="alignment",
                        help="Output title ([TITLE].clustal_num, [TITLE].pim).")

    return parser.parse_args()

def main():
    """
    Aligns an input FASTA file with Clustal Omega.

        Outputs:
            An alignment (.clustal_num) and percent identity matrix (.pim)
            created from the input FASTA file.
    """

    args = parse_args()

    infile = find_path(args.infile, action="r").replace("\\", "/")
    print(f"Processing sequences from {infile} \n", flush=True)
    with open(infile, "r", encoding="utf-8") as f_temp:
        seqs = "".join(f_temp.readlines())

    out_directory = find_path(args.out_directory, action="w").replace("\\", "/")
    print(f"Storing outputs in {out_directory}\n", flush=True)

    stype = args.stype
    if stype not in ["dna", "protein"]:
        print("Given stype not found. Please specify " +
              "\"-stype dna\" OR \"-stype protein\n", flush=True)
        print("Defaulting to stype=\"protein\"\n", flush=True)
        stype = "protein"

    title = args.title

    aid = af.align(args.email, stype, title, seqs)
    alignment = af.get_alignment(aid)
    pim = af.get_pim(aid)

    alignment_file = f"{out_directory}/{title}.clustal_num"
    pim_file = f"{out_directory}/{title}.pim"

    print(f"Saving {title} alignment to {alignment_file}.", flush=True)
    with open(alignment_file, "w", encoding="utf-8") as f_temp:
        f_temp.write(alignment) # just overwriting it if it exists
    with open(pim_file, "w", encoding="utf-8") as p_temp:
        p_temp.write(pim) # just overwriting if it exists

if __name__ == "__main__":
    main()
