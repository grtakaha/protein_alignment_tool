"""
Converts a given Clustal alignment to a reformatted Inkscape SVG.

Extra imports:
    pandas
        https://pandas.pydata.org/docs/getting_started/install.html
    helpers
        https://github.com/grtakaha/protein_alignment_tool/blob/main/helpers.py

Functions:
    get_max_header(protein_classes.Alignment, str) -> int
    format_alignment(protein_classes.Alignment, str, str) -> list
    get_conserved_coords(protein_classes.Alignment, int, str) -> dict
    get_feature_coords(protein_classes.Alignment, int, str) -> dict
    add_conservation(str, int, int) -> str
    create_feature(protein_classes.Feature, int, int) -> str
    create_tspan(str, str, None) -> str
    split_lines(list, int, int) -> list
    create_svg(protein_classes.Alignment, str, str, str, str)
    parse_args() -> Namespace
    main()

Command-Line Arguments:
    --infile
    --out_directory
    --codes
    --nums
    --uniprot_format
    --annotations
"""

import argparse
import pandas as pd
from helpers import find_path, read_alignment

def get_max_header(alignment, nums="FALSE"):
    """
    Returns the true max header length, given the alignment's maximum and nums status.

        Parameters:
            alignment (protein_classes.Alignment): Alignment object with features.
            nums (str): If "TRUE", incorporates residue numbers into calculations.

        Returns:
            max_header (int): Maximum length of a protein display name.
    """

    # Stored as the default or manually set max header length.
    max_header = alignment.max_header

    # Accounts for space + numbers up to max digits.
    # len(str(alignment.get_max_length())) = length of the maximum integer residue number.
    # Ex. alignment.max_length = 354; max_header -= 2 + 3
    if nums == "TRUE":
        max_header -= 2 + len(str(alignment.get_max_length()))

    return max_header

def format_alignment(alignment, codes="FALSE", nums="FALSE"):
    """
    Returns a list of lines for an alignment, reformatted for Martin Lab figures.

        Parameters:
            alignment (protein_classes.Alignment): Alignment object with features.
            codes (str): If "TRUE", includes Clustal identity codes at the bottom of each block.
            nums (str): If "TRUE", includes total residue numbers at the end of each line.

        Returns:
            lines (list): List of lines for an alignment, reformatted for Martin lab figures.
    """

    # TODO: Consider storing this formatted alignment, line by line, in an alignment object

    lines = []
    max_header = get_max_header(alignment, nums=nums)
    counts = {} # count of amino acids for each protein

    # Base sequence length around this first Protein object.
    prot1 = list(alignment.proteins.values())[0]
    for i in range(0, len(prot1.sequence), 100):
        for name in alignment.proteins: # gives names only
            lines.append("") # start of a new line
            prot = alignment.proteins[name]

            # Add display name and spaces up to max header length + 2.
            lines[-1] += prot.disp_name[:max_header] + (" "*(max_header-len(prot.disp_name))) + "  "
            seq_100 = prot.sequence[i:i+100]
            for i_10 in range(0, len(seq_100), 10): # add spaces between sequence
                lines[-1] += seq_100[i_10:i_10+10] + " "
            lines[-1] = lines[-1][:-1] # remove last " "
            if nums == "TRUE":
                if counts.get(name) is None:
                    counts[name] = 0
                counts[name] += len(seq_100) - seq_100.count("-")
                lines[-1] += f"  {str(counts[name])}"
        lines.append("") # start of new line
        if codes == "TRUE":
            lines[-1] += " "*(max_header+2) # start of code line
            codes_100 = alignment.codes[i:i+100]
            for i_10 in range(0, len(codes_100), 10):
                lines[-1] += codes_100[i_10:i_10+10] + " "
        lines[-1] = lines[-1][:-1] # remove last " "
        lines.append("") # start of new line between blocks
    return lines

def get_conserved_coords(alignment, lines_per_svg=72, nums="FALSE"):
    """
    Returns a dictionary of SVG numbers and their respective conserved residue coordinates.

        Parameters:
            alignment (protein_classes.Alignment): Alignment object with features.
            lines_per_svg (int): Max number of lines allowed per SVG.
            nums (str): If "TRUE", includes residue numbers in x_coord calculations.

        Returns:
            conserved_dict (dict): Dictionary of SVG#:conserved_residue_list pairs.
                                   Ex. {0:[(residue, center_x, center_y, i, char_num)]}
    """

    x_start, y_start = 9.11, 20.274 # in px
    char_width, line_height = 5.5, 11
    width_gap = (701.576 - 127*char_width) / (127-1) # in px
    height_gap = (827.5 - 72*line_height) / (72-1) # in px

    lines_per_block = len(alignment.proteins) + 2

    # int() rounds down, so this will be <= 72.
    max_lines = int(lines_per_svg/lines_per_block) * lines_per_block
    max_blocks = int(max_lines / lines_per_block)

    # alignment.conserved_res was set by alignment.set_conserved_res() in main().
    conserved_dict = {} # {svg_num:[(f1, x, y), (f2, x, y), etc.], svg_num:[(f3, x, y), etc.], etc.}
    for res, i in alignment.conserved_res:
        # Note: i, the overall position, is 0-indexed in this.
        # This is different from feature.clust_start, which is 1-indexed.

        # I realize this is all very repetitive...it was the only way I could think about it.
        tot_block_num = int(i / 100) + min(1, i % 100) # 1-indexed

        # Necessary for svg #
        svg_num = int(tot_block_num / max_blocks - 1 + min(1, tot_block_num % max_blocks))


        block_num = int(tot_block_num % max_blocks)
        if block_num == 0:
            block_num = max_blocks

        char_num = int(i % 100) # Necessary for x
        if char_num == 0:
            char_num = 100

        # Necessary for x
        # Remove one character, then add it back if this is char_num % 10 is not 0.
        char_num += (int(char_num / 10) - 1 + min(1, char_num % 10) + 2 +
                     get_max_header(alignment, nums=nums))

        #line_num = lines_per_block * (block_num - 1) # necessary for y

        # center_x is the center of the circle.
        # 127 is the number of characters that can fit on each line.
        center_x = x_start + ((char_width / 2) * ((2 * char_num) - 1)) + (width_gap * char_num)

        # center_y is the center of the circle.
        # 72 is the number of lines that can fit in this height.
        center_y = (y_start - line_height/2 - height_gap +
                    ((line_height + height_gap) * lines_per_block*(block_num-1)))

        if not conserved_dict.get(svg_num):
            conserved_dict[svg_num] = []
        # Going to be honest, I can't remember why I returned i and char_num...
        conserved_dict[svg_num].append((res, center_x, center_y, i, char_num))

    return conserved_dict

def get_feature_coords(alignment, lines_per_svg=72, nums="FALSE"):
    """
    Returns a dictionary of SVG numbers and their respective feature coordinates.

        Parameters:
            alignment (protein_classes.Alignment): Alignment object with features.
            lines_per_svg (int): Max number of lines allowed per SVG.
            nums (str): If "TRUE", includes residue numbers in x_coord calculations.

        Returns:
            feature_dict (dict): Dictionary of SVG#:feature_list pairs.
                                 Ex. {0:[(f1, x, y, prot_name1)]}
    """

    # TODO: Consider redoing x-coord. Rectangle widths are still not exactly aligned.
    x_start, y_start = 9.11, 20.274 # In px
    char_width, line_height = 5.280, 11
    width_gap = (701.576 - 127*char_width) / (125 + 2) # In px

    height_gap = (827.5 - 72*line_height) / (72-1) # In px

    lines_per_block = len(alignment.proteins) + 2

    # int() rounds down, so this will be <= 72.
    max_lines = int(lines_per_svg/lines_per_block) * lines_per_block
    max_blocks = int(max_lines / lines_per_block)

    # Ex. {svg_num1:[(f1, x, y, prot_name1)], svg_num2:[(f3, x, y, prot_name3)]}
    feature_dict = {}
    # i is the number of this protein in an "ordered" dictionary.
    for i, name in enumerate(alignment.proteins):
        protein = alignment.proteins[name]

        # I realize this is all very repetitive...it was the only way I could think about it.
        for feature in protein.features:
            tot_block_num = int(i / 100) + min(1, i % 100) # 1-indexed

            # Necessary for svg #
            svg_num = int(tot_block_num / max_blocks - 1 + min(1, tot_block_num % max_blocks))

            tot_block_num = int(feature.clust_start / 100) + 1 # 1-indexed

            # Necessary for svg #
            svg_num = int(tot_block_num / max_blocks - 1 + min(1, tot_block_num % max_blocks))
            block_num = int(tot_block_num / (svg_num + 1))
            char_num = int(feature.clust_start % 100) # Necessary for x

            # Necessary for x
            # Remove one character, then add it back if this is char_num % 10 is not 0.
            # get_max_header() takes nums status into account.
            char_num += (int(char_num / 10) - 1 + min(1, char_num % 10) + 2 +
                         get_max_header(alignment, nums=nums))

            line_num = lines_per_block * (block_num - 1) + (i+1) # Necessary for y

            x_coord = x_start + width_gap + (char_num - 1)*(char_width + width_gap)
            y_coord = y_start + (line_num - 1)*(line_height + height_gap)

            if not feature_dict.get(svg_num):
                feature_dict[svg_num] = []
            # feature is a protein_classes.Feature object
            # x_coord, y_coord are integers
            # protein.disp_name is a string
            feature_dict[svg_num].append((feature, x_coord, y_coord, protein.disp_name))

    return feature_dict

def add_conservation(res, center_x, center_y):
    """
    Creates a transparent rectangle to annotate the given feature.

        Parameters:
            res (str): Conserved residue, used for color coding.
            center_x (int): x-coordinate for the center of an Inkscape circle.
            center_y (int): y-coordinate for the center of an Inkscape circle.

        Returns:
            circle (str): Circle line to be added to the final SVG.
    """

    color_pos = "3953a4"
    color_neg = "ed1c24"
    color_special = "f7ec13"
    color_hphobe = "3bb54a"
    color_others = "000000"

    color_codes = {"R":color_pos, "H":color_others, "K":color_pos,
                   "D":color_neg, "E":color_neg,
                   "S":color_others, "T":color_others, "N":color_others, "Q":color_others,
                   "C":color_special, "G":color_others, "P":color_others,
                   "A":color_hphobe, "V":color_hphobe, "I":color_hphobe,
                   "L":color_hphobe, "M":color_hphobe, "F":color_hphobe,
                   "Y":color_others, "W":color_hphobe, "default": "000000"}

    circle = f"<circle cx=\"{center_x}\" cy=\"{center_y}\" r=\"3px\" style=\"fill:#"
    circle += f"{color_codes.get(res)};stroke-width:0;fill-opacity:1;stroke:rgb(0,0,0)\" />"

    return circle

def create_feature(feature, left_x, top_y):
    """
    Creates a transparent rectangle to annotate the given feature.

        Parameters:
            feature (protein_classes.Feature): Feature object to be anntoated.
                                               Not implemented.
            left_x (int): x-coordinate for the left side of an Inkscape rectangle.
            top_y (int): y-coordinate for the top of an Inkscape rectangle.

        Returns:
            rect (str): Rectangle line to be added to the final SVG.
    """

    # TODO: Add in functionality for other features.
    # Create the following for a transparent rectangle.
    # Add it to the svg based on x and y of characters.
    rectangle = f"<rect x=\"{left_x}\" y=\"{top_y}\" width=\"5.5px\" height=\"11px\" "
    rectangle += "style=\"fill:rgb(0,0,255);stroke-width:0;fill-opacity:.2;stroke:rgb(0,0,0)\" />"

    return rectangle

def create_tspan(text, tspan_id, feature=None):
    """
    Creates a tspan SVG line with correct tspan ID, text, and fill.

        Parameters:
            text (str): Text to be added to the SVG tspan.
            tspan_id (str): ID used to identify protein, line, and residue numbers.
            feature (None): If set, will override default coloring with specified feature
                            coloring. Not implemented.

        Returns:
            tspan (str): tspan line to be added to the final SVG.
    """

    # TDOO: Consider making this editable for things like active sites
    # TODO: Maybe give proteins something to override this
    # TODO: Add feature markup (for things like active His) that will override default colors

    color_pos = "3953a4"
    color_neg = "ed1c24"
    color_special = "f7ec13"
    color_hphobe = "3bb54a"
    color_others = "000000"

    color_codes = {"R":color_pos, "H":color_others, "K":color_pos,
                   "D":color_neg, "E":color_neg,
                   "S":color_others, "T":color_others, "N":color_others, "Q":color_others,
                   "C":color_special, "G":color_others, "P":color_others,
                   "A":color_hphobe, "V":color_hphobe, "I":color_hphobe,
                   "L":color_hphobe, "M":color_hphobe, "F":color_hphobe,
                   "Y":color_others, "W":color_hphobe, "default": "000000"}

    code = color_codes.get(text)
    if code is None:
        code = color_codes["default"]

    return f"<tspan id=\"{tspan_id}\" style=\"fill:#{code};\">{text}</tspan>"

def split_lines(lines, lines_per_block, lines_per_svg=72):
    """
    Splits text into groups of 72 (with appropriate breaks - only between blocks).

        Parameters:
            lines (list):           List of lines in an alignment.
                                    Formatted by format_alignment().
            lines_per_block (int):  Total number of proteins in alignment + 2.
            lines_per_svg (int):    Max number of lines allowed per SVG.

        Returns:
            svg_list (list): List of lists of lines that go in each SVG.
    """

    # int() rounds down, so max_lines will be <= 72.
    max_lines = int(lines_per_svg/lines_per_block) * lines_per_block

    # List of lists of lines that go in each SVG.
    svg_list = []
    for i in range(0, len(lines), max_lines):
        svg_list.append(lines[i:i+max_lines])

    return svg_list

# TODO: Fix bug with identical display names - ID#s will add for both.
# TODO: Add numbers to "-" and " " so that they have real IDs.
def create_svg(alignment, out_directory, codes="FALSE", nums="FALSE", features="TRUE"):
    """
    Creates reformatted alignment Inkscape SVGs from a given Clustal alignment.

        Parameters:
            alignment (protein_classes.Alignment): Alignment object with features.
            out_directory (str): Directory for storage of reformatted alignments.
            codes (str): If "TRUE", includes Clustal identity codes at the bottom of each block.
            nums (str): If "TRUE", includes total residue numbers at the end of each line.
            features (str): If "TRUE", includes feature annotations.

        Outputs:
            A series of Inkscape SVGs with Martin Lab-formatted alignments.
            SVGs are named 0.svg, 1.svg, etc.
    """

    # Adjust as needed; good enough for now.
    line_height = 1.25
    # In px; I have no idea why this doesn't work with the rectangle y...figure out later
    x_start, y_start = 9.317, 28.5

    # Format lines
    text = format_alignment(alignment, codes=codes, nums=nums) # Returns list of lines.

    lines_per_block = len(alignment.proteins) + 2
    svgs = split_lines(text, lines_per_block) # List of svgs with up to 72 lines each.

    max_header = get_max_header(alignment, nums=nums)
    feature_coords = get_feature_coords(alignment, nums=nums)
    conserved_coords = get_conserved_coords(alignment, nums=nums)

    relevant_features = ["Active site"] # TODO: Add more valid features.

    aa_counts = {} # Use display names and keep track of aa#.
    line_counts = {} # Use display names and keep track of line#.

    aa_list = ["R", "H", "K", "D", "E", "S", "T",
               "N", "Q", "C", "G", "P", "A", "V",
               "I", "L", "M", "F", "Y", "W"]

    for i, lines in enumerate(svgs):

        # Inkscape SVG things.
        svg = "<?xml version=\"1.0\" encoding=\"utf-8\" ?>\n"
        svg += "<svg baseProfile=\"full\" "
        svg += "height=\"9in\" version=\"1.1\" width=\"7.5in\" "
        svg += "xmlns=\"http://www.w3.org/2000/svg\" "
        svg += "xmlns:sodipodi=\"http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd\" "
        svg += "xmlns:ev=\"http://www.w3.org/2001/xml-events\" "
        svg += "xmlns:xlink=\"http://www.w3.org/1999/xlink\"><defs />"
        svg += f"<text style=\"fill:#000000;line-height:{line_height};"
        svg += "font-family:Courier New;font-size:9.2px;font-weight:bold\" "
        svg += f"xml:space=\"preserve\" x=\"{x_start}\" y=\"{y_start}\">"

        # TODO: Consider removing need to re-split.
        for line in lines:
            # Important to add this back in, or they will not act like lines.
            line += "\n"
            if line == "\n":
                # Add a space if empty so user has something to "see".
                header = " "
                current_prot = "line"
            else:
                # Protein name and spaces.
                header = line[:max_header+2]
                current_prot = header.split(" ")[0]
                if current_prot not in aa_counts:
                    # Start at 0 and increment when valid AA is found.
                    aa_counts[current_prot] = 0

            if current_prot not in line_counts:
                line_counts[current_prot] = 0
            line_counts[current_prot] += 1

            svg += f"<tspan id=\"{current_prot}{line_counts[current_prot]}\" "
            svg += f"sodipodi:role=\"line\"><tspan>{header}</tspan>"
            for residue in line[max_header+2:]:
                if residue in aa_list:
                    aa_counts[current_prot] += 1
                    res_id = f"{current_prot}:{residue}{aa_counts[current_prot]}"
                else:
                    res_id = f"{residue}"
                svg += create_tspan(residue, res_id)
            svg += "</tspan>"
        svg += "</text>"

        # Add in features (transparent rectangles).
        if features == "TRUE" and feature_coords.get(i):
            for f_info in feature_coords[i]:
                feature = f_info[0]
                if feature.feature_type in relevant_features:
                    print(f"{feature.feature_type} found for " +
                          f"{f_info[3]}, residue {feature.start}, " +
                          f"alignment position {feature.clust_start}. " +
                          f"Adding to {i}.svg.", flush=True)
                    svg += create_feature(f_info[0], f_info[1], f_info[2]) # feature, x, y

        # Add in conserved residues (circles).
        if conserved_coords.get(i):
            for cons_info in conserved_coords[i]:
                # residue, center_x, center_y
                svg += add_conservation(cons_info[0], cons_info[1], cons_info[2])

        svg += "</svg>"

        out = f"{out_directory}{i}.svg".replace("\\", "/")
        with open(out, "w", encoding="utf-8") as out_svg:
            print(f"Writing part {i} to {out}.\n", flush=True)
            out_svg.write(svg)

    # TODO: Add non-active site features.

def parse_args():
    """
    Takes in command-line arguments and returns an argparse Namespace object.

        Returns:
            arguments (Namespace): Namespace with command-line arguments.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--out_directory", default="./",
                        help="Full path of output directory. Must end with \"/\".")
    parser.add_argument("-c", "--codes", default="FALSE",
                        help="When set to TRUE, includes Clustal identity " +
                        "codes at the bottom of each block.")
    parser.add_argument("-n", "--nums", default="FALSE",
                        help="When set to TRUE, includes total residue numbers " +
                        "at the end of each line.")
    parser.add_argument("-u", "--uniprot_format", default="FALSE",
                        help="When set to TRUE, truncates all accessions as if " +
                        "they were UniProt entries.\n" +
                        "Ex. sp|P00784|PAPA1_CARPA -> PAPA1_CARPA")
    parser.add_argument("-a", "--annotations", default="",
                        help="If an annotation file is provided, it will be " +
                        "used to annotate the resulting SVG files.")

    return parser.parse_args()

def main():
    """
    Converts a given Clustal alignment to a reformatted Inkscape SVG.

        Outputs:
            A series of Inkscape SVGs with Martin Lab-formatted alignments.
            SVGs are named 0.svg, 1.svg, etc.
    """

    args = parse_args()

    infile = find_path(args.infile, action="r").replace("\\", "/")
    print(f"Processing sequences from {infile} \n", flush=True)

    out_directory = find_path(args.out_directory, action="w").replace("\\", "/")
    print(f"Storing outputs in {out_directory}\n", flush=True)

    alignment = read_alignment(infile)
    alignment.set_disp_names(uniprot_format=args.uniprot_format)

    annotations = args.annotations
    if annotations != "":
        annotations = pd.read_csv(find_path(annotations, "r"), sep="\t")
        features = "TRUE"
    else:
        annotations = pd.DataFrame()
        features = "FALSE"
    alignment.add_features(annotations) # Add nothing if annotations aren't set.

    create_svg(alignment, find_path(out_directory, "w"),
               codes=args.codes, nums=args.nums, features=features)

if __name__ == "__main__":
    main()
