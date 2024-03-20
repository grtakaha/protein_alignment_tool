import protein_classes as pc
import argparse
import os
import pandas as pd

# this is the redesigned clustal_to_svg tool

def get_max_header(alignment, nums="FALSE"):
    mh = alignment.max_header
    # accounts for space + numbers up to max digits; e.g. max length = 354; mh -= 2 + 3
    if nums == "TRUE":
        mh -= 2 + len(str(alignment.get_max_length()))
    return mh

def format_alignment(alignment, codes="FALSE", nums="FALSE"):
    ##### consider storing this formatted alignment, line by line, in an alignment object #####
    ##### add a features list with feature objects and coordinates for each svg #####
    
    aa_count = 0
    lines = []
    mh = get_max_header(alignment, nums=nums)
    counts = {} # count of amino acids for each protein
    prot1 = list(alignment.proteins.values())[0] # base sequence length around this first Protein object
    for i in range(0, len(prot1.sequence), 100):
        for name in alignment.proteins: # gives names only
            lines.append("") # start of a new line
            prot = alignment.proteins[name]
            lines[-1] += prot.disp_name[:mh] + (" "*(mh-len(prot.disp_name))) + "  " # add display name and spaces up to max header length + 2
            s = prot.sequence[i:i+100]
            for a in range(0, len(s), 10): # add spaces between sequence
                lines[-1] += s[a:a+10] + " "
            lines[-1] = lines[-1][:-1] # remove last " "
            if nums == "TRUE":
                if counts.get(name) == None:
                    counts[name] = 0
                counts[name] += len(s) - s.count("-")
                lines[-1] += f"  {str(counts[name])}"
        lines.append("") # start of new line
        if codes == "TRUE":
            lines[-1] += " "*(mh+2) # start of code line
            c = alignment.codes[i:i+100]
            for a in range(0, len(c), 10):
                lines[-1] += c[a:a+10] + " "
        lines[-1] = lines[-1][:-1] # remove last " "
        lines.append("") # start of new line between blocks
    return lines

def get_conserved_coords(alignment, lines_per_svg=72, nums="FALSE"):
    ##### return list of tuples with x and y coordinates for each svg (given the max number of lines)
    x, y = 9.11, 20.274 # in px
    char_width, line_height = 5.5, 11
    width_gap = (701.576 - 127*char_width) / (127-1) # in px
    height_gap = (827.5 - 72*line_height) / (72-1) # in px
      
    lines_per_block = len(alignment.proteins) + 2
    print(lines_per_block)
    max_lines = int(lines_per_svg/lines_per_block) * lines_per_block # int rounds down, so this will be <= 72
    max_blocks = int(max_lines / lines_per_block)
    print(max_blocks)
    
    conserved_dict = {} # {svg_num:[(f1, x, y), (f2, x, y), etc.], svg_num:[(f3, x, y), etc.], etc.}
    for c, i in alignment.conserved_res:
        # note: i, the clustal position, is 0-indexed in this; this is different from feature.clust_start, which I believe is 1-indexed
        # I realize this is all very repetitive...it was the only way I could think about it
        tot_block_num = int(i / 100) + min(1, i % 100) # 1-indexed

        svg_num = int(tot_block_num / max_blocks - 1 + min(1, tot_block_num % max_blocks)) # necessary for svg #

        block_num = int(tot_block_num % max_blocks)
        if block_num == 0:
            block_num = max_blocks
        
        char_num = int(i % 100) # necessary for x
        if char_num == 0:
            char_num = 100 
        char_num += int(char_num / 10) - 1 + min(1, char_num % 10) + 2 + get_max_header(alignment, nums=nums) # necessary for x; remove one character, then add it back if this is char_num % 10 is not 0
        
        line_num = lines_per_block * (block_num - 1) # necessary for y
        
        
        cx = x + ((char_width / 2) * ((2 * char_num) - 1)) + (width_gap * char_num)  # cx is the center of the circle; 127 is the number of characters that can fit on each line
        cy = y - line_height/2 - height_gap + ((line_height + height_gap)*lines_per_block*(block_num-1)) # cy is the center of the circle; 72 is the number of lines that can fit in this height
        
        if not conserved_dict.get(svg_num):
            conserved_dict[svg_num] = []
        conserved_dict[svg_num].append((c, cx, cy, i, char_num))
        
    return conserved_dict

def get_feature_coords(alignment, lines_per_svg=72, nums="FALSE"):
    ##### return list of tuples with x and y coordinates for each svg (given the max number of lines)
    # COME BACK TO THIS LATER MAYBE? RECTANGLES ARE STILL NOT EXACTLY ALIGNED WITH CHARS IN X-COORD...BUT NOT URGENT OR NOTICEABLE
    x, y = 9.11, 20.274 # in px
    char_width, line_height = 5.280, 11
    width_gap = (701.576 - 127*char_width) / (125 + 2) # in px
    #width_gap = (700.866 - 127*char_width) / (127-1) # in px
    height_gap = (827.5 - 72*line_height) / (72-1) # in px 
      
    lines_per_block = len(alignment.proteins) + 2
    print(lines_per_block)
    max_lines = int(lines_per_svg/lines_per_block) * lines_per_block # int rounds down, so this will be <= 72
    max_blocks = int(max_lines / lines_per_block)
    print(max_blocks)
    
    feature_dict = {} # {svg_num:[(f1, x, y), (f2, x, y), etc.], svg_num:[(f3, x, y), etc.], etc.}
    for i, name in enumerate(alignment.proteins): # is the number of this protein in an "ordered" dictionary
        protein = alignment.proteins[name]
        for feature in protein.features:
            # I realize this is all very repetitive...it was the only way I could think about it
            
            # chcck to see if the below code works...because the original was found to be buggy for conserved_coords()
            # the below, commented-out code works for conserved_coords()
            
            tot_block_num = int(i / 100) + min(1, i % 100) # 1-indexed
    
            svg_num = int(tot_block_num / max_blocks - 1 + min(1, tot_block_num % max_blocks)) # necessary for svg #
    
            ##### THIS DOES NOT WORK #####
            #block_num = int(tot_block_num % max_blocks)
            #if block_num == 0:
                #block_num = max_blocks
            
            #char_num = int(i % 100) # necessary for x
            #if char_num == 0:
                #char_num = 100 
            #char_num += int(char_num / 10) - 1 + min(1, char_num % 10) + 2 + get_max_header(alignment, nums=nums) # necessary for x; remove one character, then add it back if this is char_num % 10 is not 0
            ##### THIS DOES NOT WORK #####
            
            tot_block_num = int(feature.clust_start / 100) + 1 # 1-indexed

            svg_num = int(tot_block_num / max_blocks - 1 + min(1, tot_block_num % max_blocks)) # necessary for svg #
            block_num = int(tot_block_num / (svg_num + 1))
            char_num = int(feature.clust_start % 100) # necessary for x
            char_num += int(char_num / 10) - 1 + min(1, char_num % 10) + 2 + get_max_header(alignment, nums=nums) # necessary for x; remove one character, then add it back if this is char_num % 10 is not 0
           
            line_num = lines_per_block * (block_num - 1) + (i+1) # necessary for y
            
            x_coord = x + width_gap + (char_num - 1)*(char_width + width_gap) # fix later; this doesn't totally align; do it like height instead
            #x_coord = x + (char_num - 1)*(char_width + width_gap)
            y_coord = y + (line_num - 1)*(line_height + height_gap)
            
            if not feature_dict.get(svg_num):
                feature_dict[svg_num] = []
            feature_dict[svg_num].append((feature, x_coord, y_coord, protein.disp_name))
            
    return feature_dict

def add_conservation(c, cx, cy):
    # will add a colored circle above conserved residues
    color_pos = "3953a4" # new color
    color_neg = "ed1c24"
    color_special = "f7ec13" # new color
    color_hphobe = "3bb54a" # new color
    color_others = "000000" # consider not using...remain black; technically has a different color in pdf
    
    color_codes = {"R":color_pos, "H":color_others, "K":color_pos,
                   "D":color_neg, "E":color_neg,
                   "S":color_others, "T":color_others, "N":color_others, "Q":color_others,
                   "C":color_special, "G":color_others, "P":color_others,
                   "A":color_hphobe, "V":color_hphobe, "I":color_hphobe,
                   "L":color_hphobe, "M":color_hphobe, "F":color_hphobe,
                   "Y":color_others, "W":color_hphobe, "default": "000000"}
    
    return f"<circle cx=\"{cx}\" cy=\"{cy}\" r=\"3px\" style=\"fill:#{color_codes.get(c)};stroke-width:0;fill-opacity:1;stroke:rgb(0,0,0)\" />"

def create_feature(feature, x, y):
    ##### create the following for a transparent rectangle, and add it to the svg based on x and y of characters #####
    return f"<rect x=\"{x}\" y=\"{y}\" width=\"5.5px\" height=\"11px\" style=\"fill:rgb(0,0,255);stroke-width:0;fill-opacity:.2;stroke:rgb(0,0,0)\" />"
    
def create_tspan(c, ID, feature=None):
    # later, consider making this editable for things like active sites - maybe give proteins something to override this
    ##### ADD FEATURE MARKUP THAT WILL OVERRIDE COLORS LATER #####
    ##### CONSIDER PUTTING NON-COLOR MARKUP IN SEPARATE FUNCTION ######
    
    color_pos = "3953a4" # new color
    color_neg = "ed1c24"
    color_special = "f7ec13" # new color
    color_hphobe = "3bb54a" # new color
    color_others = "000000" # consider not using...remain black; technically has a different color in pdf
    
    color_codes = {"R":color_pos, "H":color_others, "K":color_pos,
                   "D":color_neg, "E":color_neg,
                   "S":color_others, "T":color_others, "N":color_others, "Q":color_others,
                   "C":color_special, "G":color_others, "P":color_others,
                   "A":color_hphobe, "V":color_hphobe, "I":color_hphobe,
                   "L":color_hphobe, "M":color_hphobe, "F":color_hphobe,
                   "Y":color_others, "W":color_hphobe, "default": "000000"}
    
    code = color_codes.get(c)
    if code == None:
        code = color_codes["default"]

    return f"<tspan id=\"{ID}\" style=\"fill:#{code};\">{c}</tspan>"

def split_lines(lines, lines_per_block, lines_per_svg=72):
    ''' splits text into groups of 72 (with appropriate breaks - only between blocks)
    lines is formatted by format_alignment()
    lines_per_block is # of proteins + 2
    '''
    max_lines = int(lines_per_svg/lines_per_block) * lines_per_block # int rounds down, so this will be <= 72
    print(f"Max lines: {max_lines}.\n", flush=True)
    
    svg_list = [] #list of lines that go in each svg
    for i in range(0, len(lines), max_lines):
        svg_list.append(lines[i:i+max_lines])
        
    return svg_list # list of lists

def create_svg(alignment, out_directory, codes="FALSE", nums="FALSE", features="TRUE"):
    # alignment is an Alignment object; out is the outfile path
    
    line_height = 1.25 # Adjust as needed; not sure if 1.25 is a good idea yet, but it's good enough for now
    x, y = 9.317, 28.5 # in px; I have no idea why this doesn't match the rectangle y...figure out later

    # Format lines
    text = format_alignment(alignment, codes=codes, nums=nums) # returns list of lines

    lines_per_block = len(alignment.proteins) + 2 # consider not using this, and just checking/adding to the proper svg as we go...
    svgs = split_lines(text, lines_per_block) # list of svgs with up to 72 lines each

    mh = get_max_header(alignment, nums=nums)
    feature_coords = get_feature_coords(alignment, nums=nums)
    conserved_coords = get_conserved_coords(alignment, nums=nums)
    #print(conserved_coords)
    
    relevant_features = ["Active site"] # add more later
    
    aa_counts = {} # use display names and keep track of aa#
    line_counts = {} # use display names and keep track of line#
    clust_num = 0 # 0-indexed because it is not counting just aa... it is counting positions in the alignment
    aa_list = ["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W"]
    #print(svgs)
    for i, lines in enumerate(svgs):

        svg = f"<?xml version=\"1.0\" encoding=\"utf-8\" ?>\n<svg baseProfile=\"full\" height=\"9in\" version=\"1.1\" width=\"7.5in\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:sodipodi=\"http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd\" xmlns:ev=\"http://www.w3.org/2001/xml-events\" xmlns:xlink=\"http://www.w3.org/1999/xlink\"><defs />"
        svg += f"<text style=\"fill:#000000;line-height:{line_height};font-family:Courier New;font-size:9.2px;font-weight:bold\" xml:space=\"preserve\" x=\"{x}\" y=\"{y}\">"
        
        for line in lines:
            line += "\n" # important to add this back in, or they will not act like lines; consider removing need to re-split
            if line == "\n":
                header = " " # add a space if this is an empty line, just so that something is there for the user to "see"
                current_prot = "line"
            else:
                header = line[:mh+2] # protein name and spaces; consider making this just the protein name
                current_prot = header.split(" ")[0]
                if current_prot not in aa_counts:
                    aa_counts[current_prot] = 0 # start at 0 and increment when valid AA is found
                    
            if current_prot not in line_counts:
                line_counts[current_prot] = 0
            line_counts[current_prot] += 1
            
            ##### ADD ID #####
            svg += f"<tspan id=\"{current_prot}{line_counts[current_prot]}\" sodipodi:role=\"line\"><tspan>{header}</tspan>"
            for c in line[mh+2:]:
                if c in aa_list:
                    aa_counts[current_prot] += 1
                    ID = f"{current_prot}:{c}{aa_counts[current_prot]}"
                else:
                    ID = f"{c}"
                svg += create_tspan(c, ID)
            svg += f"</tspan>"
        svg += f"</text>"
        # add in features
        if features == "TRUE" and feature_coords.get(i):
            for f in feature_coords[i]:
                feature = f[0]
                if feature.feature_type in relevant_features:
                    print(f"{feature.feature_type} found at {f[3]}, residue {feature.start}, {feature.clust_start}. Adding to {i}.svg.", flush=True)
                    print(f"{f[1]}, {f[2]}\n", flush=True)
                    svg += create_feature(f[0], f[1], f[2]) # feature, x, y
                    
        # not saving svgs to 0...
        if conserved_coords.get(i):
            for c in conserved_coords[i]:
                svg += add_conservation(c[0], c[1], c[2]) # c, cx, cy
                #print(c)
                
        svg += f"</svg>"
        
        out = f"{out_directory}{i}.svg"
        with open(out, "w") as o:
            print(f"Writing part {i} to {out}.", flush=True)
            o.write(svg)

    ##### add active site annotation/other features
    ##### consider grouping like-formatted things; not necessary, but does save on space; make sure it still allows for individual character formatting
    ##### consider not hard-coding centering (don't set x, set center to page)
    ##### add other features back in (like input options)
    ##### add API access to UniProt/database annotations/tabular annotations
    
def read_alignment(infile, max_header=16):
    ''' max_header is the max length (before 2 spaces are added) of the protein name + spaces.
    example: max_header=10; protein name + spaces = 12
    If the given protein name exceeds this value, it will be truncated.
    '''
    alignment = pc.Alignment() # empty alignment
    protein_key = {} # does not contain protein objects; just names and aligned seqs
    codes = ""
    header_len = 1000 # set so that it just uses an empty line if there is no header_len set
    with open(infile, "r") as clust:
        header = clust.readline() # includes "\n" at the end
        for line in clust:
            if line == "\n":
                continue # skip newlines
            else:
                line = line.rstrip("\n") # remove newlines
                if line[0] != " ": # find lines with protein seqs
                    l = line.split(" ")
                    name = l[0]
                    header_len = len(name) + line.count(" ") # used for codes; depends on codes coming after a block
                    if name not in protein_key:
                        protein_key[name] = ""
                    protein_key[name] += l[-1].split("\t")[0] # sequence
                else:
                    codes += line[header_len:] # consider adding own code checker
    
    # populate alignment
    alignment.add_codes(codes)
    alignment.set_max_header(max_header) # will add 2 spaces regardless; max_header is max accession length
    for name in protein_key:
        alignment.add_protein(pc.Protein(name, protein_key[name]))
    
    alignment.set_conserved_res() # set conserved residues for the completed alignment
    return alignment

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

def parse_args():
    # add descriptions to arguments later
    parser = argparse.ArgumentParser()
    # change clust_path arguments later
    parser.add_argument("-i", "--infile")    
    parser.add_argument("-o", "--out_directory", default="./") # consider changing to out_directory
    parser.add_argument("-c", "--codes", default="FALSE")
    parser.add_argument("-n", "--nums", default="FALSE")
    parser.add_argument("-u", "--uniprot_format", default="FALSE")
    parser.add_argument("-a", "--annotations", default="") # will annotate if this is provided at all
    # consider making features separate
    
    return parser.parse_args()

def main():
    args = parse_args()
    # if running on its own, will use command-line arguments; otherwise, pass arguments for this function
    # options must be a dictionary of valid options for this module
    
    # need to change out all these options arguments
    
    alignment = read_alignment(find_path(args.infile, "r"))
    alignment.set_disp_names(uniprot_format=args.uniprot_format)
    
    annotations = args.annotations
    if annotations != "":
        annotations = pd.read_csv(find_path(annotations, "r"), sep="\t")
        features = "TRUE"
    else:
        annotations = pd.DataFrame()
        features = "FALSE"
    alignment.add_features(annotations) # add nothing if annotations aren't set
    
    create_svg(alignment, find_path(args.out_directory, "w") , codes=args.codes, nums=args.nums, features=features)

if __name__ == "__main__":
    main()

