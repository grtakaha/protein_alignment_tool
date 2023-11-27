import protein_classes as pc
import argparse
import os

# this is the redesigned clustal_to_svg tool

def get_max_header(alignment, nums="FALSE"):
    mh = alignment.max_header
    # accounts for space + numbers up to max digits; e.g. max length = 354; mh -= 2 + 3
    if nums == "TRUE":
        mh -= 2 + len(str(alignment.get_max_length()))
    return mh

def format_alignment(alignment, codes="FALSE", nums="FALSE"):
    ##### consider storing this formatted alignment, line by line, in an alignment object #####
    
    aa_count = 0
    lines = []
    mh = get_max_header(alignment, nums=nums)
    counts = {} # count of amino acids for each protein
    prot1 = list(alignment.proteins.values())[0] # base everything around this first Protein object
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
            lines[-1] += " "*(mh+2) # start of code line; only +1 because one space is added regardless
            c = alignment.codes[i:i+100]
            for a in range(0, len(c), 10):
                lines[-1] += c[a:a+10] + " "
        lines[-1] = lines[-1][:-1] # remove last " "
        lines.append("") # start of new line between blocks
    return lines

def create_tspan(c, markup=None):
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

    return f"<tspan style=\"fill:#{code};\">{c}</tspan>"

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

def create_svg(alignment, out_prefix, codes="FALSE", nums="FALSE"):
    # alignment is an Alignment object; out is the outfile path
    
    line_height = 1.25 # Adjust as needed; not sure if 1.25 is a good idea yet, but it's good enough for now
    x, y = 9.317, 28.5

    # Format lines
    text = format_alignment(alignment, codes=codes, nums=nums) # returns list of lines

    lines_per_block = len(alignment.proteins) + 2 # consider not using this, and just checking/adding to the proper svg as we go...
    svgs = split_lines(text, lines_per_block) # list of svgs with up to 72 lines each

    mh = get_max_header(alignment, nums=nums)
    
    for i, lines in enumerate(svgs):

        svg = f"<?xml version=\"1.0\" encoding=\"utf-8\" ?>\n<svg baseProfile=\"full\" height=\"9in\" version=\"1.1\" width=\"7.5in\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:sodipodi=\"http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd\" xmlns:ev=\"http://www.w3.org/2001/xml-events\" xmlns:xlink=\"http://www.w3.org/1999/xlink\"><defs />"
        svg += f"<text style=\"fill:#000000;line-height:{line_height};font-family:Courier New;font-size:9.2px;font-weight:bold\" xml:space=\"preserve\" x=\"{x}\" y=\"{y}\">"
        
        for line in lines:
            line += "\n" # important to add this back in, or they will not act like lines; consider removing need to re-split
            if line == "\n":
                header = " " # add a space if this is an empty line, just so that something is there for the user to "see"
            else:
                header = line[:mh+2] # protein name and spaces; consider making this just the protein name
            ##### ADD ID #####
            svg += f"<tspan sodipodi:role=\"line\"><tspan>{header}</tspan>"
            for c in line[mh+2:]:
                svg += create_tspan(c)
            svg += f"</tspan>"
        svg += f"</text></svg>"
        
        out = f"{out_prefix}{i}.svg"
        with open(out, "w") as o:
            print(f"Writing part {i} to {out}.", flush=True)
            o.write(svg)

    ##### add active site annotation/other features
    ##### consider grouping like-formatted things; not necessary, but does save on space; make sure it still allows for individual character formatting
    ##### consider not hard-coding centering (don't set x, set center to page)
    ##### add other features back in (like input options)
    ##### add API access to UniProt/database annotations/tabular annotations
    
def read_alignment(clust_file, max_header=16):
    ''' max_header is the max length (before 2 spaces are added) of the protein name + spaces.
    example: max_header=10; protein name + spaces = 12
    If the given protein name exceeds this value, it will be truncated.
    '''
    alignment = pc.Alignment() # empty alignment
    protein_key = {} # does not contain protein objects; just names and aligned seqs
    codes = ""
    header_len = 1000 # set so that it just uses an empty line if there is no header_len set
    with open(clust_file, "r") as clust:
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

def main():
    # add descriptions to arguments later
    parser = argparse.ArgumentParser()
    # change clust_path arguments later
    parser.add_argument("clust_file", nargs="?")    
    parser.add_argument("-out_prefix", nargs="?", default="./")
    parser.add_argument("-codes", nargs="?", default="FALSE")
    parser.add_argument("-nums", nargs="?", default="FALSE")
    parser.add_argument("-uniprot_format", nargs="?", default="FALSE")
    args = parser.parse_args()
    
    clust_file = find_path(args.clust_file, "r")
    out_prefix = find_path(args.out_prefix, "w")
    
    alignment = read_alignment(clust_file)
    alignment.set_disp_names(uniprot_format=args.uniprot_format)

    create_svg(alignment, out_prefix, codes=args.codes, nums=args.nums)

if __name__ == "__main__":
    main()

