# this is the redesigned protein_objects


##### FIRST MAKE A CLEAN ALIGNMENT STRUCTURE #####

# need to deconstruct clustal output and save it in protein object
# sequences must be aligned for this to work, regardless of input format

class Protein:
    def __init__(self, name, sequence):
        self.name = name
        self.disp_name = name # unless set otherwise
        self.sequence = sequence # in aligned FASTA format, with "-" for gaps (no >name, just sequence)
        self.seq_raw = self.sequence.replace("-", "")
        self.features = []

    def add_feature(self, feature):
        clust_start = None
        clust_end = None
        
        raw_count = 1
        clust_count = 1        
        for r in self.sequence:
            if raw_count == feature.start:
                clust_start = clust_count
            elif raw_count == feature.end:
                clust_end = clust_count
                break # end the loop once it finds the end of the feature
            clust_count += 1
            if r != "-":
                raw_count += 1
                
        feature.set_clust_positions(clust_start, clust_end)
        self.features.append(feature)
    ##### end here 2023.10.18
    
    def set_disp_name(self, length=1000, uniprot_format="FALSE"):
        if uniprot_format == "TRUE":
            self.disp_name = self.name.split("|")[-1][:length]
        else:
            self.disp_name = self.name[:length]

class Alignment:
    def __init__(self):
        self.proteins = {} # dictionary of Proteins; keys are Protein.name; consider using just a list...
        # later, consider making this editable for things like active sites - maybe give proteins something to override this
        self.codes = ""
        self.max_header = None
        
    def add_protein(self, protein):
        self.proteins[protein.name] = protein
    
    def add_codes(self, codes):
        self.codes = codes
        
    def set_max_header(self, max_header):
        self.max_header = max_header
        
    def get_max_length(self):
        length = 0
        for name in self.proteins:
            length = max(len(self.proteins[name].seq_raw), length)
        return length
    
    def add_features(self, annotations):
        for row in annotations.itertuples(index=False):
            
            feature_type = row[annotations.columns.get_loc("type")]
            start = row[annotations.columns.get_loc("location.start.value")]
            end = row[annotations.columns.get_loc("location.end.value")]
            
            feature = Feature(feature_type, start, end, row.description)

            self.proteins[row.whole_prot].add_feature(feature)
    
    def get_features(self):
        ##### set/get feature coordinates(?) #####
        return feature_coords
    
    def set_disp_names(self, length=1000, uniprot_format="FALSE"):
        # Length set to 1000 just so full length is kept; change later maybe
        for name in self.proteins:
            self.proteins[name].set_disp_name(length=length, uniprot_format=uniprot_format)
        
class Feature:
    # deal with mismatched sequences later (maybe use Carter's 5res before/after
    # start and end are 1-indexed, like a normal protein sequence
    def __init__(self, feature_type, start, end, description):
        self.feature_type = feature_type
        self.start = start # position in raw sequence
        self.end = end # position in raw sequence
        self.clust_start = None
        self.clust_end = None
        self.description = description
        # consider adding self.sequence in order to identify later...
        
    def set_clust_positions(self, clust_start, clust_end): # clustal positions based on start and end
        self.clust_start = clust_start
        self.clust_end = clust_end

def pub_format(self, nums="FALSE", codes="FALSE"):
    lines_dict = {k: v.pub_format(nums) for k, v in self.proteins.items()}
    if codes == "TRUE" or codes == "T":
        codes_pub = self.codes_pub_format() # list of publication-formatted codes
    if self.max_acc_len is not None:
        #acc_len = min(self.max_acc_len, max(len(x) for x in lines_dict.keys())) # max for everything other than spaces; add 2 spaces at print time
        acc_len = self.max_acc_len
    else:
        acc_len = max(len(x) for x in lines_dict.keys())
    lines = ""
    for i in range(0, len(list(lines_dict.values())[0])): # IMPORTANT: BUG WITH .CLUSTAL FILES WHERE THE RANGE IS ONE MORE
        for p in lines_dict.items():
            lines += p[0][:acc_len] + (" " * (acc_len-len(p[0]))) + "  " # add accession and spaces up to max accession length + 2
            lines += p[1][i] # ends in a number and a newline; NOTE: THIS ISN'T RIGHT-ALIGNED - maybe fix later
        if codes == "TRUE" or codes == "T":
            lines += (" " * acc_len) + "  " + codes_pub[i] # add spaces and ends in a newline 
        else:
            lines += "\n" # extra newline, even if codes not included
        lines += "\n" # extra newline
    # need to strip final newline (see below) think of better way to do this later
    lines = "\n".join(lines.split("\n")[:-1])
    return lines 

# Example usage
alignment = Alignment()

# Create proteins and add them to the alignment
protein1 = Protein("Protein1", "SEQUENCE1")
protein2 = Protein("Protein2", "SEQUENCE2")
alignment.add_protein(protein1)
alignment.add_protein(protein2)

# Add features to the proteins
feature1 = Feature("Active Site", 10, 15, "Important active site")
feature2 = Feature("Conserved Residue", 25, 30, "Highly conserved region")
protein1.add_feature(feature1)
protein2.add_feature(feature2)
