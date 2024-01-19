def translate_style(style_dict):
    style_str = "" 
    for s in style_dict:
        style_str += f"{s}:{style_dict[s]},"
    return style_str

import pathlib

class Dot:
    def __init__(self, x=0, y=0, fill="000000"):
        self.x = x
        self.y = y
        self.fill = fill
    
    def svg_form(self):
        # fill in later
        return
    
class Arrow:
    def __init__(self, x=0, y=0, fill="000000"): # fix this later
        self.x = x
        self.y = y
        self.fill = fill
    def __str__(self):
        with open(f"{pathlib.Path(__file__).parent.resolve()}/svg_formats/arrow_editable.txt", "r") as svg: # fix later 
            s = ""
            for line in svg:
                s += line
        s = s.replace("X_TRANS", str(self.x))
        s = s.replace("Y_TRANS", str(self.y))
        s = s.replace("REPLACE_FILL", self.fill)
        #3aae9b GREEN
        return s
        
    