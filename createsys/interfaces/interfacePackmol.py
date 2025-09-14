
import os



class packmolBuilder:
    def __init__(self, outfile="packmol_out.pdb"):
        self.header_dict = {
            'tolerance': 2.0,
            'filetype': 'pdb',
        }
        self.structures = []
        self.outfile = outfile

    def add_structure(self, structure):
        self.structures.append(structure)
    
    def set_seed(self, seed):
        self.header_dict["seed"] = seed

    def set_pbc(self, *args):
        args = self.convert_box_args(args)
        self.header_dict["pbc"] = ' '.join(args)

    def write_input(self, inputfilename="packmol.inp", outfile=None):
        self.header_dict['output']  = outfile or self.outfile

        input_header = "\n".join([f"{kw} {val}" for kw, val in self.header_dict.items()])

        content = input_header + "\n\n"
        content += "\n".join([struc.get_text() for struc in self.structures])
        with open(inputfilename, 'w') as file:
            file.write(content)

    def execute_packmol(self, inputfilename="packmol.inp"): 
        os.system(f"packmol < {inputfilename}")
    
    def remove_inputfile(self, inputfilename="packmol.inp"):
        os.remove(inputfilename)
    
    def run(self, outfile=None, inputfilename="packmol.inp"):
        self.write_input(inputfilename, outfile)
        self.execute_packmol(inputfilename)
        #self.remove_inputfile(inputfilename)

    @staticmethod 
    def convert_box_args(*args):
        if len(args) == 3:
            args = [0] * 3 + list(args)
        elif len(args) == 6:
            pass
        elif len(args) == 1:
            while isinstance(args, tuple):
                args = args[0]
            args = [0] * 3 + [args] * 3
        else:
            raise ValueError("Please provide 1, 3 or 6 arguments for box dimensions")
        args = [str(float(arg)) for arg in args]
        return args
    
    @staticmethod
    def convert_sphere_args(radius, *args):
        if len(args) == 1:
            args *= 3
        elif len(args) != 3:
            raise ValueError("Please provide 1 or 3 arguments for the center of the sphere.")
        args = [str(float(arg)) for arg in args]
        args += [str(float(radius))]
        return args

class packmolStructure:
    def __init__(self, filename, number=1):
        self.filename = filename
        self.number = number
        self.text = f"""    number {int(number)}
"""

    
    def add_constraint(self, where, type, *args):
        if where not in ['inside', 'outside']:
            raise NotImplementedError
        if type not in ['sphere', 'box']:
            raise NotImplementedError

        if type == 'sphere':
            dimensions = packmolBuilder.convert_sphere_args(*args)
        elif type == 'box':
            dimensions = packmolBuilder.convert_box_args(*args)
        self.text += f"    {where} {type} {' '.join(dimensions)}\n"

    def get_text(self):
        return f"structure {self.filename}\n{self.text}end structure\n"
    
    def inside_sphere(self, radius, *center):
        self.add_constraint("inside", 'sphere', radius, *center)
    def inside_box(self, *args):
        self.add_constraint("inside", 'box', *args)
    def outside_sphere(self, radius, *center):
        self.add_constraint("outside", 'sphere', radius, *center)
    def outside_box(self, *args):
        self.add_constraint("outside", 'box', *args)


