import numpy as np
import math
import re

# Creates a single or double helix, well-suited towards nucleic acid structures.
# collect_name: the name of the collection to be created with these objects
# objects: a list of Blender objects that each rung of the helix is chosen from
# reference: a list of strings, each corresponding to the name of an object in 'objects'
# comp_fun: a function that takes an object name as input and returns the name of a complementary object
# comp: true if the sequence should be complemented in advance
# strand: number of strands, 1 or 2
# angle: positive for right hand rule, negative for left hand rule, determines rotation per base
# radius: the distance of the center of each base from the z-axis
# tilt: to create the major / minor grooves, each base must be tilted
# rotate: the degree of rotation of each base of the complement strand relative to the reference strand
# height: the height between each pair of bases
# location: the location of the center of the bottom rung, default is cursor location
def makeHelix(collect_name, objects, reference, comp_fun, comp, strand, angle, 
radius, tilt, rotate, height, location = bpy.context.scene.cursor.location):
    # avoid making duplicate collections
    for collection in bpy.data.collections:
        if collect_name == collection.name:
            print("Error: A collection with the same name already exists.")
            return
    # make a list of the names of all candidate objects
    object_names = []
    for i in range(0, len(objects)):
        object_names.append(objects[i].name)
    # creates the complementary sequence
    seq_length = len(reference)
    complement = [None] * seq_length
    for i in range(0, seq_length):
        complement[i] = comp_fun(reference[i])
        # check for non-present elements    
        if not ((reference[i] in object_names) and (complement[i] in object_names)):
            print("Error: Elements of this sequence cannot be recognized.")
            return
    # location abbreviations
    loc_x = location[0]
    loc_y = location[1]
    loc_z = location[2]
    # create the new collection
    collection = bpy.data.collections.new(collect_name)
    bpy.context.scene.collection.children.link(collection)
    # deselect all existing objects
    for sel_obj in bpy.context.selected_objects:
        sel_obj.select_set(False)
    # if necessary, swap the reference and complement sequences
    if comp:
        temp = reference
        reference = complement
        complement = temp
    # create the first strand
    curHeight = 0
    for i in range(0, seq_length):
        ref_name = reference[i]
        format_name = '{}_R{}'.format(collect_name, i)
        ref_obj = objects[object_names.index(ref_name)]
        new_obj = ref_obj.copy()
        new_obj.data = ref_obj.data.copy()
        new_obj.name = format_name
        adjust = angle*i
        new_obj.rotation_euler[1] = -tilt
        new_obj.rotation_euler[2] = adjust - math.pi
        new_obj.location = [radius*math.cos(adjust)+loc_x, radius*math.sin(adjust)+loc_y, curHeight+loc_z]
        bpy.data.collections[collect_name].objects.link(new_obj)
        new_obj.select_set(True)
        curHeight += height
    # if only one strand was requested, quit
    if strand < 2:
        return
    # create the second strand
    curHeight = 2*radius*math.sin(tilt)
    for i in range(0, seq_length):
        com_name = complement[i]
        format_name = '{}_C{}'.format(collect_name, i)
        ref_obj = objects[object_names.index(com_name)]
        new_obj = ref_obj.copy()
        new_obj.data = ref_obj.data.copy()
        new_obj.name = format_name
        adjust = angle*i + rotate
        new_obj.rotation_euler[0] = math.pi
        new_obj.rotation_euler[1] = tilt
        new_obj.rotation_euler[2] = adjust - math.pi
        new_obj.location = [radius*math.cos(adjust)+loc_x, radius*math.sin(adjust)+loc_y, curHeight+loc_z]
        bpy.data.collections[collect_name].objects.link(new_obj)
        new_obj.select_set(True)
        curHeight += height
            
# makes a molecule
# name: the name of the final molecules formed
# attachments: an array of groups to be attached to the molecule.
#   usually of the form [group_name, distance between centers, down-euler, z-euler, bond-rotate-angle]
#   also, bond-rotate-angle is right-hand rule, pointing towards parent
# location: the location of the center of the molecule, default is cursor location
def makeMolecule(name, attachments, location = bpy.context.scene.cursor.location):
    objects = bpy.data.collections["Molecules"].all_objects
    parent_obj = bpy.data.objects.new("empty", None)
    parent_obj.location = location
    for sel_obj in bpy.context.selected_objects:
        sel_obj.select_set(False)
    object_names = []
    for i in range(0, len(objects)):
        obj_name = objects[i].name
        if name == obj_name:
            print("Error: A molecule with the same name already exists")
            return
        object_names.append(obj_name)
    format_name = name
    for i in range(0, len(attachments)):
        o = bpy.data.objects.new("empty", None)
        ref_name = attachments[i][0]
        if i > 0:
            format_name = "{}_{}".format(name, i)
        ref_obj = objects[object_names.index(ref_name)]
        new_obj = ref_obj.copy()
        new_obj.data = ref_obj.data.copy()
        new_obj.name = format_name
        new_obj.rotation_euler = [0,0,0]
        new_obj.parent = o
        o.parent = parent_obj
        new_obj.location = (0, 0, attachments[i][1])
        bpy.data.collections["Molecules"].objects.link(new_obj)
        o.rotation_euler[1] = attachments[i][2]
        o.rotation_euler[2] = attachments[i][3]
        new_obj.rotation_euler[2] = attachments[i][4]
        new_obj.select_set(True)
        if i == 0:
            bpy.context.view_layer.objects.active = new_obj
    bpy.ops.object.join()

# complement function for DNA
def complementDNA(name):
    if name == 'A':
        return 'T'
    if name == 'T':
        return 'A'
    if name == 'G':
        return 'C'
    return 'G'

# complement function for RNA
def complementRNA(name):
    if name == 'A':
        return 'U'
    if name == 'U':
        return 'A'
    if name == 'G':
        return 'C'
    return 'G'

# reverses a list, usually the reference sequence
def reverse(reference):
    return reference[::-1]

# checks if seq is a valid sequence of nucleotides
# Aa, Cc, Gg represent deoxyribose and ribose sugars, respectively
def isValidSeq(seq):
    return bool(re.match("^[AaCcGgTU]*$", seq))

# converts seq, a string, to either DNA or RNA
# accepted forms for DNA: d, DNA, dna
# anything else is converted to RNA
def convertType(seq, type):
    if (type == "DNA" or type == "d" or type == "dna"):
        return seq.replace("U", "T")
    return seq.replace("T", "U")

# H (Hydrogen): 120 pm
# O (Oxygen): 152 pm
# N (Nitrogen): 155 pm
# C (Carbon): 170 pm
# P (Phosphorus): 195 pm
h_rad = 0.120
o_rad = 0.152
n_rad = 0.155
c_rad = 0.170
p_rad = 0.195

# frequent bond angles
tetra = radians(109.47)
trigo = radians(120.00)
penta = radians(108.00)
oppos = radians(180.00)

# let's make some molecules!
makeMolecule("Phosphate", [
    ["Phosphorus", 0, 0, 0, 0],
    ["Oxygen", p_rad + o_rad, tetra, 0, 0],
    ["Oxygen", p_rad + o_rad, tetra, +trigo, 0],
    ["Oxygen", p_rad + o_rad, tetra, -trigo, 0]
])

makeMolecule("OPhosphate", [
    ["Oxygen", 0, 0, 0, 0],
    ["Phosphate", p_rad + o_rad, radians(104.45), 0, trigo/4]
])

makeMolecule("CarbPho", [
    ["Carbon", 0, 0, 0, 0],
    ["OPhosphate", c_rad + o_rad, tetra, 0, 0],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0],
    ["Hydrogen", c_rad + h_rad, tetra, -trigo, 0]
])

makeMolecule("CarbPen", [
    ["CarbPho", c_rad + c_rad, penta, 0, -trigo],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0],
    ["Oxygen", c_rad + o_rad, tetra, -trigo, 0]
])

makeMolecule("CarbP2", [
    ["CarbPen", c_rad + c_rad, penta, 0, trigo/2],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0]
]) 

makeMolecule("CarbP3", [
    ["CarbP2", c_rad + c_rad, penta, 0, oppos],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0],
    ["Hydrogen", c_rad + h_rad, tetra, -trigo, 0]
])

makeMolecule("CarbP4", [
    ["CarbP3", c_rad + c_rad, 150*math.pi/180, 0, -trigo/2],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0]
])

makeMolecule("Methyl", [
    ["Hydrogen", c_rad + h_rad, tetra, 0, 0],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0],
    ["Hydrogen", c_rad + h_rad, tetra, -trigo, 0]
])

makeMolecule("Nitro", [
    ["Hydrogen", n_rad + h_rad, trigo, 0, 0],
    ["Hydrogen", n_rad + h_rad, trigo, oppos, 0]
])

makeMolecule("HC_Vert", [
    ["Hydrogen", c_rad + h_rad, oppos - penta/2, 0, 0]
])

makeMolecule("NHC_Vert", [
    ["HC_Vert", c_rad + n_rad, penta, 0, 0]
])

makeMolecule("CNHC_Vert", [
    ["NHC_Vert", c_rad + n_rad, trigo*1.115, 0, 0]
])

makeMolecule("A_CNHC", [
    ["CNHC_Vert", c_rad + c_rad, trigo, oppos, 0],
    ["Nitro", c_rad + n_rad, trigo, 0, 0]
])

makeMolecule("A1N", [
    ["A_CNHC", n_rad + c_rad, trigo, oppos, 0]
])

makeMolecule("A1CN", [
    ["A1N", n_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", n_rad + h_rad, trigo, 0, 0]
])

makeMolecule("A1NCN", [
    ["A1CN", n_rad + c_rad, trigo*1.035, oppos, 0]
])

makeMolecule("A2C", [
    ["A1NCN", n_rad + c_rad, trigo, oppos, oppos]
])

makeMolecule("HC", [
    ["Hydrogen", c_rad + h_rad, trigo, 0, 0]
])

makeMolecule("HCC", [
    ["HC", c_rad + c_rad, trigo, oppos, 0],
    ["Methyl", c_rad + c_rad, trigo, 0, trigo/2]
])

makeMolecule("HCCC", [
    ["HCC", c_rad + c_rad, trigo, oppos, 0],
    ["Oxygen", c_rad + o_rad, trigo, 0, 0]
])

makeMolecule("HCCCN", [
    ["HCCC", n_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", n_rad + h_rad, trigo, 0, 0]
])

makeMolecule("H2CN", [
    ["HCCCN", n_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", c_rad + h_rad, trigo, 0, 0]
])

makeMolecule("UCC", [
    ["HC", c_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", c_rad + h_rad, trigo, 0, trigo/2]
])

makeMolecule("UCCC", [
    ["UCC", c_rad + c_rad, trigo, oppos, 0],
    ["Oxygen", c_rad + o_rad, trigo, 0, 0]
])

makeMolecule("UCCCN", [
    ["UCCC", n_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", n_rad + h_rad, trigo, 0, 0]
])

makeMolecule("U2CN", [
    ["UCCCN", n_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", c_rad + h_rad, trigo, 0, 0]
])

makeMolecule("YC", [
    ["HC", c_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", c_rad + h_rad, trigo, 0, 0]
])

makeMolecule("YCC", [
    ["YC", c_rad + c_rad, trigo, oppos, 0],
    ["Nitro", c_rad + o_rad, trigo, 0, 0]
])

makeMolecule("YN", [
    ["YCC", n_rad + c_rad, trigo, oppos, 0]
])

makeMolecule("YCN", [
    ["YN", n_rad + c_rad, trigo, oppos, 0],
    ["Oxygen", c_rad + o_rad, trigo, 0, 0]
])

makeMolecule("G_CNHC", [
    ["CNHC_Vert", c_rad + c_rad, trigo, oppos, 0],
    ["Oxygen", c_rad + o_rad, trigo, 0, 0]
])

makeMolecule("G1N", [
    ["G_CNHC", n_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", n_rad + h_rad, trigo, 0, 0]
])

makeMolecule("G1CN", [
    ["G1N", c_rad + n_rad, trigo, oppos, 0],
    ["Nitro", c_rad + n_rad, trigo, 0, 0]
])

makeMolecule("G1NCN", [
    ["G1CN", n_rad + c_rad, trigo*1.035, oppos, 0]
])

makeMolecule("G2C", [
    ["G1NCN", n_rad + c_rad, trigo, oppos, oppos]
])

makeMolecule("A2", [
    ["A2C2", n_rad + c_rad, oppos/2, 0-10*math.pi/180, oppos/2],
    ["test2", n_rad + c_rad, oppos*3/4, oppos, -oppos]
])

makeMolecule("G2", [
    ["G2C2", n_rad + c_rad, oppos/2, 0-12*math.pi/180, oppos/2],
    ["test2", n_rad + c_rad, oppos*3/4, oppos, -oppos]
])

makeMolecule("T2", [
    ["H2CN2", n_rad + c_rad, oppos/2, -trigo/2+33*math.pi/180, oppos/2+oppos],
    ["test2", n_rad + c_rad, oppos*3/4, oppos, oppos]
])

makeMolecule("U2", [
    ["U2CN2", n_rad + c_rad, oppos/2, -trigo/2+33*math.pi/180, oppos/2+oppos],
    ["test2", n_rad + c_rad, oppos*3/4, oppos, oppos]
])

makeMolecule("C2", [
    ["YCN2", n_rad + c_rad, oppos/2, -trigo/2+31*math.pi/180, oppos/2+oppos],
    ["test2", n_rad + c_rad, oppos*3/4, oppos, oppos]
])

pen_attach = [
    ["test", 0.3 + c_rad, 0.75*math.pi, math.pi, oppos]
]

makeMolecule("A", pen_attach)
makeMolecule("C", pen_attach)
makeMolecule("G", pen_attach)
makeMolecule("T", pen_attach)
makeMolecule("U", pen_attach)

# let's make some helices!
DNAseq = "CCGATCGCGGGGAGGGGTCAAGCGCGACGTTCGCCACGCAAAACCGGGCACCCCGACTAACCGTACCGGCGCTACGCCCATTTACGAACGTAAGCGGTAA"
DNA_list = list(DNAseq)
RNAseq = convertType(DNAseq, "RNA")
RNA_list = list(RNAseq)
backbone = 'D' * len(DNA_list)

# cool looking B-DNA that reflects reality
# 10.5 bases per rotation
# 10 angstrom (1 nm) in radius, solve r^2 + r^2 - 2^2 = 2*r*r*cos 135 to get 1.1
# right-handed
# 2 strands
# 1.2 degree base tilt
# 135 degrees rotated from reference to complement
# 3x vertical scaling to make viewing easier, 0.332 nm height/base pair
makeHelix("Stylish DNA Test 2", bpy.context.selected_objects, DNA_list, complementDNA, False, 2, 
2*math.pi/10.5, 1.10, math.pi/180*1.2, math.pi/180*135, 0.332*3, [0,0,0])

makeHelix("Stylish DNA Backbone", backbone, complementBack, False, False, 2, 
2*math.pi/10.5, 1.15, math.pi/180*1.2, math.pi/180*135, 0.332*3, True)


# cool looking, not necessarily accurate
makeHelix("Stylish RNA", RNA_list, complementRNA, False, False, 1, -2*math.pi/10.5/3, 3*3, math.pi/180*1.2*18)
makeHelix("Stylish RNA 2", RNA_list, complementRNA, False, False, 1, -2*math.pi/10.5/3, 3*3, 0)
# unraveled RNA
makeHelix("Stylish RNA 3", RNA_list, complementRNA, False, False, 1, 0, 3*3, 0, 0)

#    
# = []
