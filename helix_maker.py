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
    # location abbreviations
    loc_x = location[0]
    loc_y = location[1]
    loc_z = location[2]
    # create the new collection
    collection = bpy.data.collections.new(collect_name)
    bpy.context.scene.collection.children.link(collection)
    # make a list of the names of all candidate objects
    object_names = []
    for i in range(0, len(objects)):
        object_names.append(objects[i].name)
    # creates the complementary sequence
    seq_length = len(reference)
    complement = [None] * seq_length
    for i in range(0, seq_length):
        complement[i] = comp_fun(reference[i])
    # if necessary, swap the reference and complement sequences
    if comp:
        temp = reference
        reference = complement
        complement = temp
    # create the first strand
    curHeight = 0
    for i in range(0, seq_length):
        ref_name = reference[i]
        format_name = 'Reference_{}_{}'.format(i, ref_name)
        ref_obj = objects[object_names.index(ref_name)]
        new_obj = ref_obj.copy()
        new_obj.data = ref_obj.data.copy()
        new_obj.name = format_name
        adjust = angle*i
        new_obj.rotation_euler[1] = -tilt
        new_obj.rotation_euler[2] = adjust - math.pi
        new_obj.location = [radius*math.cos(adjust)+loc_x, radius*math.sin(adjust)+loc_y, curHeight+loc_z]
        bpy.data.collections[collect_name].objects.link(new_obj)
        curHeight += height
    # if only one strand was requested, quit
    if strand < 2:
        return
    # create the second strand
    curHeight = 2*radius*math.sin(tilt)
    for i in range(0, seq_length):
        com_name = complement[i]
        format_name = 'Complement_{}_{}'.format(i, com_name)
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
    for i in range(0, len(bpy.context.selected_objects)):
        bpy.context.selected_objects[i].select_set(False)
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
            format_name = "Attachment_{}".format(i)
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
        return 'U_C'
    if name == 'U':
        return 'A_C'
    if name == 'G':
        return 'C'
    return 'G'

# reverses a list, usually the reference sequence
def reverse(reference):
    return reference[::-1]

# checks if seq is a valid sequence of nucleotides
def isValidSeq(seq):
    return bool(re.match("^[aAgGtTcCuU]*$", seq))

# converts seq, a string, to either DNA or RNA
# accepted forms for DNA: d, DNA, dna
# anything else is converted to RNA
def convertType(seq, type):
    if (type == "DNA" or type == "d" or type == "dna"):
        return seq.replace("U", "T")
    return seq.replace("T", "U")

# van der Waals radii are used
# H (Hydrogen): 120 pm
h_rad = 0.120
# O (Oxygen): 152 pm
o_rad = 0.152
# N (Nitrogen): 155 pm
n_rad = 0.155
# C (Carbon): 170 pm
c_rad = 0.170
# P (Phosphorus): 195 pm
p_rad = 0.195

# frequent bond angles
tetra = 109.47 * math.pi / 180
trigo = 120.00 * math.pi / 180
penta = 108.00 * math.pi / 180
o_ang = 104.45 * math.pi / 180
oppos = 180 * math.pi / 180

# let's make some molecules!
makeMolecule("Phosphate", [
    ["Phosphorus", 0, 0, 0, 0],
    ["Oxygen", p_rad + o_rad, tetra, 0, 0],
    ["Oxygen", p_rad + o_rad, tetra, +trigo, 0],
    ["Oxygen", p_rad + o_rad, tetra, -trigo, 0]
])

makeMolecule("OPhosphate", [
    ["Oxygen", 0, 0, 0, 0],
    ["Phosphate", p_rad + o_rad, o_ang, 0, trigo/4]
])

c_to_opho = [
    ["OPhosphate", c_rad + o_rad, tetra, 0, 0],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0],
    ["Hydrogen", c_rad + h_rad, tetra, -trigo, 0]
]

makeMolecule("CarbPho", c_to_opho)

c_to_pen = [
    ["CarbPho", c_rad + c_rad, penta, 0, -trigo],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0],
    ["Oxygen", c_rad + o_rad, tetra, -trigo, 0]
]

makeMolecule("CarbPen", c_to_pen)

pen2 = [
    ["CarbPen", c_rad + c_rad, penta, 0, trigo/2],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0]
]

makeMolecule("CarbP2", pen2)

pen3 = [
    ["CarbP2", c_rad + c_rad, penta, 0, oppos],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0],
    ["Hydrogen", c_rad + h_rad, tetra, -trigo, 0]
]

makeMolecule("CarbP3", pen3)

pen4 = [
    ["CarbP3", c_rad + c_rad, 150*math.pi/180, 0, -trigo/2],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0]
]

makeMolecule("CarbP4", pen4)

methyl_attach = [
    ["Hydrogen", c_rad + h_rad, tetra, 0, 0],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0],
    ["Hydrogen", c_rad + h_rad, tetra, -trigo, 0]
]

makeMolecule("Methyl", methyl_attach)

nitro_attach = [
    ["Hydrogen", n_rad + h_rad, trigo, 0, 0],
    ["Hydrogen", n_rad + h_rad, trigo, oppos, 0]
]

makeMolecule("Nitro", nitro_attach)

hc_vert = [
    ["Hydrogen", c_rad + h_rad, oppos - penta/2, 0, 0]
]

makeMolecule("HC_Vert", hc_vert)

nhc_vert = [
    ["HC_Vert", c_rad + n_rad, penta, 0, 0]
]

makeMolecule("NHC_Vert", nhc_vert)

cnhc_vert = [
    ["NHC_Vert", c_rad + n_rad, trigo*1.115, 0, 0]
]

makeMolecule("CNHC_Vert", cnhc_vert)

a_cnhc = [
    ["CNHC_Vert", c_rad + c_rad, trigo, oppos, 0],
    ["Nitro", c_rad + n_rad, trigo, 0, 0]
]

makeMolecule("A_CNHC", a_cnhc)

a1n = [
    ["A_CNHC", n_rad + c_rad, trigo, oppos, 0]
]

makeMolecule("A1N", a1n)

a1cn = [
    ["A1N", n_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", n_rad + h_rad, trigo, 0, 0]
]

makeMolecule("A1CN", a1cn)

a1ncn = [
    ["A1CN", n_rad + c_rad, trigo*1.035, oppos, 0]
]

makeMolecule("A1NCN", a1ncn)

a2c = [
    ["A1NCN", n_rad + c_rad, trigo, oppos, oppos]
]

makeMolecule("A2C", a2c)

hc = [
    ["Hydrogen", c_rad + h_rad, trigo, 0, 0]
]

makeMolecule("HC", hc)

hcc = [
    ["HC", c_rad + c_rad, trigo, oppos, 0],
    ["Methyl", c_rad + c_rad, trigo, 0, trigo/2]
]

makeMolecule("HCC", hcc)

hccc = [
    ["HCC", c_rad + c_rad, trigo, oppos, 0],
    ["Oxygen", c_rad + o_rad, trigo, 0, 0]
]

makeMolecule("HCCC", hccc)

hcccn = [
    ["HCCC", n_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", n_rad + h_rad, trigo, 0, 0]
]

makeMolecule("HCCCN", hcccn)

h2cn = [
    ["HCCCN", n_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", c_rad + h_rad, trigo, 0, 0]
]

makeMolecule("H2CN", h2cn)

ucc = [
    ["HC", c_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", c_rad + h_rad, trigo, 0, trigo/2]
]

makeMolecule("UCC", ucc)

uccc = [
    ["UCC", c_rad + c_rad, trigo, oppos, 0],
    ["Oxygen", c_rad + o_rad, trigo, 0, 0]
]

makeMolecule("UCCC", uccc)

ucccn = [
    ["UCCC", n_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", n_rad + h_rad, trigo, 0, 0]
]

makeMolecule("UCCCN", ucccn)

u2cn = [
    ["UCCCN", n_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", c_rad + h_rad, trigo, 0, 0]
]

makeMolecule("U2CN", u2cn)

yc = [
    ["HC", c_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", c_rad + h_rad, trigo, 0, 0]
]

makeMolecule("YC", yc)

ycc = [
    ["YC", c_rad + c_rad, trigo, oppos, 0],
    ["Nitro", c_rad + o_rad, trigo, 0, 0]
]

makeMolecule("YCC", ycc)

yn = [
    ["YCC", n_rad + c_rad, trigo, oppos, 0]
]

makeMolecule("YN", yn)

ycn = [
    ["YN", n_rad + c_rad, trigo, oppos, 0],
    ["Oxygen", c_rad + o_rad, trigo, 0, 0]
]

makeMolecule("YCN", ycn)

g_cnhc = [
    ["CNHC_Vert", c_rad + c_rad, trigo, oppos, 0],
    ["Oxygen", c_rad + o_rad, trigo, 0, 0]
]

makeMolecule("G_CNHC", g_cnhc)

g1n = [
    ["G_CNHC", n_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", n_rad + h_rad, trigo, 0, 0]
]

makeMolecule("G1N", g1n)

g1cn = [
    ["G1N", c_rad + n_rad, trigo, oppos, 0],
    ["Nitro", c_rad + n_rad, trigo, 0, 0]
]

makeMolecule("G1CN", g1cn)

g1ncn = [
    ["G1CN", n_rad + c_rad, trigo*1.035, oppos, 0]
]

makeMolecule("G1NCN", g1ncn)

g2c = [
    ["G1NCN", n_rad + c_rad, trigo, oppos, oppos]
]

makeMolecule("G2C", g2c)

wholeA = [
    ["A2C2", n_rad + c_rad, oppos/2, 0-10*math.pi/180, oppos/2],
    ["test2", n_rad + c_rad, oppos*3/4, oppos, -oppos]
]

makeMolecule("A2", wholeA)

wholeG = [
    ["G2C2", n_rad + c_rad, oppos/2, 0-12*math.pi/180, oppos/2],
    ["test2", n_rad + c_rad, oppos*3/4, oppos, -oppos]
]

makeMolecule("G2", wholeG)

wholeT = [
    ["H2CN2", n_rad + c_rad, oppos/2, -trigo/2+33*math.pi/180, oppos/2+oppos],
    ["test2", n_rad + c_rad, oppos*3/4, oppos, oppos]
]

makeMolecule("T2", wholeT)

wholeU = [
    ["U2CN2", n_rad + c_rad, oppos/2, -trigo/2+33*math.pi/180, oppos/2+oppos],
    ["test2", n_rad + c_rad, oppos*3/4, oppos, oppos]
]

makeMolecule("U2", wholeU)

wholeC = [
    ["YCN2", n_rad + c_rad, oppos/2, -trigo/2+31*math.pi/180, oppos/2+oppos],
    ["test2", n_rad + c_rad, oppos*3/4, oppos, oppos]
]

makeMolecule("C2", wholeC)

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
makeHelix("Stylish DNA Test", DNA_list, complementDNA, False, False, 2, 
2*math.pi/10.5, 1.10, math.pi/180*1.2, math.pi/180*135, 0.332*3)

makeHelix("Stylish DNA Backbone", backbone, complementBack, False, False, 2, 
2*math.pi/10.5, 1.15, math.pi/180*1.2, math.pi/180*135, 0.332*3, True)


# cool looking, not necessarily accurate
makeHelix("Stylish RNA", RNA_list, complementRNA, False, False, 1, -2*math.pi/10.5/3, 3*3, math.pi/180*1.2*18)
makeHelix("Stylish RNA 2", RNA_list, complementRNA, False, False, 1, -2*math.pi/10.5/3, 3*3, 0)
# unraveled RNA
makeHelix("Stylish RNA 3", RNA_list, complementRNA, False, False, 1, 0, 3*3, 0, 0)

#    
# = []
