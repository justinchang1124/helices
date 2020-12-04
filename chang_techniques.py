import numpy as np
import math
import re
import mathutils
import bpy

# ----------------
# HELPER FUNCTIONS
# ----------------

# complement function for DNA
def complementDNA(name):
    if name == 'A1':
        return 'T1'
    if name == 'T1':
        return 'A1'
    if name == 'G1':
        return 'C1'
    if name == 'C1':
        return 'G1'
    if name == 'A2':
        return 'T2'
    if name == 'T2':
        return 'A2'
    if name == 'G2':
        return 'C2'
    if name == 'C2':
        return 'G2'
    if name == 'A3':
        return 'T3'
    if name == 'T3':
        return 'A3'
    if name == 'G3':
        return 'C3'
    if name == 'C3':
        return 'G3'

# complement function for RNA
def complementRNA(name):
    if name == 'a1':
        return 'U1'
    if name == 'U1':
        return 'a1'
    if name == 'g1':
        return 'c1'
    if name == 'c1':
        return 'g1'
    if name == 'a2':
        return 'U2'
    if name == 'U2':
        return 'a2'
    if name == 'g2':
        return 'c2'
    if name == 'c2':
        return 'g2'
    if name == 'a3':
        return 'U3'
    if name == 'U3':
        return 'a3'
    if name == 'g3':
        return 'c3'
    if name == 'c3':
        return 'g3'

# checks if seq is a valid sequence of nucleotides
# Aa, Cc, Gg represent the corresponding bases with deoxyribose and ribose sugars, respectively
# A: Adenine
# C: Cytosine
# G: Guanine
# T: Thymine
# U: Uracil
def isValidSeq(seq):
    return bool(re.match("^[AaCcGgTU]*$", seq))

# adds a suffix to each term of a sequence
def addEnd(seq, suffix):
    s = [None] * len(seq)
    for i in range(0,len(seq)):
        s[i] = "{}{}".format(seq[i], suffix)
    return s

# converts a DNA sequence to an RNA sequence
def DNAtoRNA(seq):
    return seq.replace("A", "a").replace("C", "c"
    ).replace("G", "g").replace("T", "U")

# reverses a list, usually the reference sequence (not used in this project)
def reverse(reference):
    return reference[::-1]

# the identity function
def identity(x):
    return x

# -----------
# TECHNIQUE 1
# -----------

# Creates a single or double helix, well-suited towards nucleotide biopolymers.
# collect_name: the name of the collection to be created with these objects
# objects: a list of Blender objects that each rung of the helix is chosen from
# reference: a list of strings, each corresponding to the name of an object in 'objects'
# comp_fun: a function that takes an object name as input and returns the name of a complementary object
# comp: true if the sequence should be complemented in advance
# strand: number of strands, 1 or 2
# radius: the distance of the center of each base from the z-axis
# height: the height between each pair of bases
# base: angular rotation per base about z-axis, positive for right hand helix, negative for left hand helix
# tilt: to create the major / minor grooves, each base must be tilted by this angle measure
# disp: the angular displacement of each base of the complement strand relative to the reference strand
# location: the location of the center of the bottom rung, default is cursor location
def makeHelix(collect_name, objects, reference, comp_fun, comp, strand, radius, height, 
base, tilt, disp, location = bpy.context.scene.cursor.location):
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
        adjust = base*i
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
        adjust = base*i + disp
        new_obj.rotation_euler[0] = math.pi
        new_obj.rotation_euler[1] = tilt
        new_obj.rotation_euler[2] = adjust - math.pi
        new_obj.location = [radius*math.cos(adjust)+loc_x, radius*math.sin(adjust)+loc_y, curHeight+loc_z]
        bpy.data.collections[collect_name].objects.link(new_obj)
        new_obj.select_set(True)
        curHeight += height

# -----------
# TECHNIQUE 2
# -----------

# joins all selected objects, with the final location being the name of the parent obj
def joinSelected(name = None):
    selected = bpy.context.selected_objects
    parent_obj = None
    if len(selected) < 1:
            print("Error: No objects have been selected.")
            return
    if name != None:
        for i in range(0, len(selected)):
            if selected[i].name == name:
                parent_obj = selected[i]
    bpy.context.view_layer.objects.active = parent_obj
    bpy.ops.object.join()

# makes a molecule by appending several attachments to a parent atom
# note: joining 
# name: the name of the final molecule formed
# parent_name: the name of the atom at the center of the molecule
# attachments: an array of groups to be attached to the molecule.
#   usually of the form [group_name, distance between centers, down-euler, z-euler, bond-rotate-angle]
#   note: bond-rotate-angle obeys the right-hand rule, pointing towards the parent atom
# location: the location of the center of the molecule, default is cursor location
def makeMolecule(name, parent_name, attachments, location, join = True):
    objects = bpy.data.collections["Molecules"].all_objects
    for sel_obj in bpy.context.selected_objects:
        sel_obj.select_set(False)
    # checking object names
    object_names = []
    for i in range(0, len(objects)):
        obj_name = objects[i].name
        if name == obj_name:
            print("Error: A molecule with the same name already exists")
            return
        object_names.append(obj_name)
    if not (parent_name in object_names):
        print("Error: No parent group can be found.")
        return
    # create parent object
    parent_ref_obj = objects[object_names.index(parent_name)]
    parent_obj = copy_rec(parent_ref_obj, None)
    parent_obj.location = location
    parent_obj.name = name
    parent_obj.select_set(True)
    # add attachments
    for i in range(0, len(attachments)):
        ref_name = attachments[i][0]
        format_name = "{}_{}".format(name, i)
        ref_obj = objects[object_names.index(ref_name)]
        new_obj = copy_rec(ref_obj, None)
        new_obj.location = (0, 0, attachments[i][1])
        new_obj.name = format_name
        new_obj.rotation_euler = [math.pi,0,0]
        o = bpy.data.objects.new("empty", None)
        o.empty_display_size = 3*attachments[i][1]
        o.empty_display_type = 'SINGLE_ARROW'
        new_obj.parent = o
        o.parent = parent_obj
        bpy.data.collections["Molecules"].objects.link(o)
        o.rotation_euler[1] = attachments[i][2]
        o.rotation_euler[2] = attachments[i][3]
        new_obj.rotation_euler[2] = attachments[i][4]
        new_obj.select_set(True)
    if join:
        joinSelected(name)

def copy_rec(ref, parent):
    copy = ref.copy()
    copy.parent = parent
    bpy.data.collections["Molecules"].objects.link(copy)
    for child in ref.children:
        copy_rec(child, copy)
    return copy

tree_copy(bpy.context.selected_objects[0], None)

# -------------------
# MOLECULE GENERATION
# -------------------

# assigns radians to the math.radians function
radians = math.radians

# Note: atoms are represented by colored spheres, with radius equal to their Van der Waals radius
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

# Phosphate-Sugar Backbone
makeMolecule("Phosphate", "Phosphorus", [
    ["Oxygen", p_rad + o_rad, tetra, 0, 0],
    ["Oxygen", p_rad + o_rad, tetra, +trigo, 0],
    ["Oxygen", p_rad + o_rad, tetra, -trigo, 0]
], [0,48,0], False)

makeMolecule("OPhosphate", "Oxygen", [
    ["Phosphate", p_rad + o_rad, radians(104.45), 0, trigo/4]
], [0,52,0], False)

makeMolecule("CarbPho", "Carbon", [
    ["OPhosphate", c_rad + o_rad, tetra, 0, 0],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0],
    ["Hydrogen", c_rad + h_rad, tetra, -trigo, 0]
], [0,56,0], False)

makeMolecule("CarbPen", "Carbon", [
    ["CarbPho", c_rad + c_rad, penta, 0, -trigo],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0],
    ["Oxygen", c_rad + o_rad, tetra, -trigo, 0]
], [0,60,0], False)

makeMolecule("CarbP2", "Carbon", [
    ["CarbPen", c_rad + c_rad, penta, 0, trigo/2],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0]
], [0,64,0], False) 

makeMolecule("CarbP3", "Carbon", [
    ["CarbP2", c_rad + c_rad, penta, 0, oppos],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0],
    ["Hydrogen", c_rad + h_rad, tetra, -trigo, 0]
], [0,68,0], False)

makeMolecule("CarbP4", "Carbon", [
    ["CarbP3", c_rad + c_rad, 150*math.pi/180, 0, -trigo/2],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0]
], [0,72,0], False)

makeMolecule("Hydroxyl", "Oxygen", [
    ["Hydrogen", o_rad+h_rad, radians(104.45), 0, 0]
], [0,76,0], False)

makeMolecule("RibP3", "Carbon", [
    ["CarbP2", c_rad + c_rad, penta, 0, oppos],
    ["Hydroxyl", c_rad + o_rad, tetra, trigo, 0],
    ["Hydrogen", c_rad + h_rad, tetra, -trigo, 0]
], [0,80,0], False)

makeMolecule("RibP4", "Carbon", [
    ["RibP3", c_rad + c_rad, 150*math.pi/180, 0, -trigo/2],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0]
], [0,84,0], False)

# Adenine
makeMolecule("Methyl", "Carbon", [
    ["Hydrogen", c_rad + h_rad, tetra, 0, 0],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0],
    ["Hydrogen", c_rad + h_rad, tetra, -trigo, 0]
], [4,48,0])

makeMolecule("Nitro", "Nitrogen", [
    ["Hydrogen", n_rad + h_rad, trigo, 0, 0],
    ["Hydrogen", n_rad + h_rad, trigo, oppos, 0]
], [4,52,0])

makeMolecule("HC_Vert", "Carbon", [
    ["Hydrogen", c_rad + h_rad, oppos - penta/2, 0, 0]
], [4,56,0])

makeMolecule("NHC_Vert", "Nitrogen", [
    ["HC_Vert", c_rad + n_rad, penta, 0, 0]
], [4,60,0])

makeMolecule("CNHC_Vert", "Carbon", [
    ["NHC_Vert", c_rad + n_rad, trigo*1.115, 0, 0]
], [4,64,0])

makeMolecule("A_CNHC", "Carbon", [
    ["CNHC_Vert", c_rad + c_rad, trigo, oppos, 0],
    ["Nitro", c_rad + n_rad, trigo, 0, 0]
], [4,68,0])

makeMolecule("A1N", "Nitrogen", [
    ["A_CNHC", n_rad + c_rad, trigo, oppos, 0]
], [4,72,0])

makeMolecule("A1CN", "Carbon", [
    ["A1N", n_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", n_rad + h_rad, trigo, 0, 0]
], [4,76,0])

makeMolecule("A1NCN", "Nitrogen", [
    ["A1CN", n_rad + c_rad, trigo*1.035, oppos, 0]
], [4,80,0])

makeMolecule("A2C", "Carbon", [
    ["A1NCN", n_rad + c_rad, trigo, oppos, oppos]
], [4,84,0])

# Thymine
makeMolecule("HC", "Carbon", [
    ["Hydrogen", c_rad + h_rad, trigo, 0, 0]
], [8,48,0])

makeMolecule("HCC", "Carbon", [
    ["HC", c_rad + c_rad, trigo, oppos, 0],
    ["Methyl", c_rad + c_rad, trigo, 0, trigo/2]
], [8,52,0])

makeMolecule("HCCC", "Carbon", [
    ["HCC", c_rad + c_rad, trigo, oppos, 0],
    ["Oxygen", c_rad + o_rad, trigo, 0, 0]
], [8,56,0])

makeMolecule("HCCCN", "Nitrogen", [
    ["HCCC", n_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", n_rad + h_rad, trigo, 0, 0]
], [8,60,0])

makeMolecule("H2CN", "Carbon", [
    ["HCCCN", n_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", c_rad + h_rad, trigo, 0, 0]
], [8,64,0])

# Uracil
makeMolecule("UCC", "Carbon", [
    ["HC", c_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", c_rad + h_rad, trigo, 0, trigo/2]
], [12,48,0])

makeMolecule("UCCC", "Carbon", [
    ["UCC", c_rad + c_rad, trigo, oppos, 0],
    ["Oxygen", c_rad + o_rad, trigo, 0, 0]
], [12,52,0])

makeMolecule("UCCCN", "Nitrogen", [
    ["UCCC", n_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", n_rad + h_rad, trigo, 0, 0]
], [12,56,0])

makeMolecule("U2CN", "Carbon", [
    ["UCCCN", n_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", c_rad + h_rad, trigo, 0, 0]
], [12,60,0])

# Cytosine
makeMolecule("YC", "Carbon", [
    ["HC", c_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", c_rad + h_rad, trigo, 0, 0]
], [16,48,0])

makeMolecule("YCC", "Carbon", [
    ["YC", c_rad + c_rad, trigo, oppos, 0],
    ["Nitro", c_rad + o_rad, trigo, 0, 0]
], [16,52,0])

makeMolecule("YN", "Nitrogen", [
    ["YCC", n_rad + c_rad, trigo, oppos, 0]
], [16,56,0])

makeMolecule("YCN", "Carbon", [
    ["YN", n_rad + c_rad, trigo, oppos, 0],
    ["Oxygen", c_rad + o_rad, trigo, 0, 0]
], [16,60,0])

# Guanine
makeMolecule("G_CNHC", "Carbon", [
    ["CNHC_Vert", c_rad + c_rad, trigo, oppos, 0],
    ["Oxygen", c_rad + o_rad, trigo, 0, 0]
], [20,48,0])

makeMolecule("G1N", "Nitrogen", [
    ["G_CNHC", n_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", n_rad + h_rad, trigo, 0, 0]
], [20,52,0])

makeMolecule("G1CN", "Carbon", [
    ["G1N", c_rad + n_rad, trigo, oppos, 0],
    ["Nitro", c_rad + n_rad, trigo, 0, 0]
], [20,56,0])

makeMolecule("G1NCN", "Nitrogen", [
    ["G1CN", n_rad + c_rad, trigo*1.035, oppos, 0]
], [20,60,0])

makeMolecule("G2C", "Carbon", [
    ["G1NCN", n_rad + c_rad, trigo, oppos, oppos]
], [20,64,0])

# Put them together 
makeMolecule("A1", "Nitrogen", [
    ["Adenine1", n_rad + c_rad, oppos/2, 0-10*math.pi/180, oppos/2],
    ["Deoxy1", n_rad + c_rad, oppos*3/4, oppos, -oppos]
], [32,48,0])

makeMolecule("a1", "Nitrogen", [
    ["Adenine1", n_rad + c_rad, oppos/2, 0-10*math.pi/180, oppos/2],
    ["Ribose1", n_rad + c_rad, oppos*3/4, oppos, -oppos]
], [32,52,0])

makeMolecule("C1", "Nitrogen", [
    ["Cytosine1", n_rad + c_rad, oppos/2, -trigo/2+31*math.pi/180, oppos/2+oppos],
    ["Deoxy1", n_rad + c_rad, oppos*3/4, oppos, oppos]
], [32,56,0])

makeMolecule("c1", "Nitrogen", [
    ["Cytosine1", n_rad + c_rad, oppos/2, -trigo/2+31*math.pi/180, oppos/2+oppos],
    ["Ribose1", n_rad + c_rad, oppos*3/4, oppos, oppos]
], [32,60,0])

makeMolecule("G1", "Nitrogen", [
    ["Guanine1", n_rad + c_rad, oppos/2, 0-12*math.pi/180, oppos/2],
    ["Deoxy1", n_rad + c_rad, oppos*3/4, oppos, -oppos]
], [32,64,0])

makeMolecule("g1", "Nitrogen", [
    ["Guanine1", n_rad + c_rad, oppos/2, 0-12*math.pi/180, oppos/2],
    ["Ribose1", n_rad + c_rad, oppos*3/4, oppos, -oppos]
], [32,68,0])

makeMolecule("T1", "Nitrogen", [
    ["Thymine1", n_rad + c_rad, oppos/2, -trigo/2+33*math.pi/180, oppos/2+oppos],
    ["Deoxy1", n_rad + c_rad, oppos*3/4, oppos, oppos]
], [32,72,0])

makeMolecule("U1", "Nitrogen", [
    ["Uracil1", n_rad + c_rad, oppos/2, -trigo/2+33*math.pi/180, oppos/2+oppos],
    ["Ribose1", n_rad + c_rad, oppos*3/4, oppos, oppos]
], [32,76,0])

# coloring was done by hand
makeMolecule("A2", "Nitrogen", [
    ["Adenine2", n_rad + c_rad, oppos/2, 0-10*math.pi/180, oppos/2],
    ["Deoxy2", n_rad + c_rad, oppos*3/4, oppos, -oppos]
], [36,48,0])

makeMolecule("a2", "Nitrogen", [
    ["Adenine2", n_rad + c_rad, oppos/2, 0-10*math.pi/180, oppos/2],
    ["Ribose2", n_rad + c_rad, oppos*3/4, oppos, -oppos]
], [36,52,0])

makeMolecule("C2", "Nitrogen", [
    ["Cytosine2", n_rad + c_rad, oppos/2, -trigo/2+31*math.pi/180, oppos/2+oppos],
    ["Deoxy2", n_rad + c_rad, oppos*3/4, oppos, oppos]
], [36,56,0])

makeMolecule("c2", "Nitrogen", [
    ["Cytosine2", n_rad + c_rad, oppos/2, -trigo/2+31*math.pi/180, oppos/2+oppos],
    ["Ribose2", n_rad + c_rad, oppos*3/4, oppos, oppos]
], [36,60,0])

makeMolecule("G2", "Nitrogen", [
    ["Guanine2", n_rad + c_rad, oppos/2, 0-12*math.pi/180, oppos/2],
    ["Deoxy2", n_rad + c_rad, oppos*3/4, oppos, -oppos]
], [36,64,0])

makeMolecule("g2", "Nitrogen", [
    ["Guanine2", n_rad + c_rad, oppos/2, 0-12*math.pi/180, oppos/2],
    ["Ribose2", n_rad + c_rad, oppos*3/4, oppos, -oppos]
], [36,68,0])

makeMolecule("T2", "Nitrogen", [
    ["Thymine2", n_rad + c_rad, oppos/2, -trigo/2+33*math.pi/180, oppos/2+oppos],
    ["Deoxy2", n_rad + c_rad, oppos*3/4, oppos, oppos]
], [36,72,0])

makeMolecule("U2", "Nitrogen", [
    ["Uracil2", n_rad + c_rad, oppos/2, -trigo/2+33*math.pi/180, oppos/2+oppos],
    ["Ribose2", n_rad + c_rad, oppos*3/4, oppos, oppos]
], [36,76,0])

# ----------------
# HELIX GENERATION
# ----------------

DNAseq = "ACTGACTGAAAATTTTCCCCGGGG"
DNA_list = list(DNAseq)
DNA1 = addEnd(DNA_list, "1")
DNA2 = addEnd(DNA_list, "2")
DNA3 = addEnd(DNA_list, "3")
RNAseq = DNAtoRNA(DNAseq)
RNA_list = list(RNAseq)
RNA1 = addEnd(RNA_list, "1")
RNA2 = addEnd(RNA_list, "2")
RNA3 = addEnd(RNA_list, "3")

# B-DNA that reflects reality
# 2 strands
# 10 angstrom (1 nm) in radius, solve r^2 + r^2 - 2^2 = 2*r*r*cos 135 to get 1.10
# 3x vertical scaling to make viewing easier, 0.332 nm height/base pair
# right-handed, 10.5 bases per rotation
# 1.2 degree base tilt
# 135 degrees rotated from reference to complement
makeHelix("DNA1 Test", bpy.data.collections["Molecules"].all_objects, DNA1, complementDNA, False, 2, 
1.10, 0.332*3, 2*math.pi/10.5, math.pi/180*1.2, math.pi/180*135, [48,48,0])

# same as above, but coloring
makeHelix("DNA2 Test", bpy.data.collections["Molecules"].all_objects, DNA2, complementDNA, False, 2, 
1.10, 0.332*3, 2*math.pi/10.5, math.pi/180*1.2, math.pi/180*135, [48,56,0])

# RNA corresponding to DNA
makeHelix("RNA1 Test", bpy.data.collections["Molecules"].all_objects, RNA1, complementRNA, False, 1, 
1.10, 0.332*3, 2*math.pi/10.5, math.pi/180*1.2, math.pi/180*135, [48,64,0])

# unraveled
makeHelix("RNA2 Test", bpy.data.collections["Molecules"].all_objects, RNA2, complementRNA, False, 1, 
1.10, 0.332*3, 0, math.pi/180*1.2, math.pi/180*135, [48,72,0])

# left-handed DNA is not realistic, since the backbone doesn't line up
makeHelix("DNA3 Test", bpy.data.collections["Molecules"].all_objects, DNA1, complementDNA, False, 1, 
1.10, 0.332*3, -2*math.pi/10.5, math.pi/180*1.2*20, math.pi/180*135, [48,80,0])

# unrealistic DNA - extreme twisting
makeHelix("DNA4 Test", bpy.data.collections["Molecules"].all_objects, DNA2, complementDNA, False, 2, 
1.10, 0.332*3, 2*math.pi/10.5, math.pi/180*1.2*20, math.pi/180*135, [48,88,0])

# very basic, exaggerated double-stranded left-handed RNA model for fun
makeHelix("RNA3 Test", bpy.data.collections["Molecules"].all_objects, RNA3 * 2, complementRNA, False, 2, 
1.10*3, 1.2, -2*math.pi/10.5/3, math.pi/180*1.2*20, math.pi/180*120, [48,104,0])

# double-stranded basic DNA, but now every base is paired with itself and no major-minor is present
makeHelix("DNA5 Test", bpy.data.collections["Molecules"].all_objects, DNA3, identity, False, 2, 
1.10, 1.1, 2*math.pi/10.5, math.pi/180*1.2, math.pi, [48,120,0])

class CURVE_PT_ChangFinalPanel(Panel):
    bl_label = "Nucleic Acid Techniques"
    bl_idname = "CURVE_PT_ChangFinalPanel"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Create"
    bl_context = "objectmode"
    bl_options = {"DEFAULT_OPEN"}

    def draw(self, context):
        layout = self.layout
        wm = context.window_manager
        col = layout.column(align=True)

        prop_new = col.operator("curve.ivy_gen", text="Add New Ivy", icon="OUTLINER_OB_CURVE")
        prop_new.defaultIvy = False
        prop_new.updateIvy = True

        prop_def = col.operator("curve.ivy_gen", text="Add New Default Ivy", icon="CURVE_DATA")
        prop_def.defaultIvy = True
        prop_def.updateIvy = True

        col = layout.column(align=True)
        col.label(text="Generation Settings:")
        col.prop(wm.ivy_gen_props, "randomSeed")
        col.prop(wm.ivy_gen_props, "maxTime")

        col = layout.column(align=True)
        col.label(text="Size Settings:")
        col.prop(wm.ivy_gen_props, "maxIvyLength")
        col.prop(wm.ivy_gen_props, "ivySize")
        col.prop(wm.ivy_gen_props, "maxFloatLength")
        col.prop(wm.ivy_gen_props, "maxAdhesionDistance")

        col = layout.column(align=True)
        col.label(text="Weight Settings:")
        col.prop(wm.ivy_gen_props, "primaryWeight")
        col.prop(wm.ivy_gen_props, "randomWeight")
        col.prop(wm.ivy_gen_props, "gravityWeight")
        col.prop(wm.ivy_gen_props, "adhesionWeight")

        col = layout.column(align=True)
        col.label(text="Branch Settings:")
        col.prop(wm.ivy_gen_props, "branchingProbability")
        col.prop(wm.ivy_gen_props, "ivyBranchSize")

        col = layout.column(align=True)
        col.prop(wm.ivy_gen_props, "growLeaves")

        if wm.ivy_gen_props.growLeaves:
            col = layout.column(align=True)
            col.label(text="Leaf Settings:")
            col.prop(wm.ivy_gen_props, "ivyLeafSize")
            col.prop(wm.ivy_gen_props, "leafProbability")


class ChangFinalParams(PropertyGroup):
    maxIvyLength: FloatProperty(
        name="Max Ivy Length",
        description="Maximum ivy length in Blender Units",
        default=1.0,
        min=0.0,
        soft_max=3.0,
        subtype='DISTANCE',
        unit='LENGTH'
    )
    primaryWeight: FloatProperty(
        name="Primary Weight",
        description="Weighting given to the current direction",
        default=0.5,
        min=0.0,
        soft_max=1.0
    )
    randomWeight: FloatProperty(
        name="Random Weight",
        description="Weighting given to the random direction",
        default=0.2,
        min=0.0,
        soft_max=1.0
    )
    gravityWeight: FloatProperty(
        name="Gravity Weight",
        description="Weighting given to the gravity direction",
        default=1.0,
        min=0.0,
        soft_max=1.0
    )
    adhesionWeight: FloatProperty(
        name="Adhesion Weight",
        description="Weighting given to the adhesion direction",
        default=0.1,
        min=0.0,
        soft_max=1.0
    )
    branchingProbability: FloatProperty(
        name="Branching Probability",
        description="Probability of a new branch forming",
        default=0.05,
        min=0.0,
        soft_max=1.0
    )
    leafProbability: FloatProperty(
        name="Leaf Probability",
        description="Probability of a leaf forming",
        default=0.35,
        min=0.0,
        soft_max=1.0
    )
    ivySize: FloatProperty(
        name="Ivy Size",
        description="The length of an ivy segment in Blender"
                    " Units",
        default=0.02,
        min=0.0,
        soft_max=1.0,
        precision=3
    )
    ivyLeafSize: FloatProperty(
        name="Ivy Leaf Size",
        description="The size of the ivy leaves",
        default=0.02,
        min=0.0,
        soft_max=0.5,
        precision=3
    )
    growLeaves: BoolProperty(
        name="Grow Leaves",
        description="Grow leaves or not",
        default=True
    )

def register(classes):
    # add all classes
    for cls in classes:
        bpy.utils.register_class(cls)
    # create a parameter holder in WindowManager
    bpy.types.WindowManager.chang_final_params = PointerProperty(type=ChangFinalParams)

def unregister(classes):
    # remove parameters from WindowManager
    del bpy.types.WindowManager.chang_final_params
    # remove all classes
    for cls in classes:
        bpy.utils.unregister_class(cls) 

my_classes = (
    # actual operator
    ChangTechniques,
    # parameter holder
    ChangFinalParams,
    # panel
    CURVE_PT_ChangFinalPanel
)
     
register(my_classes)
# unregister(my_classes)