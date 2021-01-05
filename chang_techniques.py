# ---------
# LIBRARIES
# ---------

import numpy as np
import math
import re
import mathutils
import bpy
from bpy.types import (
    Operator,
    Panel,
    PropertyGroup,
)
from bpy.props import (
    BoolProperty,
    FloatProperty,
    IntProperty,
    PointerProperty,
    StringProperty
)

# assigns rad to the math.radians function
rad = math.radians

# --------------------------------
# TECHNIQUE 1: MOLECULE GENERATION
# --------------------------------

# returns the collection named "Molecules"
def molCol():
    return bpy.data.collections["Molecules"]

# deselects all provided objects
# objects: the objects to be deselected
def deselect(objects = None):
    if objects == None:
        objects = bpy.context.selected_objects
    for obj in objects:
        obj.select_set(False)
        
# toggles the visibility of all empties in the provided objects
# objects: the objects to be toggled
def hideEmpties(hide = True):
    for obj in bpy.context.scene.view_layers[0].objects:
        if obj.type == 'EMPTY':
            obj.hide_set(hide)

# copies an object and all descendants recursively
# obj: the object to be copied
# parent: the object to be assigned as the parent of this object
# collection: the collection that the copy should be added to
def copyRec(obj, parent = None, collection = None):
    if collection == None:
        collection = molCol()
    copy = obj.copy()
    if not obj.data == None:
        # this step is really important to ensure the reference isn't mutated!
        copy.data = obj.data.copy()
    copy.parent = parent
    collection.objects.link(copy)
    for child in obj.children:
        copyRec(child, copy, collection)
    return copy

# selects an object and its descendants
# obj: the object to be selected
def seleRec(obj):
    obj.select_set(True)
    for child in obj.children:
        seleRec(child)
        
# assigns a material to an object and its descendants
def coloRec(obj, material, clearExtra = True):
    if obj.type != 'EMPTY':
        if clearExtra:
            obj.data.materials.clear()
        obj.data.materials.append(material)
    for child in obj.children:
        coloRec(child, material, clearExtra)

# joins an object and its descendants, 'promoting' empties
# parent: the object that will act as the parent 
# tossEmpties: whether empties should be tossed from the join
def joinRec(parent = None):
    # by default, join from the first selected object
    if parent == None and len(bpy.context.selected_objects) > 0:
        parent = bpy.context.selected_objects[0]
    # if no children, ignore
    if len(parent.children) == 0:
        print("Error: This object has no children.")
        return
    # deselect all objects
    deselect()
    # select the parent and all children
    seleRec(parent)
    # 'promote' if you select an empty
    parent = parent.parent
    parent.select_set(True)
    # make the parent the active object
    bpy.context.view_layer.objects.active = parent
    # remove all empties from selection
    empty_list = []
    for obj in bpy.context.selected_objects:
        if obj.type == 'EMPTY':
            obj.select_set(False)
            empty_list.append(obj)
    # actually join
    bpy.ops.object.join()
    # destroy all empties
    for empty in empty_list:
        bpy.data.objects.remove(empty, do_unlink=True)

# makes a molecule by appending several attachments to a parent atom
# name: the name of the final molecule formed
# parent_name: the name of the group at the center of the molecule
# attachments: an array of groups to be attached to the molecule.
#   [group_name, distance between centers, down-euler, z-euler, bond-rotate-angle]
#   note: bond-rotate-angle obeys the right-hand rule, pointing towards the parent atom
# location: the location of the center of the molecule
# objects: the group of objects that this molecules is formed from and added to
def makeMolecule(name, parent_name, attachments = [], location = None, objects = None):
    if location == None:
        location = bpy.context.scene.cursor.location
    if objects == None:
        objects = molCol().objects
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
    parent_ref_obj = objects[parent_name]
    parent_obj = copyRec(parent_ref_obj)
    parent_obj.location = location
    parent_obj.name = name
    # add attachments
    for i in range(0, len(attachments)):
        format_name = "{}_{}".format(name, i)
        # create an empty to link the parent and its attachment
        o = bpy.data.objects.new("empty", None)
        objects.link(o)
        o.empty_display_size = 3*attachments[i][1]
        o.empty_display_type = 'SINGLE_ARROW'
        o.name = "EMPTY_{}".format(format_name)
        o.rotation_euler[1] = attachments[i][2]
        o.rotation_euler[2] = attachments[i][3]
        # copy the reference attachment to make a new object
        new_obj = copyRec(objects[attachments[i][0]])
        new_obj.location = (0, 0, attachments[i][1])
        new_obj.rotation_euler = [math.pi,0,0]
        new_obj.name = format_name
        new_obj.rotation_euler[2] = attachments[i][4]
        # make the appropriate relation
        new_obj.parent = o
        o.parent = parent_obj
    return parent_obj
    
# --------------------------
# TECHNIQUE 1: DEMONSTRATION
# --------------------------

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
tetra = rad(109.47)
trigo = rad(120.00)
penta = rad(108.00)
oppos = rad(180.00)

# Phosphate-Sugar Backbone
makeMolecule("Phosphate", "Phosphorus", [
    ["Oxygen", p_rad + o_rad, tetra, 0, 0],
    ["Oxygen", p_rad + o_rad, tetra, +trigo, 0],
    ["Oxygen", p_rad + o_rad, tetra, -trigo, 0]
], [0,48,0])

joinRec(molCol().objects["Phosphate"])

makeMolecule("OPhosphate", "Oxygen", [
    ["Phosphate", p_rad + o_rad, rad(104.45), 0, trigo/4]
], [0,52,0])

makeMolecule("CarbPho", "Carbon", [
    ["OPhosphate", c_rad + o_rad, tetra, 0, 0],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0],
    ["Hydrogen", c_rad + h_rad, tetra, -trigo, 0]
], [0,56,0])

joinRec(molCol().objects["EMPTY_CarbPho_1"])
joinRec(molCol().objects["EMPTY_CarbPho_2"])

makeMolecule("CarbPen", "Carbon", [
    ["CarbPho", c_rad + c_rad, penta, 0, -trigo],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0],
    ["Oxygen", c_rad + o_rad, tetra, -trigo, 0]
], [0,60,0])

joinRec(molCol().objects["EMPTY_CarbPen_1"])

makeMolecule("CarbP2", "Carbon", [
    ["CarbPen", c_rad + c_rad, penta, 0, trigo/2],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0]
], [0,64,0]) 

joinRec(molCol().objects["EMPTY_CarbP2_1"])

makeMolecule("CarbP3", "Carbon", [
    ["CarbP2", c_rad + c_rad, penta, 0, oppos],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0],
    ["Hydrogen", c_rad + h_rad, tetra, -trigo, 0]
], [0,68,0])

joinRec(molCol().objects["EMPTY_CarbP3_1"])
joinRec(molCol().objects["EMPTY_CarbP3_2"])

makeMolecule("CarbP4", "Carbon", [
    ["CarbP3", c_rad + c_rad, rad(150), 0, -trigo/2],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0]
], [0,72,0])

joinRec(molCol().objects["EMPTY_CarbP4_1"])

makeMolecule("Hydroxyl", "Oxygen", [
    ["Hydrogen", o_rad+h_rad, rad(104.45), 0, 0]
], [0,76,0])

joinRec(molCol().objects["EMPTY_Hydroxyl_0"])

makeMolecule("RibP3", "Carbon", [
    ["CarbP2", c_rad + c_rad, penta, 0, oppos],
    ["Hydroxyl", c_rad + o_rad, tetra, trigo, 0],
    ["Hydrogen", c_rad + h_rad, tetra, -trigo, 0]
], [0,80,0])

joinRec(molCol().objects["EMPTY_RibP3_2"])

makeMolecule("RibP4", "Carbon", [
    ["RibP3", c_rad + c_rad, 150*math.pi/180, 0, -trigo/2],
    ["Hydrogen", c_rad + h_rad, tetra, trigo, 0]
], [0,84,0])

joinRec(molCol().objects["EMPTY_RibP4_1"])

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

makeMolecule("AttachA", "Carbon", [
    ["A1NCN", n_rad + c_rad, trigo, oppos, oppos]
], [4,84,0])

joinRec(bpy.data.collections["Molecules"].objects["AttachA"])

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

makeMolecule("AttachT", "Carbon", [
    ["HCCCN", n_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", c_rad + h_rad, trigo, 0, 0]
], [8,64,0])

joinRec(bpy.data.collections["Molecules"].objects["AttachT"])

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

makeMolecule("AttachU", "Carbon", [
    ["UCCCN", n_rad + c_rad, trigo, oppos, 0],
    ["Hydrogen", c_rad + h_rad, trigo, 0, 0]
], [12,60,0])

joinRec(bpy.data.collections["Molecules"].objects["AttachU"])

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

makeMolecule("AttachC", "Carbon", [
    ["YN", n_rad + c_rad, trigo, oppos, 0],
    ["Oxygen", c_rad + o_rad, trigo, 0, 0]
], [16,60,0])

joinRec(bpy.data.collections["Molecules"].objects["AttachC"])

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

makeMolecule("AttachG", "Carbon", [
    ["G1NCN", n_rad + c_rad, trigo, oppos, oppos]
], [20,64,0])

joinRec(bpy.data.collections["Molecules"].objects["AttachG"])

# make copies of ribose, deoxy and color
# Note that you can use makeMolecule to copy molecules as well!
r1 = makeMolecule("Ribose1", "RibP4", location=[24, 48, 0])
r1.scale.y = -1

r2 = makeMolecule("Ribose2", "RibP4", location=[24, 52, 0])
r2.scale.y = -1
coloRec(r2, bpy.data.materials['Ribose'])

d1 = makeMolecule("Deoxy1", "CarbP4", location=[24, 56, 0])
d1.scale.y = -1

d2 = makeMolecule("Deoxy2", "CarbP4", location=[24, 60, 0])
d2.scale.y = -1
coloRec(d2, bpy.data.materials['Deoxy'])

# make copies of bases, color
makeMolecule("Adenine1", "AttachA", location=[28, 48, 0])
new_a = makeMolecule("Adenine2", "AttachA", location=[28, 52, 0])
coloRec(new_a, bpy.data.materials['Adenine'])

makeMolecule("Cytosine1", "AttachC", location=[28, 56, 0])
new_c = makeMolecule("Cytosine2", "AttachC", location=[28, 60, 0])
coloRec(new_c, bpy.data.materials['Cytosine'])

makeMolecule("Guanine1", "AttachG", location=[28, 64, 0])
new_g = makeMolecule("Guanine2", "AttachG", location=[28, 68, 0])
coloRec(new_g, bpy.data.materials['Guanine'])

makeMolecule("Thymine1", "AttachT", location=[28, 72, 0])
new_t = makeMolecule("Thymine2", "AttachT", location=[28, 76, 0])
coloRec(new_t, bpy.data.materials['Thymine'])

makeMolecule("Uracil1", "AttachU", location=[28, 80, 0])
new_u = makeMolecule("Uracil2", "AttachU", location=[28, 84, 0])
coloRec(new_u, bpy.data.materials['Uracil'])

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

# -----------------------------
# TECHNIQUE 2: HELIX GENERATION
# -----------------------------

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

# the identity function
def identity(x):
    return x

# reverses a list, usually the reference sequence (not used in this project)
def reverse(reference):
    return reference[::-1]

# checks if seq is a valid sequence of nucleotides
# Aa, Cc, Gg represent the corresponding bases with deoxyribose and ribose sugars, respectively
# A: Adenine
# C: Cytosine
# G: Guanine
# T: Thymine
# U: Uracil
def isValidSeq(seq):
    return bool(re.match("^[AaCcGgTU]*$", seq))

# converts a DNA sequence to an RNA sequence
def DNAtoRNA(seq):
    return seq.replace("A", "a").replace("C", "c"
    ).replace("G", "g").replace("T", "U")

# adds a suffix to each term of a sequence
def addEnd(seq, suffix):
    s = [None] * len(seq)
    for i in range(0,len(seq)):
        s[i] = "{}{}".format(seq[i], suffix)
    return s

# Gets the complementary sequence of a reference sequence
# reference: a list of strings, each corresponding to the name of an object in 'objects'
# comp_fun: a function that takes an object name as input and returns the name of a complementary object
# object_names: a list of valid object names that can be optionally provided for validation
def getComp(reference, comp_fun, object_names = None):
    seq_length = len(reference)
    complement = [None] * seq_length
    for i in range(0, seq_length):
        complement[i] = comp_fun(reference[i])
        # check for non-present elements
        if object_names != None:    
            if not (reference[i] in object_names):
                print("Error: Elements of the reference sequence cannot be recognized.")
                return
            if not (complement[i] in object_names):
                print("Error: Elements of the complement sequence cannot be recognized.")
                return
    return complement

# Makes a collection of empties for a future helix, well-suited towards nucleotide biopolymers.
# collect_name: the name of the collection to make
# num_objects: the number of empties to make
# height: the height between each pair of bases
# base: angular rotation per base about z-axis, positive for right handed, negative for left handed
# location: the location of the center of the bottom rung, default is cursor location
def makeSkeleton(collect_name, num_objects, height, base = 0, location = None):
    # fallback location
    if location == None:
        location = bpy.context.scene.cursor.location
    # avoid making duplicate collections
    for collection in bpy.data.collections:
        if collect_name == collection.name:
            print("Error: A collection with the same name already exists.")
            return
    # create the new collection
    collection = bpy.data.collections.new(collect_name)
    bpy.context.scene.collection.children.link(collection)
    for i in range(0, num_objects):
        format_name = '{}_{}'.format(collection.name, i)
        rotate_name = 'Rotate_{}'.format(format_name)
        # reference empty
        o = bpy.data.objects.new("empty", None)
        bpy.data.collections[collect_name].objects.link(o)
        o.name = format_name
        o.location = location
        o.empty_display_size = height/4
        o.empty_display_type = 'PLAIN_AXES'
        # rotator empty
        r = bpy.data.objects.new("empty", None)
        bpy.data.collections[collect_name].objects.link(r)
        r.name = rotate_name
        r.location = [0,0,0]
        r.empty_display_size = height/8
        r.empty_display_type = 'CUBE'
        r.rotation_euler[2] = base*i
        r.parent = o
        location[2] += height

# Clears a collection of its helical decorations, handy for efficient reuse of skeleton
def clearChildren(collect_name):
    deselect()
    collection = bpy.data.collections[collect_name]
    # create object names
    object_names = []
    for i in range(0, len(collection.objects)):
        object_names.append(collection.objects[i].name)
    # select children of all top-level skeletal rungs
    i = 0
    while 'Rotate_{}_{}'.format(collection.name, i) in object_names:
        obj = collection.objects['Rotate_{}_{}'.format(collection.name, i)]
        seleRec(obj)
        obj.select_set(False)
        i += 1
    bpy.ops.object.delete()

# Decorates a skeleton with a single or double helix, well-suited towards nucleotide biopolymers.
# collection: the collection to be decorated with bases
# objects: a list of Blender objects that each rung of the helix is chosen from
# reference: a list of strings, each corresponding to the name of an object in 'objects'
# comp_fun: a function that takes an object name as input and returns the name of a complementary object
# strand: number of strands, 1 or 2
# radius: the distance of the center of each base from the z-axis
# tilt: for accuracy, tilt each base (reference upwards, complement downwards)
# disp: the angular displacement of each base of the complement strand relative to the reference strand
def makeHelix(collect_name, objects, reference, comp_fun, strand, radius, tilt, disp):
    if not (collect_name in bpy.data.collections):
        print("Error: No collection with this name exists.")
        return
    collection = bpy.data.collections[collect_name]
    seq_length = len(reference)
    # make a list of the names of all candidate objects
    object_names = []
    for i in range(0, len(objects)):
        object_names.append(objects[i].name)
    # check that all rungs exist
    for i in range(0, seq_length):
        rotate_name = 'Rotate_{}_{}'.format(collect_name, i)
        if not (rotate_name in collection.objects):
            print("Error: Skeleton length does not equal reference sequence length.")
            return
    # creates the complementary sequence
    complement = getComp(reference, comp_fun, object_names)
    # create the strands
    for i in range(0, seq_length):
        rotate_name = 'Rotate_{}_{}'.format(collect_name, i)
        if len(collection.objects[rotate_name].children) != 0:
            print("Error: Adornment has already occurred.")
            return
        # reference strand
        new_ref = copyRec(objects[reference[i]], None, collection)
        new_ref.name = '{}R'.format(rotate_name)
        new_ref.rotation_euler[1] = -tilt
        new_ref.rotation_euler[2] = -math.pi
        new_ref.location = [radius, 0, 0]
        new_ref.parent = collection.objects[rotate_name]
        # complement strand
        if strand == 2:
            new_com = copyRec(objects[complement[i]], None, collection)
            new_com.name = '{}C'.format(rotate_name)
            new_com.rotation_euler[0] = math.pi
            new_com.rotation_euler[1] = tilt
            new_com.rotation_euler[2] = disp-math.pi
            new_com.location = [radius*math.cos(disp), radius*math.sin(disp), 2*radius*math.sin(tilt)]
            new_com.parent = collection.objects[rotate_name]

# --------------------------
# TECHNIQUE 2: DEMONSTRATION
# --------------------------

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
# 3x vertical scaling to make viewing easier, 0.332 nm height/base pair
# right-handed, 10.5 bases per rotation
# 2 strands
# 10 angstrom (1 nm) in radius, solve r^2 + r^2 - 2^2 = 2*r*r*cos 135 to get 1.10
# 1.2 degree base tilt
# 135 degrees rotated from reference to complement
makeSkeleton("DNA1 Test", 24, 0.332*3, 2*math.pi/10.5, [48,48,0])
makeHelix("DNA1 Test", molCol().objects, DNA1, complementDNA, 2, 1.10, rad(1.2), rad(135))

# same as above, but coloring
makeSkeleton("DNA2 Test", 24, 0.332*3, 2*math.pi/10.5, [48,56,0])
makeHelix("DNA2 Test", molCol().objects, DNA2, complementDNA, 2, 1.10, rad(1.2), rad(135))

# RNA corresponding to DNA
makeSkeleton("RNA1 Test", 24, 0.332*3, 2*math.pi/10.5, [48,64,0])
makeHelix("RNA1 Test", molCol().objects, RNA1, complementRNA, 1, 1.10, rad(1.2), rad(135))

# unraveled
makeSkeleton("RNA2 Test", 24, 0.332*3, 0, [48,72,0])
makeHelix("RNA2 Test", molCol().objects, RNA2, complementRNA, 1, 1.10, rad(1.2), rad(135))

# left-handed DNA is not realistic, since the backbone doesn't line up
makeSkeleton("DNA3 Test", 24, 0.332*3, -2*math.pi/10.5, [48,80,0])
makeHelix("DNA3 Test", molCol().objects, DNA1, complementDNA, 1, 1.10, 0, rad(135))

# unrealistic DNA - extreme twisting
makeSkeleton("DNA4 Test", 24, 0.332*3, 2*math.pi/10.5, [48,88,0])
makeHelix("DNA4 Test", molCol().objects, DNA2, complementDNA, 2, 1.10, rad(1.2*20), rad(135))

# very basic, exaggerated double-stranded left-handed RNA model for fun
makeSkeleton("RNA3 Test", 48, 1.2, -2*math.pi/10.5/3, [48,104,0])
makeHelix("RNA3 Test", molCol().objects, 2*RNA3, complementRNA, 2, 1.10*3, rad(1.2*20), rad(90))

# double-stranded basic DNA, but now every base is paired with itself and no major-minor is present
makeSkeleton("DNA5 Test", 24, 1.1, 2*math.pi/10.5, [48,120,0])
makeHelix("DNA5 Test", molCol().objects, DNA3, identity, 2,  1.10, rad(1.2), math.pi)

# clears the skeleton of t1 as a test
#clearChildren("DNA1 Test")
# restore the helix
#makeHelix("DNA1 Test", molCol().objects, DNA1, complementDNA, 2, 1.10, rad(1.2), rad(135))

# --------------------------
# TECHNIQUE 3: BOND ROTATION
# --------------------------

# sets the interpolation method for an object
def setInterpolation(obj, inter="LINEAR"):
    for fcurve in obj.animation_data.action.fcurves:
        for keyframe in fcurve.keyframe_points:
            keyframe.interpolation = inter

# checks if two object names are 'equal' in the context of a helix
# the reason why we split by . is that Blender adds dots to duplicates
def checkNameEqual(name1, name2):
    return name1.split('.')[0] == name2.split('.')[0]

# given an object (usually a bond) that is part of a rung
# select the corresponding object in all other bonds
def selectSim(object = None):
    if object == None:
        object = bpy.context.selected_objects[0]
    collection = object.users_collection[0]
    for obj in collection.objects:
        if checkNameEqual(obj.name, object.name):
            obj.select_set(True)

# obj: the attachment or bond to be rotated about
# n: the number of keyframes to add
# vf: the number of frames between each rotation
# vt: the angle of rotation per frame, obeying right-hand rule along bond
# f: the initial frame number
# inter: the interpolation type (CONSTANT, LINEAR, BEZIER) 
def rotateBond(obj, n, vf, vt, f = 1, inter='LINEAR', subEmpty=True):
    # if you select a bond, make sure you go to the descendant!
    if subEmpty and obj.type == 'EMPTY':
        obj = obj.children[0]
    ref_z = obj.rotation_euler.z
    for i in range(0, n):
        bpy.context.scene.frame_set(f)
        obj.rotation_euler.z = ref_z
        obj.keyframe_insert(data_path = "rotation_euler", index = -1)
        f += vf
        ref_z += vt
    # change all interpolation types
    setInterpolation(obj, inter)

# ---------------------------
# TECHNIQUE 4: HELIX ROTATION
# ---------------------------

# resets the animations of an object, leaving parameters as found at the frame
def resetAnimation(obj, default_frame = 1):
    bpy.context.scene.frame_set(default_frame)
    obj.animation_data_clear()
    
# resets the animations of a group of objects
def resetAnimations(objects, default_frame = 1):
    for obj in objects:
        resetAnimation(obj, default_frame)

# assigns animations to a singular object
# Note: Not used in this project
# f: the starting frame
# vf: the distance between each frame
# pos_rot: stands for positions, rotations
#   note: of the form [[x,y,z], [rot_x, rot_y, rot_z]] * n
def animateObject(obj, f, vf, pos_rot):
    for i in range(0, len(frames_pr)):
        bpy.context.scene.frame_set(f)
        obj.location = pos_rot[i][0]
        obj.rotation_euler = pos_rot[i][1]
        obj.keyframe_insert(data_path = "location", index = -1)
        obj.keyframe_insert(data_path = "rotation_euler", index = -1)
        f += vf

# given a single selected bond, rotates all similar bonds in the same way    
# n: the number of keyframes to add
# vf: the number of frames between each rotation
# vt: the angle of rotation per frame, obeying right-hand rule along bond
# f: the initial frame number
# t: the initial angle theta of rotation
# inter: the interpolation type (CONSTANT, LINEAR, BEZIER) 
def rotateSimilars(n, vf, vt, f = 1, inter='LINEAR', clear=False):
    if len(bpy.context.selected_objects) != 1:
        print("Error: multiple bonds selected.")
        return
    selectSim(bpy.context.selected_objects[0])
    if clear:
        resetAnimations(bpy.context.selected_objects)
    for obj in bpy.context.selected_objects:
        rotateBond(obj, n, vf, vt, f, inter, subEmpty=True)

# assigns coiling rotation animations to a collection of objects comprising a helix
# collect_name: name of the collection
# n: the number of keyframes to add, starting at f and increasing by vf
# vf: the number of frames between each rotation
# vt: the angle of rotation per frame, obeying right-hand rule along bond
# f: the initial frame number
# t: the initial angle theta of rotation
def coilCollection(collect_name, n, vf, vt, f = 1, inter="LINEAR", uncoil=True):
    # get all objects in the collection
    objects = bpy.data.collections[collect_name].objects
    # get object names
    object_names = []
    for i in range(0, len(objects)):
        object_names.append(objects[i].name)
    # count the number of rungs in the collection   
    num_rungs = 0
    while '{}_{}'.format(collect_name, num_rungs) in object_names:
        num_rungs += 1
    for i in range(0, num_rungs):
        obj = objects[object_names.index('Rotate_{}_{}'.format(collect_name, i))]
        resetAnimation(obj, f)
        rotateBond(obj, n, vf, vt*i, f = f, inter=inter, subEmpty=False)
    if uncoil:
        for i in range(0, num_rungs):
            obj = objects[object_names.index('Rotate_{}_{}'.format(collect_name, i))]
            rotateBond(obj, n, vf, -vt*i, f = f+(n-1)*vf, inter=inter, subEmpty=False)

# -----------------------------
# TECHNIQUES 3, 4 DEMONSTRATION
# -----------------------------

# Example 1: Rotating / Uncoiling DNA that also follows a path
if False: # prevent accidental activation when sourcing
    coilCollection("DNA1 Test", 5, 30, -math.pi/8)
    collect_name = "DNA1 Test"
    objects = bpy.data.collections[collect_name].objects
    for i in range(0, 24):
        obj = objects['{}_{}'.format(collect_name, i)]
        resetAnimation(obj)
        for j in range(0, 21):
            r = 15
            angle = math.pi/60 * i + math.pi/10 * j
            bpy.context.scene.frame_set(12*j)
            obj.location = [r*math.cos(angle), 0, r*math.sin(angle)]
            obj.keyframe_insert(data_path = "location", index = -1)
            obj.rotation_euler.y = angle
            obj.keyframe_insert(data_path = "rotation_euler", index = -1)
        
# Example 2: DNA / RNA Rotating in Place
if False: # prevent accidental activation when sourcing
    coilCollection("DNA2 Test", 5, 30, -math.pi/8)


# Example 3
if False: # prevent accidental activation when sourcing
coilCollection("RNA1 Test", 2, 120, -2*math.pi/10.5)
# note: need to select a bond pointing out of phosphate (purple sphere) in RNA Test 1
rotateSimilars(2, 120, -math.pi/6, f = 1, clear=True)
rotateSimilars(2, 120, math.pi/6, f = 121)
collect_name = "RNA1 Test"
objects = bpy.data.collections[collect_name].objects
for i in range(0, 24):
    obj = objects['{}_{}'.format(collect_name, i)]
    resetAnimation(obj, 0)
    ref_z = obj.location.z
    for j in range(0, 10):
        bpy.context.scene.frame_set(12*j)
        obj.location.z = ref_z
        obj.keyframe_insert(data_path = "location", index = -1)
        ref_z += 1
    for j in range(10, 21):
        bpy.context.scene.frame_set(12*j)
        obj.location.z = ref_z
        obj.keyframe_insert(data_path = "location", index = -1)
        ref_z -= 1
        
# Example 4
# select some random bonds in the duplicate adenine and make random movements
if False: # prevent accidental activation when sourcing
    rotateBond(bpy.context.selected_objects[0], 9, 30, math.pi)
    resetAnimation(bpy.context.selected_objects[0])
    bpy.context.selected_objects[0].animation_data_clear()

# NOT USED: animateObject
#collect_name = "DNA1 Test"
#objects = bpy.data.collections[collect_name].objects
#animate_data_1 = [None] * 24
#vf = -math.pi/6
#h = 2

#object_names = []
#for i in range(0, len(objects)):
#    object_names.append(objects[i].name)

#for rung_num in range(0, 24):
#    animate_data_1[rung_num] = [None] * 11
#    obj = objects[object_names.index('{}_{}'.format(collect_name, rung_num))]
#    print(obj)
#    rot_start = obj.rotation_euler.copy()
#    pos_start = obj.location.copy()
#    print(rot_start)
#    print(pos_start)
#    print("---")
#    for keyframe in range(0, 6): # 6 fps for 10 s
#        animate_data_1[rung_num][keyframe] = [30*keyframe, 
#        [20*math.cos(math.pi/12*keyframe), 
#        20*math.sin(math.pi/12*keyframe), 
#        pos_start[2]], rot_start.copy()]
#        pos_start[2] += h
#        rot_start[2] += vf*rung_num
#    for keyframe in range(6, 11):
#        animate_data_1[rung_num][keyframe] = [30*keyframe,
#        [20*math.cos(math.pi/12*keyframe), 
#        20*math.sin(math.pi/12*keyframe), 
#        pos_start[2]], rot_start.copy()]
#        pos_start[2] -= h
#        rot_start[2] -= vf*rung_num

#animateCollection(collect_name, animate_data_1)

# ----------------
# BONUS: INTERFACE
# ----------------

# please see hideEmpties
class HideEmpties(Operator):
    bl_idname = "curve.hide_empties"
    bl_label = "HideEmpties"
    hide: BoolProperty(
        name="Hide",
        description="Whether empties will be hidden by this",
        default=True
    )
    def execute(self, context):
        hideEmpties(self.hide)
        return {'FINISHED'}

# please see seleRec, copyRec, joinRec
class OperatorRec(Operator):
    bl_idname = "curve.operator_rec"
    bl_label = "OperatorRec"
    type: FloatProperty(
        name="Type",
        description="select=1, copy=2, join=3",
        default=1
    )
    def execute(self, context):
        if len(bpy.context.selected_objects) != 1:
            print("Error: number of selected objects does not equal 1.")
            return {'FINISHED'}
        obj = bpy.context.selected_objects[0]
        if self.type == 1:
            seleRec(obj)
        if self.type == 2:
            copyRec(obj)
        if self.type == 3:
            joinRec(obj)
        return {'FINISHED'}

# please see resetAnimations
class ResetAnimations(Operator):
    bl_idname = "curve.reset_animations"
    bl_label = "ResetAnimations"
    def execute(self, context):
        params = context.window_manager.chang_final_params
        resetAnimations(bpy.context.selected_objects, params.default_frame)
        return {'FINISHED'}

# please see rotateBond
class RotateBond(Operator):
    bl_idname = "curve.rotate_bond"
    bl_label = "RotateBond"
    subEmpty: BoolProperty(
        name="subEmpty",
        description="Whether empties will have their children rotated.",
        default=True
    )
    def execute(self, context):
        params = context.window_manager.chang_final_params
        if len(bpy.context.selected_objects) != 1:
            print("Error: number of selected objects does not equal 1.")
            return {'FINISHED'}
        obj = bpy.context.selected_objects[0]
        rotateBond(obj, params.n4, params.vf4, params.vt4, params.f4, "LINEAR", self.subEmpty)
        return {'FINISHED'}

class CURVE_PT_ChangFinalPanel(Panel):
    bl_label = "Nucleic Acid Techniques"
    bl_idname = "CURVE_PT_ChangFinalPanel"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Create"
    bl_context = "objectmode"
    bl_options = {"DEFAULT_CLOSED"}
    def draw(self, context):
        layout = self.layout
        wm = context.window_manager
        col = layout.column(align=True)
        col.label(text="Technique Set 1")
        prop_hide = col.operator("curve.hide_empties", text="Hide Empties")
        prop_hide.hide = True
        prop_show = col.operator("curve.hide_empties", text="Show Empties")
        prop_show.hide = False
        prop_sele = col.operator("curve.operator_rec", text="Recursive Select")
        prop_sele.type = 1
        prop_copy = col.operator("curve.operator_rec", text="Recursive Copy")
        prop_copy.type = 2
        prop_join = col.operator("curve.operator_rec", text="Recursive Join")
        prop_join.type = 3
        col = layout.column(align=True)
        col.label(text="Technique Set 2")
        col = layout.column(align=True)
        col.label(text="Technique Set 3")
        col = layout.column(align=True)
        col.label(text="Technique Set 4")
        prop_reset = col.operator("curve.reset_animations", text="Reset Animations")
        col.prop(wm.chang_final_params, "default_frame")
        col = layout.column(align=True)
        prop_sube = col.operator("curve.rotate_bond", text="Rotate Bond")
        prop_sube.subEmpty=True
        prop_subn = col.operator("curve.rotate_bond", text="Rotate Group")
        prop_subn.subEmpty=False
        col.prop(wm.chang_final_params, "n4")
        col.prop(wm.chang_final_params, "vf4")
        col.prop(wm.chang_final_params, "vt4")
        col.prop(wm.chang_final_params, "f4")

class ChangFinalParams(PropertyGroup):
    default_frame: FloatProperty(
        name="Default Frame",
        default=1,
        min=0
    )
    n4: FloatProperty(
        name="Number of Intervals",
        default=2,
        min=0
    )
    vf4: FloatProperty(
        name="Frames per Interval",
        default=30,
        min=0
    )
    vt4: FloatProperty(
        name="Angular Velocity",
        default=0
    )
    f4: FloatProperty(
        name="Starting Frame",
        default=1
    )

my_classes = (
    HideEmpties,
    OperatorRec,
    ResetAnimations,
    RotateBond,
    # parameter holder
    ChangFinalParams,
    # panel
    CURVE_PT_ChangFinalPanel
)

# registers classes, a general method
def register(classes):
    # add all classes
    for cls in classes:
        bpy.utils.register_class(cls)
    # create a parameter holder in WindowManager
    bpy.types.WindowManager.chang_final_params = PointerProperty(type=ChangFinalParams)

# unregisters classes, a general method
def unregister(classes):
    # remove parameters from WindowManager
    del bpy.types.WindowManager.chang_final_params
    # remove all classes
    for cls in my_classes:
        bpy.utils.unregister_class(cls)
     
unregister(my_classes)
register(my_classes)