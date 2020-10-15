import numpy as np
import math
import re

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

# collect_name: the name of the collection to be created with these objects
# objects: a list of Blender objects that each rung of the helix is chosen from
#   not a function parameter, is obtained from all user-selected objects
# reference: a list of strings, each corresponding to the name of an object in 'objects'
# comp_fun: a function that takes an object name as input and returns the name of a complementary object 
# rev: true if the sequence should be reversed in advance
# comp: true if the sequence should be complemented in advance
# strand: number of strands, 1 or 2
# angle: positive for right hand rule, negative for left hand rule, determines rotation per base
# radius: the distance of the center of each base from the z-axis
# tilt: to create the major / minor grooves, each base must be tilted
def makeHelix(collect_name, reference, comp_fun, comp, rev, strand, angle, radius, tilt, rotate):
    objects = bpy.context.selected_objects
    collection = bpy.data.collections.new(collect_name)
    bpy.context.scene.collection.children.link(collection)
    object_names = []
    for i in range(0, len(objects)):
        object_names.append(objects[i].name)
    print(object_names)
    if rev:
        reference = reference[::-1]       
    complement = []
    for i in range(0, len(reference)):
        complement.append(comp_fun(reference[i])) 
    if comp:
        temp = reference
        reference = complement
        complement = temp       
    commonHeight = 0
    for i in range(0, len(reference)):
        ref_name = reference[i]
        format_name = 'Reference_{}_{}'.format(i, ref_name)
        ref_obj = objects[object_names.index(ref_name)]
        new_obj = bpy.data.objects.new(format_name, ref_obj.data)
        new_obj.location = [radius*math.cos(angle*i), radius*math.sin(angle*i), commonHeight]
        commonHeight += new_obj.dimensions[2] / math.cos(tilt)
        new_obj.rotation_euler[1] = -tilt
        new_obj.rotation_euler[2] = angle*i - math.pi
        bpy.data.collections[collect_name].objects.link(new_obj)
    if strand > 1:
        commonHeight = 2*radius*math.sin(tilt)
        for i in range(0, len(reference)):
            com_name = complement[i]
            format_name = 'Complement_{}_{}'.format(i, com_name)
            ref_obj = objects[object_names.index(com_name)]
            new_obj = bpy.data.objects.new(format_name, ref_obj.data)
            new_obj.location = [radius*math.cos(angle*i + rotate), 
                radius*math.sin(angle*i + rotate), commonHeight]
            commonHeight += new_obj.dimensions[2] / math.cos(tilt)
            new_obj.rotation_euler[1] = tilt
            new_obj.rotation_euler[2] = angle*i + rotate - math.pi
            bpy.data.collections[collect_name].objects.link(new_obj)

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

DNAseq = "ACCGATCGCGGGGAGGGGTCAAGCGCGACGTTCGCCACGCAAAACCGGGCACCCCGACTAACCGTACCGGCGCTACGCCCATTTACGAACGTAAGCGGTA"
DNA_list = list(DNAseq)

# pretty accurate
makeHelix("Realistic DNA", DNA_list, complementDNA, False, False, 2, 2*math.pi/10.5, 3, math.pi/180*1.2, math.pi)

# cool looking, not necessarily accurate
makeHelix("Stylish DNA", DNA_list, complementDNA, False, False, 2, -2*math.pi/10.5/3, 3*3, math.pi/180*1.2*18, math.pi)

# cool looking and reflects reality to a degree
makeHelix("Stylish DNA 2", DNA_list, complementDNA, 
False, False, 2, -2*math.pi/10.5/3, 3*3, math.pi/180*1.2, math.pi*2/3)

RNAseq = convertType(DNAseq, "RNA")
RNA_list = list(RNAseq)

# cool looking, not necessarily accurate
makeHelix("Stylish RNA", RNA_list, complementRNA, False, False, 1, -2*math.pi/10.5/3, 3*3, math.pi/180*1.2*18)
makeHelix("Stylish RNA 2", RNA_list, complementRNA, False, False, 1, -2*math.pi/10.5/3, 3*3, 0)

#    
# = []

tetConstant = 109.47 * math.pi / 180
triConstant = 120.00 * math.pi / 180

objects = bpy.context.selected_objects
object_names = []
for i in range(0, len(objects)):
    object_names.append(objects[i].name)
    

phosphate_attachments = [
    ["Oxygen", [0, tetConstant, 0]],
    ["Oxygen", [0, tetConstant, triConstant]],
    ["Oxygen", [0, tetConstant, -triConstant]]
]

impossible_attachments = [
    ["Phosphate", [0, tetConstant, 0]],
    ["Phosphate", [0, tetConstant, triConstant]],
    ["Phosphate", [0, tetConstant, -triConstant]]
]

parent_obj = bpy.context.selected_objects[0]

for i in range(0, len(phosphate_attachments)):
    o = bpy.data.objects.new( "empty", None )
    ref_name = phosphate_attachments[i][0]
    format_name = 'Phosphate_{}_{}'.format(i, ref_name)
    ref_obj = objects[object_names.index(ref_name)]
    new_obj = bpy.data.objects.new(format_name, ref_obj.data)
    new_obj.rotation_euler = [math.pi,0,0]
    new_obj.parent = o
    o.parent = parent_obj
    new_obj.location = (0, 0, 5)
    bpy.data.collections["Collection"].objects.link(new_obj)
    o.rotation_euler = phosphate_attachments[i][1]
    
for i in range(0, len(impossible_attachments)):
    o = bpy.data.objects.new( "empty", None )
    ref_name = impossible_attachments[i][0]
    format_name = 'Impossible_{}_{}'.format(i, ref_name)
    ref_obj = objects[object_names.index(ref_name)]
    new_obj = bpy.data.objects.new(format_name, ref_obj.data)
    new_obj.rotation_euler = [math.pi,0,0]
    new_obj.parent = o
    o.parent = parent_obj
    new_obj.location = (0, 0, 5)
    bpy.data.collections["Collection"].objects.link(new_obj)
    o.rotation_euler = impossible_attachments[i][1]
