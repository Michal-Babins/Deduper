#!/usr/bin/env python



import cairo
import math
import re
import random
import numpy as np
import argparse

#Argparse takes in arguments that are callable when running the code.

def args():
    ''' Argparse takes various arguments that will be specific to the users input. Argparse requires
        a fasta file, and motif file. '''
    parser=argparse.ArgumentParser(description = "visualize your motifs!")
    parser.add_argument("-i","--input", help="fasta file with sequences", required = True, type = str)
    parser.add_argument("-m","--motif", help="File containing motifs", required = True, type = str)
    return parser.parse_args()
args = args()

#set argparse 
input_file = args.input
motif_txt = args.motif


''' Dictionaty containing regex expression for IUPAC values'''
IUPAC = {
    "A":"[Aa]",
    "C":"[Cc]",
    "G":"[Gg]",
    "T":"[TtUu]",
    "U":"[UuTt]",
    "W":"[AaTtUu]",
    "S":"[CcGg]",
    "M":"[AaCc]",
    "K":"[GgTtUu]",
    "R":"[AaGg]",
    "Y":"[CcTtUu]",
    "B":"[CcGgTt]",
    "D":"[AaGgTtUu]",
    "H":"[AaCcTtUu]",
    "V":"[AaCcGg]",
    "N":"[AaCcGgTtUu]",
    "Z":"[-]",
}




def parse_motif(motif_txt):
    ''' Extract motifs from motif file '''
    with open(motif_txt, "r") as motif:
        motif_list = [] 
        for line in motif:
            line = line.strip() 
            motif_list.append(line)

    return motif_list


def get_exon():
    ''' Iterate through seqeunces to find exons in seqs'''
    exon_dict = {}
    for i in my_dict.keys():
        exon_seq = my_dict[i] #match exon sequences in main seq dictionary 
        iterate = re.finditer("[A-Z]+", exon_seq) 
        for match in iterate:
            exon_dict[i] = match.span() #finding start and stop place of exon

    return exon_dict



def get_motif(motif_list):
    ''' Convert motifs to possible IUPAC values'''
    mtfL = []
    for i in motif_list:
        new_motif = "" #set empty string
        i = i.upper() #convert all to upper
        for char in i:
            new_motif += str(IUPAC[char]) #add mtf values in IUPAC terms before appendign to list
        mtfL.append(new_motif)
    return(mtfL)



def long_gene():
    ''' Find the longest gene value for reference to plot width  '''
    long_johns = 0
    for i in my_dict:
        val = my_dict[i]
        if len(val) > long_johns:
            long_johns = len(val)

    return long_johns



my_dict = {} #set main dict
''' Retrieving header and sequence information from the fasta '''
with open(input_file, "r") as fa:
    my_dict = {}
    for line in fa:
        line = line.strip()

        if line[0] == '>':
            header = line
            seq = ''
        elif line[0] != '>':
            seq += line
            my_dict[header] = (seq) 


#set the longest gene length 
width_value = long_gene()

#set width and heigth values
width = int(width_value)*1.25
height = int(len(my_dict.values()) * 120) + 20

surface = cairo.SVGSurface("motif_plot.svg", width, height)
context = cairo.Context(surface)

#set vertical and horizontal positons 
cord_start = 25
position_vert = 25

#extract exons
exon_d = get_exon()

#extract motifs
motif_list = parse_motif(motif_txt)
mtf = get_motif(motif_list)


def col():
    ''' set values for rgb input'''
    colors = {}
    count = 0
    for i in range(50):
        count += 1
        x,y,z = random.random(),random.random(), random.random()
        colors[count] = [x,y,z]
    return colors

#get color dictionary 
col = col()



#find the genes im looking for
for i in my_dict:
    gene_seq = my_dict[i] #set gene_seqs
    #get physical, linear placement of the genes
    context.set_line_width(5)
    context.set_source_rgb(0,0,0)
    context.move_to(cord_start,position_vert +25)
    context.line_to(cord_start+len(gene_seq), position_vert + 25) #plot line of gene length
    context.stroke()

    #place exon span as rectangle
    exon_cord = exon_d[i] #get exon pos
    exon_cord = list(exon_cord)
    context.rectangle(exon_cord[0], position_vert, exon_cord[1] - exon_cord[0],50) #place exon as large rectangle
    context.set_source_rgb(0.8,0.8,0.8) #set gray scale
    context.fill()

    #find and plot motifs
    motif_count = 1 #start motif counter
    for j in mtf:
        yutter = re.finditer(j,gene_seq) #find mtf in seq
        for matchseq in yutter:
            mtf_finder = matchseq.span() #gather location of motif
            mtf_finder = list(mtf_finder) #turn into accessible list

            x,y,z = col[motif_count] #extract values for rgb 

            context.set_line_width(5)
            context.set_source_rgb(x,y,z) #insert rbg values
            context.rectangle(mtf_finder[0] + 25, position_vert,mtf_finder[1] - mtf_finder[0], 50) #place vertical motif rectangle
            context.fill()
        motif_count += 1 #increase motif counter

    #gene name
    context.move_to(cord_start + 10, (position_vert) - 5) #place gene pos
    context.set_source_rgb(0,0,0) #set color to black
    context.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    context.set_font_size(12)
    context.show_text(i) #show gene name
    context.stroke()


    position_vert += 100


#Plot legend
context.move_to(25, height -100)
context.set_font_size(16)
context.show_text("Legend:")


#get visual key for the motifs
motif_count = 1
lenlg = 415
for mot in motif_list:

    context.set_source_rgb(0,0,0)
    context.move_to(25, lenlg + 5)
    context.show_text(mot) #place motif
    context.set_line_width(1) #draw line between motifs in key
    context.line_to(width - 567, lenlg + 5)
    context.stroke()

    x,y,z = col[motif_count] #extract rgb equivalents in 

    context.set_source_rgb(x,y,z) #set legend colors
    context.rectangle(width - 580, lenlg - 10 , 10, 10) #show motif color
    context.fill()

    motif_count += 1
    lenlg += 22

surface.finish() #fin
