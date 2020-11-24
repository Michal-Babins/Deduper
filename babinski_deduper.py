#!/usr/bin/env python


import gzip
import re
import argparse
import pysam
#Argparse takes in arguments that are callable when running the code.

def args():
    ''' Argparse takes various arguments that will be specific to the users input. Argparse requires
        a sam file, and umi file. '''
    parser=argparse.ArgumentParser(description = "Deduper code")
    parser.add_argument("-i","--input", help="sam file", required = True, type = str)
    parser.add_argument("-u","--umi", help="File containing Umi's", required = True, type = str)
    parser.add_argument("-p","--paired", help="Paired end files not supported", required = False, type = str)
    return parser.parse_args()
args = args()

#Setting the argparse argurment 
input_file = args.input
known_umi_file = args.umi
paired_end = args.paired

if paired_end == "Paired" or paired_end == "paired":
    print("Paired-End feature not included in package, exiting")
    exit()

#Sort sam file to be in order
pysam.sort("-o" , input_file, "-O" , "sam" , input_file)

#Create a dictionary to store umi information in
umi_dict = {}

#Open umit file and get known umi's from the list
with open(known_umi_file, "r") as index_file:
    for line in index_file:
        line = line.strip()
        umi_dict[line] = () #Fill dictionary with Umi's as key
index_file.close()

# deduped = open('{}.deduped'.format(input_file),'w')

with open(input_file, "r") as sam:
    '''Starting chromosome information needs to be stored for downstream comparison of umi and postions
        to determine if a file should be written to the deduped file or passed. '''
    for string in sam:
        if not string.startswith("@") :
            header = string.split("\t")
            chromosome = header[1]
            break
sam.close()

#Initialize dictionary that will store keys:True Position and values:Umi's
iter_pos_dict = {}

#Open file for final deduped file to be placed into
deduped = open('{}.deduped'.format(input_file),'w')

#Initiate Direction value
direction = 0

with open(input_file, "r") as sam:
    '''Main Chunk. The purpose is to grab all the relevant iniformation from the header, adjust 
    for directionality and store true position and umi in a dictionary needed for comparison 
    in pcr deduping. Chromosomes will be iterated through and reset to move through the reads in 
    same file.'''
    for string in sam:
        if not string.startswith("@"):
            #splittig the header to grab each section
            header_string = string.split("\t")
            #specifically spliting first groupt to grab the UMI
            umi_index = header_string[0].split(":")
            umi_key = umi_index[7] #Store UMI
            chrom = header_string[1] #Store Chromosome
            tru_pos = int(header_string[3]) #Store postion (will later be adjusted)
            flag = int(header_string[4]) #Store bitwise
            cigar = header_string[5] #Store cigar string

            #Iterate through cigar string to grab numerical sets to configure true postion
            cig = re.search(r'([0-9]+)([MIDNSHPX=])',cigar)[0]

            #Make sure to reset dictioinary on every new chromosome
            if chrom != chromosome:
                iter_pos_dict = {}
            #Set direciton of reads
            if umi_key in umi_dict:
                if ((flag & 16) == 16):
                    direction = 1 #Reverse Direction
                else:
                    direction = 0 #Forward Direction
            #Adjust for position for each read, main worry is if soft clipping occurred
                if direction == 0:
                    if "S" in cig:
                        soft_clip_forward = cigar.split("S")[0] #Grab numierical value before that first soft_clip letter
                        tru_pos = tru_pos - int(soft_clip_forward)
                    else:
                        tru_pos = tru_pos
                #Adjust positon in cigar string where soft clipping occurs, after first S
                elif direction == 1: 
                    if "S" in cigar:
                        #Assumption made soft clipping will not be higher than a 3 digit number
                        for i in range(3):
                            if cigar[i] == "S":
                                #adjust to after first soft clipping occurence (accounted for in forward direction)
                                cigar = cigar[i+1:]
                    #Adjust to disclude I from final summary 
                    if "I" in cigar:
                        I = re.findall(r'([0-9]+[I])',cigar) #Grab only #I
                        I = I[0].split("I")[0] #Split for numerical value associated with I
                        all_numerical_values = re.findall(r"([0-9]+)", cigar) #grab all num
                        iterated_numbers = [int(i) for i in all_numerical_values]
                        tru_pos = tru_pos + sum(iterated_numbers) - int(I) #subtract I from final sum
                    else:                
                        #If those identifiers arent there, then only get numerical values to add to final sum
                        all_numerical_values = re.findall(r"([0-9]+)",cigar)
                        iterated_numbers = [int(i) for i in all_numerical_values]
                        tru_pos = tru_pos + sum(iterated_numbers)


                #If location has an exact match than pass
                if tru_pos in iter_pos_dict:
                    pass 
                #If location not repeat, read = deduped, write file. 
                else:
                    iter_pos_dict[tru_pos] = (umi_key)
                    deduped.write(string)
#Close file
sam.close()
deduped.close()

