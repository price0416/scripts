#!/usr/bin/python
import sys, os
import csv
import math
import StringIO
from optparse import OptionParser
from Bio import SeqIO
#import string

#def get_protein_id(description):
 #   first_cut_pos = string.find(description,"protein_id")
  #  first_trim = description[first_cut_pos:]
   # second_cut_pos = string.find(first_trim,"=")
    #second_trim = first_trim[second_cut_pos+1:]
#    third_cut_pos = string.find(second_trim, "]")
 #   third_cut = second_trim[:third_cut_pos]
  #  return third_cut

#def get_protein_location(description):
 #   first_cut_pos = string.find(description,"location")
  #  first_trim = description[first_cut_pos:]
   # second_cut_pos = string.find(first_trim,"=")
#    second_trim = first_trim[second_cut_pos+1:]
 #   third_cut_pos = string.find(second_trim, ".")
  #  third_cut = second_trim[:third_cut_pos]
   # return third_cut


if __name__ == '__main__':

    usage = "\n\nParses read count file into separate files representing coding regions based on genbank annotations.\nUsage: %prog arg"
    parser = OptionParser(usage)

    parser.add_option("-a", "--anno", dest="anno_file", help="Genbank annotation file (.gbk)")
    parser.add_option("-c", "--counts", dest="count_file", help="File representing # reads for each position. (.csv)")
    parser.add_option("-o", "--output", dest="output_file", help="Output file.")

    (options, args) = parser.parse_args()
    
    if len(args) != 0  or not options.anno_file or not options.count_file or not options.output_file:
                parser.error("incorrect number of arguments.\n\t\tUse -h to get more information")

    genome=SeqIO.read(options.anno_file, 'genbank')
    #fasta_file = open(options.fasta_file, "rU")
    #fasta_genes = list(SeqIO.parse(fasta_file,"fasta")
    num_genes = len(genome.features)-1
    gene_start = []
    gene_end   = []
    gene_name  = []


    try:
        countfile = open(options.count_file, "r")
        try:
            count_string = countfile.read()
        finally:
            countfile.close()
    except IOError:
        pass

    li = count_string.split(',')[:-1]
    count_list = [int(x) for x in li]

    count_len = len(count_list)
    print "There are " + str(count_len) + " items in the count list."

    outfile = open(options.output_file, "w")

    for i,feature in enumerate(genome.features):
        if feature.type != "CDS" and feature.type != "rRNA":
            continue
        gene_start.append(int(feature.location.start))
        gene_end.append(int(feature.location.end))
        gene_name.append(feature.qualifiers['locus_tag'][0])


    for cur_gene in range(len(gene_start)):
        startval = gene_start[cur_gene]
        if startval == 1:
            continue
        endval   = gene_end[cur_gene]
        value_list = count_list[startval:endval]

#        outfile.write("Gene_"+str(startval) + '\t')
        outfile.write(gene_name[cur_gene] + '\t')
        for x in enumerate(value_list[:-1]):
            outfile.write(str(x[1]) + ';')

        outfile.write(str(value_list[-1]));
        outfile.write('\n')

    outfile.close()

