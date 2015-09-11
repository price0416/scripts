#!/usr/bin/python
import sys, os
import csv
import math
import StringIO
from optparse import OptionParser

if __name__ == '__main__':

    usage = "\n\nAccepts 2 .csv files represneting read counts for each position.\nUsage: %prog arg"
    parser = OptionParser(usage)

    parser.add_option("-1", "--file1", dest="file1", help="Count file #1")
    parser.add_option("-2", "--dile2", dest="file2", help="Count file #2")
    parser.add_option("-o", "--output", dest="output_file", help="Output file destination")

    (options, args) = parser.parse_args()
    
    if len(args) != 0  or not options.file1 or not options.file2 or not options.output_file:
                parser.error("incorrect number of arguments.\n\t\tUse -h to get more information")


    f1_string = ""
    f1_string = ""
    
    try:
        f1 = open(options.file1, "r")
        f2 = open(options.file2, "r")
        try:
            f1_string = f1.read()
            f2_string = f2.read() 
        finally:
            f1.close()
            f2.close()
    except IOError:
        pass

    li = f1_string.split(',')[:-1]
    f1_list = [int(x) for x in li]
    li = f2_string.split(',')[:-1]
    f2_list = [int(x) for x in li]

    f1_len = len(f1_list)
    f2_len = len(f2_list)

    add2length = max(len(f1_list),len(f2_list)) - min(len(f1_list),len(f2_list))
    if f1_len > f2_len:
        print "f1 bigger"
        f2_list += [0] * add2length
    if f1_len < f2_len:
        print "f2 bigger"
        f1_list += [0] * add2length

    print str(len(f1_list)) + " items in f1_list."
    print str(len(f2_list)) + " items in f2_list."

    f1_numReads = 0
    f2_numReads = 0
    for position in range(len(f1_list)):
        f1_numReads = f1_numReads + f1_list[position]
        f2_numReads = f2_numReads + f2_list[position]

    print "There are " + str(f1_numReads) + " reads in f1 file."
    print "There are " + str(f2_numReads) + " reads in f2 file."

    print "Genome consists of " + str(len(f1_list)) + " positions."

    outfile = open(options.output_file, "w")
    output_list = []

    for position in range(len(f1_list)):
        sum_value = 0
        sum_value = f1_list[position] + f2_list[position]
        outfile.write(str(sum_value));
        output_list.append(sum_value)
        
        if position == len(f1_list)-1:
            break;
        else:
            outfile.write(",");
            outfile.write("\n");
        
    outfile.close()

    out_numReads = 0
    for position in range(len(f1_list)):
        out_numReads = out_numReads + output_list[position]

    print "There are " + str(out_numReads) + " in the combined file." 
