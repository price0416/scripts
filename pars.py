#!/usr/bin/python
import sys, os
import csv
import math
import StringIO
from optparse import OptionParser

if __name__ == '__main__':

    usage = "\n\nAccepts 2 .csv files represneting read counts for each position.\nUsage: %prog arg"
    parser = OptionParser(usage)

    parser.add_option("-s", "--single", dest="s1_file", help="CSV representing single stranded cut positions.")
    parser.add_option("-v", "--double", dest="v1_file", help="CSV representing double stranded cut positions.")
    parser.add_option("-p", "--outputPars", dest="pars_output", help="PARS score output file destination.")
    parser.add_option("-t", "--outputStruct", dest="struct_output", help="Structure output file destination.")

    (options, args) = parser.parse_args()
    
    if len(args) != 0  or not options.s1_file or not options.v1_file or not options.pars_output or not options.struct_output:
                parser.error("incorrect number of arguments.\n\t\tUse -h to get more information")


    s1_string = ""
    v1_string = ""
    pars_list = [];
    
    try:
        s1 = open(options.s1_file, "r")
        v1 = open(options.v1_file, "r")
        try:
            s1_string = s1.read()
            v1_string = v1.read() 
        finally:
            s1.close()
            v1.close()
    except IOError:
        pass

    li = s1_string.split(',')[:-1]
    s1_list = [int(x) for x in li]
    li = v1_string.split(',')[:-1]
    v1_list = [int(x) for x in li]

    s1_len = len(s1_list)
    v1_len = len(v1_list)

    add2length = max(len(s1_list),len(v1_list)) - min(len(s1_list),len(v1_list))
    if s1_len > v1_len:
#        print "s1 bigger"
        v1_list += [0] * add2length
    if s1_len < v1_len:
#        print "v1 bigger"
        s1_list += [0] * add2length

#    print str(len(s1_list)) + " items in s1_list."
#    print str(len(v1_list)) + " items in v1_list."

    s1_numReads = 0
    v1_numReads = 0
    for position in range(len(s1_list)):
        s1_numReads = s1_numReads + s1_list[position]
        v1_numReads = v1_numReads + v1_list[position]

 #   print "s1 reads = " + str(s1_numReads)
 #   print "v1 reads = " + str(v1_numReads)

#    print "There are " + str(s1_numReads) + " reads in s1 file."
#    print "There are " + str(v1_numReads) + " reads in v1 file."

    pars_file = open(options.pars_output, "w")
    struct_file = open(options.struct_output, "w")



    #2 = Secondary Structure; 1 = single stranded; 0 = N/A
    for position in range(len(s1_list)):
        if s1_list[position] == 0 and v1_list[position] == 0:
 #           print "both are zero"
            struct_file.write("0," + '\n')
            continue
        
        if s1_list[position] == 0 or v1_list[position] == 0:
#            print "one is zero"
            if s1_list[position] == 0:
                struct_file.write("2," + '\n')
            else:
                struct_file.write("1," + '\n')

            continue
        
 #       print "position number " + str(position)
  #      print "s1 val is " + str(s1_list[position])
   #     print "v1 val is " + str(v1_list[position])
    #    print "s1_numreads is " + str(s1_numReads)
        score_s1 = float(s1_list[position]/(s1_numReads+.0))
        score_v1 = float(v1_list[position]/(v1_numReads+.0))
#        print "s1 score: " + str(score_s1)
 #       print "v1 score: " + str(score_v1)
            
        pars_score = float(score_v1/score_s1)
  #      print "pars score: " + str(pars_score)
        pars_ratio = math.log(pars_score, 2) # put it in log2
   #     print "pars ratio: " + str(pars_ratio)

        pars_file.write(str(pars_ratio) + '\n')

        if pars_ratio > 0:
            struct_file.write("2," + '\n')
        else:
            struct_file.write("1," + '\n')

    pars_file.close()
    struct_file.close()
            
