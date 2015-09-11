#!/usr/bin/python
import sys, os
import csv
import math
import StringIO
import time
from optparse import OptionParser


if __name__ == '__main__':

    usage = "\n\nGiven S1, V1, and Control files for high quality genes, this script will split each gene into individual files, generate a dataframe for R, and instruct R to perform statistics on each gene. \nUsage: %prog arg"
    parser = OptionParser(usage)

    parser.add_option("-s", "--single", dest="s1_file", help="CDS representing single stranded HQ genes.")
    parser.add_option("-v", "--double", dest="v1_file", help="CDS representing double stranded HQ genes.")
    parser.add_option("-c", "--control", dest="control_file", help="CDS representing control HQ genes.")
    parser.add_option("-g", "--gc_contents", dest="gc_file", help="File containing GC content for HQ genes.")
    parser.add_option("-o", "--output", dest="output_dir", help="Output directory.")

    (options, args) = parser.parse_args()
    
    if len(args) != 0  or not options.s1_file or not options.v1_file or not options.output_dir or not options.control_file or not options.gc_file:
                parser.error("incorrect number of arguments.\n\t\tUse -h to get more information")


    s1_string = ""
    v1_string = ""
    control_string = ""
    gc_string = ""

    if options.output_dir[-1] != '/':
        options.output_dir = options.output_dir + "/"

    try:
        s1 = open(options.s1_file, "r")
        v1 = open(options.v1_file, "r")
        cnt = open(options.control_file, "r")
	gc = open(options.gc_file, "r")
        try:
            s1_string = s1.read()
            v1_string = v1.read()
            control_string = cnt.read()
	    gc_string = gc.read()
            #s1_string = s1_string.replace('\t',';')
            #v1_string = v1_string.replace('\t',';')
            #control_string = control_string.replace('\t',';')
        finally:
            s1.close()
            v1.close()
            cnt.close()
	    gc.close()
    except IOError,e:
        print e
        sys.exit()

    li_s1 = s1_string.split("\n")#[:-1]
    li_v1 = v1_string.split("\n")#[:-1]
    li_ctrl = control_string.split("\n")#[:-1]
    li_gc = gc_string.split("\n")#[:-1]

#    print li_gc[25]
#    print li_ctrl[24]
#    sys.exit()

    #Make sure we have the same number of genes in each file.
    if len(li_s1) != len(li_v1) or len(li_v1) != len(li_ctrl) or len(li_ctrl) != len(li_gc):
	print li_gc
	print "\n\n"
	#print li_s1
	#print "\n\n"
	print str(len(li_gc))
	print str(len(li_s1))
        print "Gene number mismatch error."
        sys.exit()

    #print li_s1

    gene_names = []

    for i in range(len(li_s1)):
        if li_s1[i] is None or li_s1[i] == "":
            print "NONE"
            break
#        print "'" + str(li_s1[i]) + "'"
#        print "--_-_---___"
        cur_li_s1 = li_s1[i].replace("\t",";")
#        print "'" + str(cur_li_s1) + "'"
#        print "--_-_---___"
#        print "'" + str(cur_li_s1[0]) + "'"
        
        cur_li_v1 = li_v1[i].replace("\t",";")
        cur_li_ctrl = li_ctrl[i].replace("\t",";")
        cur_s1 = cur_li_s1.split(";")
		
	cur_li_gc = li_gc[i].replace("\t",";")
	cur_gc = cur_li_gc.split(";")
		
#        print "--_-_---___"
#        print cur_s1
       # sys.exit()
        cur_v1 = cur_li_v1.split(";")
        cur_ctrl = cur_li_ctrl.split(";")

#        print cur_li_s1
#        print "_---------------------------------_"
#        print cur_s1
        
        s1_out = ""
        v1_out = ""
        ctrl_out = ""
	gc_out = ""

        #Make sure gene order is matching up and save the set of gene names.
        if cur_s1[0] != cur_v1[0] or cur_v1[0] != cur_ctrl[0] or cur_ctrl[0] != cur_gc[0]:
            print "Gene name mismatch error."
            print "\tS1: " + cur_s1[0] + " ||| V1: " + cur_v1[0] + " ||| Ctrl: " + cur_ctrl[0] + " ||| GC: " + cur_gc[0]
	    print cur_gc
	    print cur_s1
            sys.exit()
        else:
            gene_names.append(cur_s1[0])
            if cur_s1[0] == "":
                continue

            #Set output file paths for current genes.
            s1_out = options.output_dir + cur_s1[0] + "_S1.csv"
            v1_out = options.output_dir + cur_s1[0] + "_V1.csv"
            ctrl_out = options.output_dir + cur_s1[0] + "_control.csv"
            dataframe_out = options.output_dir + cur_s1[0] + "_dataframe.out"
	    gc_out = options.output_dir + cur_s1[0] + "_gc.out"

            #Remove gene names from each list.
            del cur_s1[0]
            del cur_v1[0]
            del cur_ctrl[0]
	    del cur_gc[0]

        #Make sure genes are the same length.
        if len(cur_s1) != len(cur_v1) or len(cur_v1) != len(cur_ctrl):
            print "Gene length mismatch error."
            sys.exit()
            
        
        try:
            s1_outfile = open(s1_out, "w")
            v1_outfile = open(v1_out, "w")
            ctrl_outfile = open(ctrl_out, "w")
	    gc_outfile = open(gc_out, "w")
        except IOError,e:
            print e
            sys.exit()


        print "Writing gene info for gene " + cur_s1[0] + "."
        s1_outfile.write(li_s1[i])
        v1_outfile.write(li_v1[i])
        ctrl_outfile.write(li_ctrl[i])
	gc_outfile.write(li_gc[i])

        s1_outfile.close()
        v1_outfile.close()
        ctrl_outfile.close()
	gc_outfile.close()

        print "Generating dataframe for gene " + cur_s1[0]
        commandString = "python /users/aprice67/scripts/hq_dataframer.py -s " + s1_out + " -v " + v1_out + " -c " + ctrl_out + " -g " + gc_out + " -o " + dataframe_out + " -f /users/aprice67/Vienna2/Progs &"
        print "\tCommand string: " + commandString
        os.system(commandString)

