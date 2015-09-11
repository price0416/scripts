#!/usr/bin/python
import sys, os
import csv
import math
import StringIO
from optparse import OptionParser

if __name__ == '__main__':


    usage = "\n\nGiven cut site count files for coding regions for both S1 and V1, this script will identify genes that meet minimum cutoff levels and produce\nnew S1 and V1 count files for the genes that pass the cutoff criteria.\nUsage: %prog arg"
    parser = OptionParser(usage)

    parser.add_option("-s", "--single", dest="s1_file", help="file representing single stranded cut positions.")
    parser.add_option("-v", "--double", dest="v1_file", help="file representing double stranded cut positions.")
    parser.add_option("-c", "--control", dest="control_file", help="file repesenting control read counts for this experiment.")
    parser.add_option("-1", "--s1_output", dest="s1_output", help="S1 output file destination.")
    parser.add_option("-2", "--v1_output", dest="v1_output", help="V1 output file destination.")
    parser.add_option("-3", "--control_output", dest="control_output", help="Control file output file destination.")
    parser.add_option("-4", "--fasta_input", dest="fasta_input", help="Fasta file with sequences for coding regions.")
    parser.add_option("-f", "--fasta_output", dest="fasta_output", help="Fasta file representing coding sequences for selected genes.")
    parser.add_option("-g", "--gc_output", dest="gc_output", help="Output file for calculated GC content of selected genes.")
    parser.add_option("-t", "--strict", dest="strict", help="Usage: strict=true; strict=false; Use this flag if you want to select genes that have good coverage across the whole gene with minimal unknown positions.\nNot using this flag will default to selecting genes that overall have average pars scores above the given threshold.")
    parser.add_option("-d", "--threshold", dest="thresh", help="Threshold for average acceptable PARS score.")

    (options, args) = parser.parse_args()
    
    if len(args) != 0  or not options.s1_file or not options.v1_file or not options.s1_output or not options.v1_output or not options.control_output or not options.control_file or not options.strict or not options.thresh or not options.fasta_input or not options.fasta_output:
                parser.error("incorrect number of arguments.\n\t\tUse -h to get more information")
                

    s1_string = ""
    v1_string = ""
    control_string = ""
    strict_mode = options.strict.lower()
    strict_mode = strict_mode.strip()

    if(strict_mode != "true" and strict_mode != "false"):
        print("Strict mode must be true or false!")
        sys.exit()

    threshold = options.thresh
    threshold = float(threshold)
    

    try:
        control = open(options.control_file, "r")
        print "Opened " + options.control_file
        s1 = open(options.s1_file, "r")
        print "Opened " + options.s1_file
        v1 = open(options.v1_file, "r")
        print "Opened " + options.v1_file
        fasta = open(options.fasta_input, "r")
        print "Opened " + options.fasta_input
        s1_outfile = open(options.s1_output,"w")
        v1_outfile = open(options.v1_output, "w")
        control_outfile = open(options.control_output, "w")
        fasta_outfile = open(options.fasta_output, "w")
	gc_outfile = open(options.gc_output, "w")
        try:
            control_string = control.read()
            print "Read " + options.control_file + ". Length is " + str(len(control_string))
            s1_string = s1.read()
            print "Read " + options.s1_file + ". Length is " + str(len(s1_string))
            v1_string = v1.read()
            print "Read " + options.v1_file + ". Length is " + str(len(v1_string))
            fasta_string = fasta.read()
            print "Read " + options.fasta_input + ". Length is " + str(len(fasta_string))
            s1_string = s1_string.replace('\t',';')
            print "Replace tab with ; in  " + options.s1_file + ". Length is " + str(len(s1_string))
            v1_string = v1_string.replace('\t',';')
            print "Replace tab with ; in  " + options.v1_file + ". Length is " + str(len(v1_string))
            control_string = control_string.replace('\t',';')
            print "Replace tab with ; in  " + options.control_file + ". Length is " + str(len(control_string))
            #fasta_string = fasta_string.replace('_','')  #Dont do this in ML
            #print "Replace _ with blank  " + options.fasta_input + ". Length is " + str(len(fasta_string))
        finally:
            s1.close()
            v1.close()
            control.close()
            fasta.close()
            print "Closed all files."
    except IOError, e:
        print "IO Error. Halting execution."
        print e[0], e[1]

    print "Length of current s1,v1,control before split: " + str(len(s1_string)) + " \t" + str(len(v1_string)) + " \t" + str(len(control_string))

    li_s1 = s1_string.split("\n")[:-1]
    li_v1 = v1_string.split("\n")[:-1]
    li_ctrl = control_string.split("\n")[:-1]
    fasta_list = fasta_string.split(">")[1:]
    #print fasta_list

    print "Length of current s1,v1,control after split: " + str(len(li_s1)) + " \t" + str(len(li_v1)) + " \t" + str(len(li_ctrl))

    for i in range(len(li_s1)):
        li_s1[i] = li_s1[i].split(";")
        li_v1[i] = li_v1[i].split(";")
        li_ctrl[i] = li_ctrl[i].split(";")
        
    for i in range(len(fasta_list)):    
        fasta_list[i] = fasta_list[i].split("\n")

    s1_totalReads = 0
    v1_totalReads = 0

    print "\n\nChecking input files for correspondence..." 
    for i in range(len(li_s1)):
        #Check to be sure gene names match.
        if(li_s1[i][0] != li_v1[i][0] or li_s1[i][0] != li_ctrl[i][0]):
            print "Gene names don't match up!\n\t" + li_s1[i][0] + " |\t| " + li_v1[i][0] + " |\t| " + li_ctrl[i][0]
            sys.exit()

        #Check to be sure length of gene matches.
        if(len(li_s1[i]) != len(li_v1[i]) or len(li_s1[i]) != len(li_ctrl[i])):
            print "s1: " + str(len(li_s1[i])) + "\tv1: " + str(len(li_v1[i])) + "\tctrl: " + str(len(li_ctrl[i]))
            print "\tLengths of gene " + li_s1[i][0] + " don't match!\n" + li_v1[i][0] + "\t" +li_ctrl[i][0]
            del li_s1[i]
            del li_v1[i]
            del li_ctrl[i]
            print "DELELTED THIS GENE-----"
            #sys.exit()

    print "\tControl, S1, and V1 files have correct number of gene and gene lengths.\n"

    print "Calculating total reads to use for normalization."
    #Count total number of reads for normalization.
    for i in range(len(li_s1)):
        for j in range(len(li_s1[i])-1):
            s1_totalReads = s1_totalReads + int(li_s1[i][j+1])
            v1_totalReads = v1_totalReads + int(li_v1[i][j+1])

    print "\tS1 input file has " + str(s1_totalReads) + " total reads."
    print "\tV1 input file has " + str(v1_totalReads) + " total reads.\n"


    print "Calculating PARS scores..."
    #Calculate a set of PARS scores that correspond to S1/V1 positions.
    pars = [[] for x in xrange(len(li_s1))]
    for i in range(len(li_s1)):
        pars[i].append(li_s1[i][0])
        
        for j in range(len(li_s1[i])-1):
            cur_s1 = float(li_s1[i][j+1])
            cur_v1 = float(li_v1[i][j+1])

            #Ignore positions where either file shows a 0 read count value.
            if(cur_s1 == 0.0 or cur_v1 == 0.0):
                pars[i].append("0.0")
                continue

            #Calculate each position and normalize for total number of reads.
            score_s1_norm = float((cur_s1)/(s1_totalReads+.0))
            score_v1_norm = float((cur_v1)/(v1_totalReads+.0))
            pars_score_norm = float(score_v1_norm/score_s1_norm)
            pars_ratio_norm = math.log(pars_score_norm, 2)
            pars[i].append(str(pars_ratio_norm))

    print "\tPARS values calculated."
    print "\tPerforming sanity check on computed PARS data."


    #Check for agreement between pars and initial data.
    for i in range(len(li_s1)):
        #Check to be sure gene names match.
        if(li_s1[i][0] != pars[i][0]):
            print "Gene names don't match up!\n\t" + li_s1[i][0] + " |\t| " + li_pars[i][0]
            sys.exit()

        #Check to be sure length of gene matches.
        if(len(li_s1[i]) != len(pars[i])):
            print "\tLengths of gene " + li_s1[i][0] + " don't match!\n"
            sys.exit()

    print "\t\tPARS data list corresponds to S1, V1, and control input files in number of genes and positions.\n"

    print "Computing average PARS score for each gene."
    if(strict_mode == "true"):
        print "\tStrict mode is enabled."
        
    #Compute a positive average pars score for each gene.
    averages = []
    passFail = []
    overThresh = 0
    for i in range(len(pars)):
        curSum = 0
        curLen = len(pars[i])-1
        for j in range(len(pars[i])-1):
            curScore = float(pars[i][j+1])

            #Don't consider 0 values as positions when computing the average score, unless in strict mode.
            if(curScore == 0.0):
                if(strict_mode == "true"):
                    continue
                else:
                    curLen = curLen - 1 
                    continue
                
            #Convert negative scores to their absolute values for average score calculation.
            if(curScore < 0.0):
                curScore = abs(curScore)

            curSum = curSum + curScore

        #Length becomes zero if there are no reads for a whole gene and we are not in strict mode.    
        if(curLen > 0.0):
            avgScore = curSum/float(curLen)
            if(avgScore > threshold):
                overThresh = overThresh + 1
                averages.append(avgScore)
                passFail.append(1)
            else:
                averages.append(avgScore)
                passFail.append(0)
        else:
            passFail.append(0)
            averages.append(0)

    print "\tChecking computed values for sanity..."
    for i in range(len(pars)):
        #Check to be sure length of gene matches.
        if(len(pars) != len(averages) or len(pars) != len(passFail)):
            print "\tLengths of gene " + pars[i][0] + " don't match!\n"
            sys.exit()

    print "\t\tComputed values appear sane.\n"
    

    print str(overThresh) + " genes passed the threshold..."
    print "\tWriting accepted genes to " + str(options.s1_output)
    print "\tWriting accepted genes to " + str(options.v1_output)
    print "\tWriting accepted genes to " + str(options.control_output)

    tmpCounter = 0 

    for index in range(len(passFail)):
#    for index in range(6):
#	print "index: " + str(index)
        if(passFail[index] == 1):
#	    if(tmpCounter > 500):
#	        sys.exit()
#	    else:
#		tmpCounter = tmpCounter + 1
#		print "-----counter------" + str(tmpCounter)

            #Write sequence info to output fasta file.
	    print "Length of fasta list is: " + str(len(fasta_list))
            for pos in range(len(fasta_list)):
                print fasta_list[pos][0] + "\t "  + li_ctrl[index][0]
                if fasta_list[pos][0] == li_ctrl[index][0]:
		    print "\t\tMATCH"
                    fasta_outfile.write(">" + fasta_list[pos][0] + "\n" + fasta_list[pos][1] + "\n")
		    gc_content = 0.0
		    gc_percent = 0.0
		    cur_gene_len = len(fasta_list[pos][1])
		    for nuc in range(len(fasta_list[pos][1])):
		    	if fasta_list[pos][1][nuc] == "G" or fasta_list[pos][1][nuc] == "C":
			    gc_content = gc_content + 1.0
							
#		    print "GC CONTENT: " + str(gc_content) + " LENGTH: " + str(cur_gene_len) + "\n"
		    gc_percent = gc_content / cur_gene_len
					
		    gc_outfile.write(fasta_list[pos][0] + "\t" + str(gc_percent) + "\n")
					
                    #break
            
            #Write selected gene to S1 outfile.
            s1_outfile.write(li_s1[index][0] + "\t")
#	    print "BEGIN GENE: " + li_s1[index][0]

#	    print "LEN OF LiS1 is " + str(len(li_s1))

            for pos in range(len(li_s1[index])-1):
#	    for pos in range(6):
#		print "\t-----Pos: " + str(pos) + "\t-----Len(li_s1[index]): " + str(len(li_s1[index])) + " if(pos != len(li_s1[index])-2) " 
                if(pos != len(li_s1[index])-2):
                    s1_outfile.write(li_s1[index][pos+1] + ";")
#		    print "li_s1[" + str(index) + "][" + str(pos+1) + "]: " + li_s1[index][pos+1] + "/" + str(len(li_s1[index])) 
#		    print "\t\tpos is != to " + str(len(li_s1[index])-2)
                else:
#		    print "LAST ITEM!---------------------"
#		    print "li_s1[" + str(index) + "][" + str(pos+1) + "]: " + li_s1[index][pos+1] 
                    s1_outfile.write(li_s1[index][-1] + "\n")
                    break

            #Write selected gene to V1 outfile.
            v1_outfile.write(li_v1[index][0] + "\t")
            for pos in range(len(li_v1[index])-1):
                if(pos != len(li_v1[index])-2):
                    v1_outfile.write(li_v1[index][pos+1] + ";")
                else:
                    v1_outfile.write(li_v1[index][-1] + "\n")
                    break

            #Write selected gene to control outfile.
            control_outfile.write(li_ctrl[index][0] + "\t")
            for pos in range(len(li_ctrl[index])-1):
                if(pos != len(li_ctrl[index])-2):
                    control_outfile.write(li_ctrl[index][pos+1] + ";")
                else:
                    control_outfile.write(li_ctrl[index][-1] + "\n")
                    break

            #Write sequence info to output fasta file.
            

    print "\n\nSelected genes written successfully.\n\n"
    
    s1_outfile.close()
    v1_outfile.close()
    control_outfile.close()
    fasta_outfile.close()
    gc_outfile.close()

