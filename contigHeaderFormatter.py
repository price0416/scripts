#!/usr/bin/python
import sys, os
import csv
import math
import StringIO
import re
from optparse import OptionParser

if __name__ == '__main__':

    usage = "\n\nGiven a text file of fasta header files from Trinity contigs, reformats them so they work with makeblastdb.\nUsage: %prog arg\nUSE: http://users-birc.au.dk/biopv/php/fabox/header_editor.php"
    parser = OptionParser(usage)

    parser.add_option("-i", "--input", dest="in_file", help="Input file containing contig headers in Trinity contig format. e.g: TR1573|c0_g1_i1 len=285 path=[525:0-284] [-1, 525, -2] ")
    parser.add_option("-o", "--output", dest="out_file", help="Output file containing modified headers that work with makeblastdb.")

    (options, args) = parser.parse_args()

    if len(args) != 0  or not options.in_file or not options.out_file:
        parser.error("incorrect number of arguments.\n\t\tUse -h to get more information.")

    infile_string = ""
    out_string = ""
    li_infile = []

    try:
        infile = open(options.in_file, "r")
        print "Opened " + options.in_file + " for reading."

        try:
            infile_string = infile.read()
            print "Read " + options.in_file + ". Length is " + str(len(infile_string))
        finally:
            infile.close()
    except IOError, e:
        print "IO Error. Halting execution."
        print e[0], e[1]

    li_infile = infile_string.split("\n")[:-1]

    reformattedList = []
    for i in range(len(li_infile)):
        curString = ""

        #Handle the TR/Id split and c0_g9_i1 format section.
        li_infile[i] = li_infile[i].split(" ")
        splitPos = re.search("\d", li_infile[i][0].split("|")[0])
        splitPos = splitPos.start()
        trSide = li_infile[i][0].split("|")[0][:splitPos]
        idSide = li_infile[i][0].split("|")[0][splitPos:]
        curString = curString + str(trSide) + "|" + str(idSide) + "|" + str(li_infile[i][0].split("|")[1])

        #Handle the len=... split.
        splitPos = re.search("\d", li_infile[i][1])
        splitPos = splitPos.start()
        lenStr = li_infile[i][1][splitPos:]
        curString = curString + "| length " + str(lenStr)


        reformattedList.append(curString)


    print reformattedList
    try:
        outfile = open(options.out_file, "w")
        print "Opened " + options.out_file + " for writing."
        for i in range(len(reformattedList)):
            outfile.write(reformattedList[i] + "\n")
    finally:
        outfile.close()


