#!/usr/bin/python
import sys, os
import StringIO
from optparse import OptionParser

if __name__ == '__main__':

    usage = "\n\nGiven a FASTA file for CDS regions, cuts the header down to simply the gene id. \nUsage: %prog arg"
    parser = OptionParser(usage)

    parser.add_option("-i", "--input", dest="infile", help="Fasta file location.")
    parser.add_option("-o", "--output", dest="outfile", help="Output location.")
    (options, args) = parser.parse_args()
    
    if len(args) != 0  or not options.infile or not options.outfile:
        parser.error("incorrect number of arguments.\n\t\tUse -h to get more information")

    in_string = ""

    try:
        fasta = open(options.infile, "r")
        outfile = open(options.outfile,"w")
        try:
            in_string = fasta.read()
        finally:
            fasta.close()
    except IOError:
        pass

    fastaFile = in_string.split(">")

    for item in range(len(fastaFile)):
        curString = fastaFile[item]
#        startPos = curString.find("gene=")
#        startPos = startPos + 5
        startPos = curString.find("locus_tag=")
	startPos = startPos + 10
        endPos = curString.find("]")
        geneName = curString[:endPos]
        geneName = geneName[startPos:]

        sequenceStart = curString.find("\n")
        sequence = curString[sequenceStart:]
        sequence = sequence.replace("\n","")
        newString = ">" + geneName + "\n" + sequence + "\n"
        fastaFile[item] = newString

        if(len(sequence) > 0):
            outfile.write(fastaFile[item])


    print fastaFile

    outfile.close()

        
