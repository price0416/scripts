#!/usr/bin/python
import sys, os
import StringIO
from optparse import OptionParser

if __name__ == '__main__':

    usage = "\n\nGiven reciprocal BLAST output for two organisms, this script will identify the set of non-ambiguous contigs/genes.\nUsage: %prog arg"
    parser = OptionParser(usage)

    parser.add_option("-c", "--contigs", dest="contig_file", help="Blast output (tab format) for contigs blasted against reference genome features database.")
    parser.add_option("-f", "--features", dest="features_file", help="Blast output (tab format) for reference genome features blasted against contigs database.")
    parser.add_option("-g", "--fasta", dest="fasta_file", help="Fasta file of reference genome features (complete) as downloaded from NCBI.")
    parser.add_option("-o", "--output", dest="output", help="Output file destination.")
  
    (options, args) = parser.parse_args()
    if len(args) != 0 or not options.contig_file or not options.features_file or not options.fasta_file or not options.output:
        parser.error("incorrect number of arguments.\n\t\tUse -h to get more information")

    contig_string = ""
    feature_string = ""
    fasta_string = ""
    contig_dict = {}
    feature_dict = {}
    
    #Open files and read to strings.
    try:
        contigs = open(options.contig_file, "r")
        print "Opened " + options.contig_file
        features = open(options.features_file, "r")
        print "Opened " + options.features_file
        fasta = open(options.fasta_file, "r")
        print "Opened " + options.fasta_file
        
        try:
            contig_string = contigs.read()
            print "\tRead " + options.contig_file + "."
            feature_string = features.read()
            print "\tRead " + options.features_file + "."
            fasta_string = fasta.read()
            print "\tRead " + options.fasta_file + "."
        finally:
            contigs.close()
            features.close()
            fasta.close()
            print "Closed all files."            
    except IOError, e:
        print "IO Error. Halting execution."
        print e[0], e[1]
        
    #Split each file into lists.
    li_contigs  = contig_string.split("\n")[:-1]
    li_features = feature_string.split("\n")[:-1]
    li_fasta    = fasta_string.split(">")[:-1]

    #Split contigs rows into individual items.
    for i in range(len(li_contigs)):
        li_contigs[i] = li_contigs[i].split("\t")
        
    for i in range(len(li_features)):
        li_features[i] = li_features[i].split("\t") 

    #Column values for BLAST input:
    #0:  query name
    #1:  subject name
    #2:  percent identities
    #3:  aligned length
    #4:  number of mismatched positions
    #5:  number of gap positions
    #6:  query sequence start
    #7:  query sequence end
    #8:  subject sequence start
    #9:  subject sequence end
    #10: e-value
    #11: bit score
    
    contig_one2one = {}
    feature_one2one = {} 
	
    #A li_contigs row now looks like this: 
    #['TR|845|c4_g7_i2|', 'NC_004459.3_gene_3035', '99.27', '1509', '11', '0', '703', '2211', '1509', '1', '0.0', ' 2726']
    for i in range(len(li_contigs)):
        #Trim up the ID's so they match in both input files. Will look like this: 845|c4_g7_i2
        splitPos = li_contigs[i][0].index("|")
        li_contigs[i][0] = li_contigs[i][0][splitPos+1:]
        if li_contigs[i][0][-1] == "|":
            li_contigs[i][0] = li_contigs[i][0][:len(li_contigs[i][0])-1]
        
        #Identify only one to one mappings first.  
        contigID = li_contigs[i][0]
        if contigID in contig_one2one:
            contig_one2one[contigID] += 1
        else:
            contig_one2one[contigID] = 1

    #A li_features row now looks like this:
    #['lcl|NC_004459.3_gene_1', 'tr|845|c4_g7_i2', '99.10', '1002', '9', '0', '1', '1002', '3253', '2252', '0.0', ' 1801']
    for i in range(len(li_features)):
        splitPos = li_features[i][0].index("|")
        li_features[i][0] = li_features[i][0][splitPos+1:]
        
        #Identify only one to one mappings first.
        featureID = li_features[i][0]
        if featureID in feature_one2one:
            feature_one2one[featureID] += 1
        else:
            feature_one2one[featureID] = 1

    contigCount = 0
    uniqueList = []
    print "Identifying contig BLAST unique hits..." 
    for key, val in contig_one2one.iteritems():
        if val == 1:
            #print "\t" + key
            contigCount = contigCount + 1

    print "\t" + str(contigCount) + " total unique genes."
    print "\n"
	
    featureCount = 0
    print "Identifying feature BLAST unique hits..." 
    for key, val in feature_one2one.iteritems():
        if val == 1:
            #print "\t" + key
            featureCount = featureCount + 1    
    
    print "\t" + str(featureCount) + " total unique genes.\n"
    
	#Create dictionaries for unique hits: {ID -> [matchID, %identity, length]}
    print "Structuring unique hits..."
    for i in range(len(li_contigs)):
        if li_contigs[i][0] in contig_one2one:
            if contig_one2one[li_contigs[i][0]] == 1:
                contig_dict[li_contigs[i][0]] = [li_contigs[i][1], li_contigs[i][2], li_contigs[i][3]]
    print contig_dict
	
    for i in range(len(li_features)):
        if li_features[i][0] in feature_one2one:
            if feature_one2one[li_features[i][0]] == 1:
                feature_dict[li_features[i][0]] = [li_features[i][1], li_features[i][2], li_features[i][3]]
    print feature_dict

    print "\tContig dictionary contains " + str(len(contig_dict)) + " unique items."
    print "\tFeature dictionary contains " + str(len(feature_dict)) + " unique items."
	
    #Identify unique reciprocal hits.
    print "Identifying unique reciprocal hits." 
    for contigKey, contigData in contig_dict.iteritems():
        for featureKey, featureData in feature_dict.iteritems():
            contigMatch = contigData[0]
            contigIdent = contigData[1]
            contigLen   = contigData[2]
            featureIdent = featureData[1]
            featureLen   = featureData[2]
            #print contigMatch + " <-----"
            if contigMatch in feature_one2one and float(contigIdent) > 95 and int(contigLen) > 200 and float(featureIdent) > 95 and int(featureLen) > 200:
                uniqueList.append([contigKey,contigMatch])
            break

    print "\tFound reciprocal unique hits: " + str(len(uniqueList)) + " matches."
    
    #Now, go through the fasta file and parse out the locus_tags and locations for these genes.
    #Each item in li_fasta currently looks like: lcl|NC_004459.3_gene_1 [locus_tag=VV1_RS00005] [location=complement(1..1002)] ...
    print "Finding locus tags and location information for unique reciprocal matches..."
    for i in range(len(li_fasta)):
        if not "|" in li_fasta[i]:
            continue
        else:
            splitPos = li_fasta[i].index("|")
            li_fasta[i] = li_fasta[i][splitPos+1:]
            splitPos = li_fasta[i].index("[")
            fasta_featureID = li_fasta[i][:splitPos-1]
            fasta_fetureID = fasta_featureID.strip()
            for j in range(len(uniqueList)):
                if uniqueList[j][1] == fasta_featureID:
                    locus_tag = li_fasta[i].index("=")
                    stopPos = li_fasta[i].index("]")
                    uniqueList[j].append(li_fasta[i][locus_tag+1:stopPos])
                    locationStart = li_fasta[i].index("n=")  #This is the end of the 'location=' part.
                    locationEnd = li_fasta[i].rindex("]")
                    location = li_fasta[i][locationStart+2:locationEnd]
                    uniqueList[j].append(location)
					
    print "Done." 
    print uniqueList

	
    #Write the output to a file.
    print "Writing output..."
    try: 
        outfile = open(options.output, "w")
        print "Opened " + options.output + " for writing."
        for i in range(len(uniqueList)):
            for item in uniqueList[i][:-1]:
                outfile.write(item + "\t")
            outfile.write(uniqueList[i][-1] + "\n")
    finally:
        outfile.close()
    print "Process complete."
	
	
	
	