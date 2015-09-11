#!/usr/bin/python
import sys, os
import csv
import math
import StringIO
import time
from optparse import OptionParser


class Position(object):
    gene_id = ""
    gene_position = -1
    ctrl_counts = -1
    s1_counts = -1
    v1_counts = -1
    pars_val = 0
    ds_prob = -1
    s1_norm = -1
    v1_norm = -1
    ctrl_norm = -1

    def __init__(self,gene_name,gene_pos,ctrl,s1,v1,pars,ds):
        self.gene_id = gene_name
        self.gene_position = gene_pos
        self.ctrl_counts = ctrl
        self.s1_counts = s1
        self.v1_counts = v1
        self.pars_val = pars
        self.ds_prob = ds

    def print_pos(self):
        print "{Gene: " + str(self.gene_id) + "\tPosition: " + str(self.gene_position) + "\tControl: " + str(self.ctrl_counts) + "\tS1: " + str(self.s1_counts) + "\tV1: " + str(self.v1_counts) + "\tPARS: " + str(self.pars_val) + "\tDS_Prob: " + str(self.ds_prob) + "}"


class Gene(object):
    gene_id = ""
    positions = []

    def __init__(self, gene_name):
        self.gene_id = gene_name

    def print_gene(self):
        print "Gene: " + self.gene_id
        for v in range(len(self.positions)):
            print self.positions[v].print_pos()

    def normalize(self):
        print "Normalizing " + self.gene_id
        total_ctrl = 0
        total_s1 = 0
        total_v1 = 0
        for v in range(len(self.positions)):
            #print "Looking at position " + str(v) + ": " + str(self.positions[v])
            total_ctrl = total_ctrl + self.positions[v].ctrl_counts
            print total_ctrl
            total_s1 = total_s1 + self.positions[v].s1_counts
            total_v1 = total_v1 + self.positions[v].v1_counts

        print "Total control counts for gene + " + self.gene_id + ": " + str(total_ctrl)
        print "Total s1 counts for gene + " + self.gene_id + ": " + str(total_s1)
        print "Total v1 counts for gene + " + self.gene_id + ": " + str(total_v1)
            
        for w in range(len(self.positions)):
            if self.positions[w].ctrl_counts <> 0:
                self.positions[w].ctrl_norm = self.positions[w].ctrl_counts/total_ctrl
                #print "Normalized control count: " + str(self.positions[w].ctrl_norm)
            if self.positions[w].s1_counts <> 0:
                self.positions[w].s1_norm = self.positions[w].s1_counts/total_s1
                #print "Normalized s1 count: " + str(self.positions[w].s1_norm)
            if self.positions[w].v1_counts <> 0:
                self.positions[w].v1_norm = self.positions[w].v1_counts/total_v1
                #print "Normalized v1 count: " + str(self.positions[w].v1_norm)


class Genome(object):
    name = ""
    allPos = []   #This is a list of all positions.
    genes = []
    outfile_pos = ""
    gc = 0.0
    
    def __init__(self, genome_name,outfile_location):
        self.name = genome_name
        self.outfile_pos = outfile_location

    def sort_ctrl(self):
        self.allPos.sort(key = lambda w: w.ctrl_counts)

    def print_genome(self):
        for h in range(len(self.allPos)):
            self.allPos[h].print_pos()

    def write_outfile(self):
        try:
            outfile = open(self.outfile_pos, "w")
        except IOError:
            pass
        
        outfile.write("gene\tposition\tcontrol\tcontrol_log\tcontrol_norm\ts1\ts1_log\ts1_norm\tv1\tv1_log\tv1_norm\tpars\tvienna\tgc\n")
        print len(self.allPos)
        for h in range(len(self.allPos)):
            outfile.write(str(self.allPos[h].gene_id) + "\t")
            outfile.write(str(self.allPos[h].gene_position) + "\t")
            outfile.write(str(self.allPos[h].ctrl_counts) + "\t")

            if self.allPos[h].ctrl_counts <> 0:
                outfile.write(str(math.log(self.allPos[h].ctrl_counts)) + "\t")
            else:
                outfile.write("0\t")

            outfile.write(str(self.allPos[h].ctrl_norm) + "\t")

            outfile.write(str(self.allPos[h].s1_counts) + "\t")

            if self.allPos[h].s1_counts <> 0:
                outfile.write(str(math.log(self.allPos[h].s1_counts)) + "\t")
            else:
                outfile.write("0\t")

            outfile.write(str(self.allPos[h].s1_norm) + "\t")

            outfile.write(str(self.allPos[h].v1_counts) + "\t")

            if self.allPos[h].v1_counts <> 0:
                outfile.write(str(math.log(self.allPos[h].v1_counts)) + "\t")
            else:
                outfile.write("0\t")

            outfile.write(str(self.allPos[h].v1_norm) + "\t")

            outfile.write(str(self.allPos[h].pars_val) + "\t")
            outfile.write(str(self.allPos[h].ds_prob) + "\t")
	    outfile.write(str(self.gc) + "\n")

        outfile.close()



if __name__ == '__main__':

    usage = "\n\nCreates an R compatable data frame given high quality CDS files for S1/V1 and control counts.\nOutput contains GeneId, position, control counts, pars score, vienna lbox scores, and vienna ubox scores. \nUsage: %prog arg"
    parser = OptionParser(usage)

    parser.add_option("-s", "--single", dest="s1_file", help="CDS representing single stranded cut positions.")
    parser.add_option("-v", "--double", dest="v1_file", help="CDS representing double stranded cut positions.")
    parser.add_option("-c", "--control",   dest="control_file", help="CDS representing control read counts.")
    parser.add_option("-g", "--gc",   dest="gc_file", help="File containg GC contents.")
    parser.add_option("-o", "--output", dest="output", help="Output file destination.")
    parser.add_option("-f", "--fold_dir", dest="fold_dir", help="Directory containing Vienna fold files. (.ps)")

    (options, args) = parser.parse_args()
    
    if len(args) != 0  or not options.s1_file or not options.v1_file or not options.output or not options.control_file or not options.fold_dir or not options.gc_file:
                parser.error("incorrect number of arguments.\n\t\tUse -h to get more information")


    s1_string = ""
    v1_string = ""
    control_string = ""
    gc_string = ""
    genome = Genome("Staphylococcus Epidermidis", options.output)

    if options.fold_dir[-1] != '/':
        options.fold_dir = options.fold_dir + "/"

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
            s1_string = s1_string.replace('\t',';')
            v1_string = v1_string.replace('\t',';')
            control_string = control_string.replace('\t',';')
	    gc_string = gc_string.replace('\t',';')
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
    li_gc = gc_string.split("\n")
    li_ubox = [[[]] for x in xrange(len(li_s1))]
#    li_lbox = [[[]] for x in xrange(len(li_s1))]

    for i in range(len(li_s1)):
        li_s1[i] = li_s1[i].split(";")
        li_v1[i] = li_v1[i].split(";")
        li_ctrl[i] = li_ctrl[i].split(";")
	li_gc[i] = li_gc[i].split(";")


    control_totalReads = 0
    s1_totalReads = 0
    v1_totalReads = 0

    print "\n\nChecking input files for correspondence..." 
    for i in range(len(li_s1)):
        print li_s1[i]
        print li_v1[i]
        print li_ctrl[i]
        #Check to be sure gene names match.
        if(li_s1[i][0] != li_v1[i][0] or li_s1[i][0] != li_ctrl[i][0] or li_ctrl[i][0] != li_gc[i][0]):
            print "Gene names don't match up!\n\t" + li_s1[i][0] + " |\t| " + li_v1[i][0] + " |\t| " + li_ctrl[i][0] + "|\t|" + li_gc[i][0]
            sys.exit()

        #Check to be sure length of gene matches.
        if(len(li_s1[i]) != len(li_v1[i]) or len(li_s1[i]) != len(li_ctrl[i])):
            print "\tLengths of gene " + li_s1[i][0] + " don't match!\n"
            sys.exit()

    print "\tControl, S1, and V1 files have correct number of gene and gene lengths.\n"

    print "Calculating total reads to use for normalization."
    print str(len(li_s1))
  #  print str(len(li_v1))
    print str(len(li_ctrl))
#    print li_s1
#    print li_v1
#    print li_ctrl
    #Count total number of reads for normalization.
    for i in range(len(li_s1)):
        for j in range(len(li_s1[i])-1):
            s1_totalReads = s1_totalReads + int(li_s1[i][j+1])
            v1_totalReads = v1_totalReads + int(li_v1[i][j+1])
            control_totalReads = control_totalReads + int(li_ctrl[i][j+1])

    print "\tS1 input file has " + str(s1_totalReads) + " total reads."
    print "\tV1 input file has " + str(v1_totalReads) + " total reads."
    print "\tControl input file has " + str(control_totalReads) + " total reads.\n"
                                                      

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


    #Compute ubox values from Vienna .ps files.
    print "Generating temporary ubox files."
    print "\tSkipping step for now, uncomment this later."
    for i in range(len(li_s1)):
        curFile = options.fold_dir + li_s1[i][0] + "_dp.ps"
        outfile = options.fold_dir + li_s1[i][0] + ".ubox"
        grep_cmd = "grep \'ubox\' " + curFile + " > " + outfile
        os.system(grep_cmd)
        time.sleep(5)
        print "\t " + li_s1[i][0] + ".ubox generated using command: " + grep_cmd
 #   sys.exit()
    
    #Remove the first line from each generated file.
    for i in range(len(li_s1)):
        curFile = options.fold_dir + li_s1[i][0] + ".ubox"
        sed_cmd = "sed -i \'1,2d\' " + curFile
        os.system(sed_cmd)
        time.sleep(5)
        print "\t\t" + li_s1[i][0] + " formatted using command: " + sed_cmd



    for i in range(len(li_s1)):
        ubox_string = ""
        li_ubox[i][0] = li_s1[i][0]
        curFile = options.fold_dir + li_s1[i][0] + ".ubox"
        print "Current target file is " + curFile
        try:
            ubox_file = open(curFile, "r")
            print "\tOpening " + curFile + "."
            try:
                ubox_string = ubox_file.read()
            finally:
                ubox_file.close()  
        except IOError:
            print "IO Error opening " + curFile
            pass

       # print li_s1
        print "UBOX STRING" 
        print ubox_string


        last = 0
        sums = 0.0
        prob_matrix = []
        prob_matrix.append(li_s1[i][0])
        counter = 0
        f = open(curFile)
        for l in f:
            counter = counter + 1
            cols = l.split()
            if int(cols[0])<>last:
                if last<>0:
                    prob_matrix.append(sums)
                    
                for j in xrange(last+1, int(cols[0])):
                   # print i, 0.0
                    prob_matrix.append(0.0)

                last = int(cols[0])
                sums =  float(cols[2]) * float(cols[2])
            else:
                sums += float(cols[2]) * float(cols[2])              
        prob_matrix.append(sums)
        
        print prob_matrix

        #Vienna doesn't usually go to the end. So trim the read information so everything's the same lenght.
        while len(prob_matrix) < len(li_s1[i]) or len(prob_matrix) < len(li_v1[i]) or len(prob_matrix) < len(li_ctrl[i]) or len(prob_matrix) < len(pars[i]):
            print "Trimming the ends."
            del li_s1[i][-1]
            del li_v1[i][-1]
            del li_ctrl[i][-1]
            del pars[i][-1]

        #Build the data into positions->genes->genome objects.
        curGene = Gene(li_ctrl[i][0])

        #Delete the name for each gene stored in each of our computed lists.
        del li_ctrl[i][0]
        del li_v1[i][0]
        del pars[i][0]
        del li_s1[i][0]
        del prob_matrix[0]

        #Sanity check to see if all of the lengths of our lists match up.
        print str(len(prob_matrix))
        print str(len(li_ctrl[i]))
        print str(len(li_s1[i]))
        print str(len(li_v1[i]))
        print str(len(pars[i]))
        if len(prob_matrix) <> len(li_ctrl[i]) or len(prob_matrix) <> len(li_s1[i]) or len(prob_matrix) <> len(li_v1[i]) or len(prob_matrix) <> len(pars[i]):
            print str(len(li_ctrl[i])) + "\t " + str(len(prob_matrix)) + "\t" + str(len(li_v1[i])) + "\t" + str(len(pars[i]))
            print "Length mismatch error.  Terminating execution."
            sys.exit()
        
        
        for n in range(len(li_ctrl[i])):
            li_ctrl[i][n] = int(li_ctrl[i][n])
            curPosition = Position(str(curGene.gene_id), int(n+1), int(li_ctrl[i][n]), int(li_s1[i][n]), int(li_v1[i][n]), float(pars[i][n]), float(prob_matrix[n]))
            curGene.positions.append(curPosition)
            genome.allPos.append(curPosition)

        curGene.print_gene()
       # curGene.normalize()
        genome.genes.append(curGene)
	genome.gc = li_gc[0][1]
	print li_gc[0][1]
        #curGene.print_gene()

        #genome.sort_ctrl()
        #genome.print_genome()
        genome.write_outfile()
            
        #print li_s1[i]
        #print pars[i]
        #print prob_matrix
        #print li_ctrl[i]
        #print sorted(li_ctrl[i])
        #print li_ctrl[i]
        #print str(len(li_s1[i]))
        #print str(len(prob_matrix))

        #li_ubox[i] = prob_matrix

    #print li_ubox
