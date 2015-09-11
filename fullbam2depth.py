#!/usr/bin/python
import sys, os
from optparse import OptionParser
import pysam


options = None

DEBUG = False


def WriteLineBreak(output_file, linebreak):
	global options
	linebreak += 1
	if options.linebreak>0 and linebreak == options.linebreak:
		output_file.write('\n')
		linebreak = 0
	return linebreak

if __name__ == '__main__':

	usage = "\n\nIt works with sorted bam files.\nyou can use samtools to sort a bamfile.\nusage: %prog arg"
	parser = OptionParser(usage)

	parser.add_option("-b", "--bamfile", dest="bamfile", help="SORTED BAM/SAM file to calculate depth graph from")
	parser.add_option("-l", "--linebreak", dest="linebreak", help="the number of items to have in every line of output file (default -1 [no line breaks])", type='int', default='-1')

	parser.add_option("-o", "--output", dest="output", help="Output file")
	
	(options, args) = parser.parse_args()

	if len(args) != 0  or not options.bamfile or not options.output:
		parser.error("incorrect number of arguments.\n\t\tUse -h to get more information")



	samfile = pysam.Samfile(options.bamfile)
	output_file = open(options.output, 'w')

	total_read_count = 0
	total_depth_count = 0
	write_order = -1
	pos = 0
	cnt = 0
	linebreak = 0
	reverse = {}
	for read in samfile.fetch():

		#if total_read_count>1000:
		#	exit(0)

		total_read_count += 1
		if read.is_reverse:
			if not read.pos+read.rlen in reverse:
				reverse[read.pos+read.rlen] = 0
			reverse[read.pos+read.rlen] += 1

		else:

			if read.pos == pos:
				cnt += 1

			else:
				# add up reverese strand read counts
				if pos in reverse:
					cnt += reverse[pos]
					del reverse[pos]

				#write the last pos count
				output_file.write(str(cnt)+',')		
				linebreak = WriteLineBreak(output_file, linebreak)

				#keep counts for checksum 
				total_depth_count += cnt		
				write_order += 1

				if DEBUG:
					print pos,',',cnt, pos in reverse
					print '------ next pos', read.pos

					if write_order != pos:
						print 'ERROR order', write_order, pos



				
				# put NA for the columns that don't exist
				for i in xrange(pos+1, read.pos):
					if i in reverse:
						output_file.write(str(reverse[i]) + ',')
						total_depth_count += reverse[i]
						del reverse[i]
					else:
						output_file.write('0,')

					write_order+=1
					linebreak = WriteLineBreak(output_file, linebreak)

				pos = read.pos
				cnt = 1

	output_file.close()


	print 'Total number of reads:', total_read_count
	print 'Sum of the depth file:', total_depth_count

	s = 0
	for p in reverse:
		s += reverse[p]
	if s>0 or total_read_count - total_depth_count>2:
		sys.stderr.write('ERROR in converting to depth\n')
		sys.stderr.write('Total number of reads (%d) is different from the sum of numbers generated (%d) by %d\n' % (total_read_count, total_depth_count, total_read_count-total_depth_count) )
		sys.stderr.write('%d positions are left in reverese list adds up to %d number of reads\n' % (len(reverse), s))

	if DEBUG:
		s = 0
		for p in reverse:
			s += reverse[p]
		print len(reverse), 'left in reverese adds up to', s
		print 'last pos was',pos
		print 'min pos in reverse is', min(reverse.keys())
		print 'max pos in reverse is', max(reverse.keys())








	# last_value = 0
	# last_start = 0
	# for pileupcolumn in samfile.pileup(chromosome_name , start, end):
	# 	#print
	# 	if pileupcolumn.pos>= start and pileupcolumn.pos <= end:
	# 		if pileupcolumn.n != last_value: 
	# 			print last_start, pileupcolumn.pos, last_value
	# 			last_value = pileupcolumn.n
	# 			last_start = pileupcolumn.pos

	samfile.close()
	output_file.close()


