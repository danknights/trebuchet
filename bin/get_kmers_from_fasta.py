# kmer count for a fasta file
# counts number of unique kmers of size k in a fasta file
# usage:
# get_kmers_from_fasta.py -i input.fasta -k 65 -o input_kmers_65.fasta

import sys
import os
from optparse import OptionParser

def make_option_parser():
    parser = OptionParser(usage="usage: %prog [options] filename",
                          version="%prog 1.0")
    parser.add_option("-i", "--input_fasta",
                      default=None,
                      type='string',
                      help="input fasta (required)")
    parser.add_option("-k", "--kmer_size",
                      default=None,
                      type='int',
                      help="Length of kmers (required)")
    parser.add_option("-t","--output_type",
                      type='string',
                      default='fasta',
                      help=("Output type: fasta, text, or count (print count only)"
                      		" (default %default)"))
    parser.add_option("-v","--verbose",
                      action="store_true",
                      default=False,
                      help="Verbose output (default %default)",)
    parser.add_option("-o","--output_file",
                      type="string",
                      default=None,
                      help="Output filename (required unless --output_type is 'count')",)
    return parser


if __name__ == '__main__':
	# make option parser and parse command line flags
    parser = make_option_parser()
    (options, args) = parser.parse_args()
    
    k = options.kmer_size
    input_fp = options.input_fasta
    
    # create and open output file if needed
    if options.output_file is not None:
        output_file = open(options.output_file,'w')


    kmers = set() # set to hold unique kmers, faster than dict or ordered list

    # print new kmers as we read them
    # this will reduce memory requirements because
    # we don't have to store the sequence IDs
    count = 0
    for line in open(input_fp,'U'):

        if line.startswith('>'):
        	# if this is a sequence header

            count += 1
            if options.verbose and count % 100000 == 0:
                print count
                
            # extract only the sequence ID (split on whitespace, first element)
            seq_id = line[1:].split()[0]
        else:
        	# if this line contains sequence

        	# remove trailing whitespace characters
            seq = line.strip()
            
            if len(seq) >= k:
	            # if this sequence is long enough to have a k-mer

				# step through sequence on base at a time
                for i in xrange(0,len(seq)-k + 1):
                	# current kmer is a substring of length k
                    kmer = seq[i:i+k]
                    
                    if not kmer in kmers:
                    	# this kmer is not in the current set of kmers; add it
                        kmers.add(kmer)
                        
                    	# print the kmer (and its sequence ID if fasta output)
                        if options.output_type == 'text' \
                        		or options.output_type == 'fasta':
                            if options.output_type == 'fasta':
                            	kmer_id = '>%s_%09d' %(seq_id,i)
                                output_file.write(kmer_id + '\n')
                            output_file.write(kmer + '\n')

    if options.output_type == 'count':
	    # print count if requested
        if options.output_file is not None:
            output_file.write(str(len(kmers)) + '\n')
        else:
            print len(kmers)
    else:
        output_file.close()

