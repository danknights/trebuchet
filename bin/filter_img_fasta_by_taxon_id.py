#!/usr/bin/env python
# usage:
# python filter_img_fasta_by_taxon_id -i input.fasta -f taxon_id_file -o output.fasta

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
    parser.add_option("-f", "--taxon_id_file",
                      default=None,
                      type='string',
                      help="file containing one taxon ID per line (this or taxon_id required)")
    parser.add_option("-t","--taxon_ids",
                      default=None,
                      type='string',
                      help="Comma-delimited list of taxon IDs to extract (this or taxon_id_file required)",)
    parser.add_option("-v","--verbose",
                      action="store_true",
                      default=False,
                      help="Verbose output (default %default)",)
    parser.add_option("-o","--output_fasta",
                      type="string",
                      default=None,
                      help="Output filename (required)",)
    return parser

if __name__ == '__main__':
    parser = make_option_parser()
    (options, args) = parser.parse_args()
    
    # load ids
    if options.taxon_ids is not None:
        taxon_ids = options.taxon_ids.split(',')
    elif options.taxon_id_file is not None:
        taxon_ids = [line.strip() for line in open(options.taxon_id_file,'U')]
    else:
        raise ValueError('Please supply --taxon_id_file or --taxon_ids')

    # read through fasta, printing only requested seqs
    output_file = open(options.output_fasta,'w')

    count = 0

    include = False
    first_seq = True
    for line in open(options.input_fasta,'U'):
        if line.startswith('>'):
            count += 1
            if options.verbose and count % 100000 == 0:
                print count
            taxon_id = line[1:line.find('_')]
            if taxon_id in taxon_ids:
                if not first_seq:
                    output_file.write('\n')
                first_seq = False
                output_file.write(line.strip() + '\n')
                include = True
            else:
                include = False
        else:
            if include:
                output_file.write(line.strip())
    output_file.write('\n')
