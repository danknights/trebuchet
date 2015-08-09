# usage
# find_distinct_taxa.py otu_by_otu_matches.txt taxonomy outfile.txt

import sys, os
from optparse import OptionParser

def make_option_parser():
    parser = OptionParser(usage="usage: %prog [options] filename",
                          version="%prog 1.0")
    parser.add_option("-i","--input_fp",
                      default=None,
                      type='string',
                      help="Path to input file: a table of OTU vs OTU pairs in the first two tab-delimited columns. [required]") 
    parser.add_option("-t","--taxonomy_file",
                      default=None,
                      type='string',
                      help="Path to taxonomy file [required]") 
    parser.add_option("-T","--taxonomy_format",
                      default='greengenes',
                      type='string',
                      help="Type of taxonomy format: greengenes or silva [default %default]") 
    parser.add_option("-v","--verbose",
                      action="store_true",
                      default=False,
                      help="Verbose output (default %default)",)
    parser.add_option("-c","--consensus",
                      action="store_true",
                      default=False,
                      help="Report consensus (lowest) taxonomy (default %default)",)
    parser.add_option("-o","--output_fp",
                      type="string",
                      default=None,
                      help="Path to output file [default os.path.splitext(os.path.basename(options.query))[0] + '-taxa.txt']",)
    return parser


if __name__ == '__main__':
    parser = make_option_parser()
    (options, args) = parser.parse_args()


    # build taxonomy map
    print "Loading taxonomy map..."
    taxa = {}
    for line in open(options.taxonomy_file,'r'):
        words = line.strip().split('\t')
        taxon_ID = words[0]
        taxonomy = words[1].replace(';',' ').split()
        taxonomy = [level.strip() for level in taxonomy]
        taxa[taxon_ID] = taxonomy
    
    # Read through otu_by_otu matches
    # For each otu, if it's in the taxonomy dict
    # add it to the list of those with species-level taxonomy
    #
    # if it matches any other otus in a different species,
    # add it to the ambiguous list
    best_labels = {} # {taxon_ID:consensus taxonomy, ...}
    
    potential_list = set()
    bad_list = set()

    for line in open(options.input_fp,'U'):
        words = line.split('\t')
        query = words[0].split()[0]
        ref = words[1].split()[0]
        consensus = os.path.commonprefix([])


        # set this query's taxonomy to the most specific consensus
        if not best_labels.has_key(query):
            best_labels[query] = taxa[query]
        else:
            best_labels[query] = os.path.commonprefix([taxa[ref], best_labels[query]])

        # only keep this as potential if it resolves to species
        if len(taxa[query]) > 6:

            if taxa[query][6] != 's__' and taxa[query][6] != 'unidentified':
                potential_list.add(query)
        
                # add to bad list if its match has a different species
                # ignore the match if it is not resolved to species level
                if len(taxa[ref]) > 6:
                    if taxa[ref][6] != taxa[query][6]:
                        bad_list.add(query)

    # print only those in the good list and not the bad list

    # good_taxa = set()
    # for key in potential_list:
    #     if not key in bad_list:
    #         good_taxa.add('; '.join(taxa[key]))

    good_taxa = set(['; '.join(taxonomy) for taxon_id, taxonomy in best_labels.iteritems()]) 

    outf = open(options.output_fp,'w')
    for taxon in sorted(good_taxa):
        outf.write(taxon + '\n')
    outf.close()

