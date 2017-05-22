# usage
# find_distinct_taxa.py otu_by_otu_matches.txt taxonomy outfile.txt
# output (tab-delimited):
# taxonomy fraction unique, number unique, number recovered, number in database
#
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
    parser.add_option("-o","--output_dir",
                      type="string",
                      default=".",
                      help="Path to output directory [default %default]")
    return parser


if __name__ == '__main__':
    parser = make_option_parser()
    (options, args) = parser.parse_args()

    # build taxonomy map
    print "Loading taxonomy map..."
    taxa = {}
    taxa_db_counts = {} # {taxonomy_str:count}
    for line in open(options.taxonomy_file,'r'):
        words = line.strip().split('\t')
        taxon_ID = words[0]
        taxonomy = words[1].split(';')
        taxonomy = [level.strip() for level in taxonomy]
        if taxonomy[0].startswith('Bacteria'):
            if len(taxonomy) == 7:
                species = taxonomy[6]
                species = species.replace(' cf.','')

                species = species.split()
                if len(species) == 2:
                    species = ' '.join(species)
                elif len(species) > 2:
                    if species[1] == 'sp.':
                        species = ' '.join(species[0:3])
                    else:
                        species = ' '.join(species[0:2])
                else:
                    species = species[0]
                taxonomy = taxonomy[:6]
                taxonomy.append(species)
        taxa[taxon_ID] = taxonomy
        taxonomy_str = '; '.join(taxonomy)
        if not taxa_db_counts.has_key(taxonomy_str):
            taxa_db_counts[taxonomy_str] = 0
        taxa_db_counts[taxonomy_str] += 1
        

    # Read through otu_by_otu matches
    # For each otu, if it's in the taxonomy dict
    # add it to the list of those with species-level taxonomy
    #
    # if it matches any other otus in a different species,
    # add it to the ambiguous list
    best_labels = {} # {taxon_ID:consensus taxonomy, ...}
    
    for line in open(options.input_fp,'U'):
        words = line.split('\t')
        query = words[0].split()[0]
        ref = words[1].split()[0]
        consensus = os.path.commonprefix([])

        # set this query's taxonomy to the most specific consensus
        if not best_labels.has_key(query):
            best_labels[query] = taxa[query]
        else:
            if taxa.has_key(ref):
                best_labels[query] = os.path.commonprefix([taxa[ref], best_labels[query]])

    good_taxa = set()
    for taxon_id, taxonomy in best_labels.iteritems():
        taxonomy_str = '; '.join(taxonomy)
        good_taxa.add(taxonomy_str)

    # For each species (defined as 7 levels with 2 words (not counting "sp.") at 7th level
    # count the fraction of its ref bugs that were ID'd
    species_counts = {} # {species:[count, total]}
    for taxon_id, taxonomy in best_labels.iteritems():
        full_taxonomy = taxa[taxon_id]
        full_taxonomy_str = '; '.join(full_taxonomy)

        if not species_counts.has_key(full_taxonomy_str):
            species_counts[full_taxonomy_str] = [0,0]
        species_counts[full_taxonomy_str][1] += 1
        # plus one for this species if it is resolved to species
        if len(taxonomy) == len(full_taxonomy):
            species_counts[full_taxonomy_str][0] += 1

    output_fp_base = os.path.splitext(os.path.basename(options.input_fp))[0]
    output_fp_base = os.path.join(options.output_dir,output_fp_base)
    output_fp_taxa = output_fp_base + '-taxa.txt'
    output_fp_species = output_fp_base + '-species-resolution.txt'
    
    outf = open(output_fp_taxa,'w')
    for taxon in sorted(good_taxa):
        outf.write(taxon + '\n')
    outf.close()

    # now save species
    outf = open(output_fp_species,'w')
    outf.write('\t'.join(['Taxon','Unique hit fraction', 'Unique hit count', 'Total hit count', 'Total in database']) + '\n')
    for taxon in sorted(species_counts):
        total = taxa_db_counts[taxon]
        a = species_counts[taxon][0]
        b = species_counts[taxon][1]
        frac = a/float(b)
        outf.write('\t'.join([str(x) for x in [taxon, frac, a, b, total]]) + '\n')
    outf.close()

