# usage
# find_distinct_taxa.py otu_by_otu_matches.txt taxonomy outfile.txt

import sys

# build taxonomy map
print "Loading taxonomy map..."
taxa = {}
for line in open(sys.argv[2],'r'):
    words = line.strip().split('\t')
    taxonomy = words[1].split('; ')
    taxa[words[0]] = taxonomy
	
# Read through otu_by_otu matches
# For each otu, if it's in the taxonomy dict
# add it to the list of those with species-level taxonomy
#
# if it matches any other otus in a different species,
# add it to the ambiguous list
potential_list = set()
bad_list = set()

for line in open(sys.argv[1],'U'):
	words = line.split('\t')
	query = words[0]
	ref = words[1]
	
	# only keep this as potential if it resolves to species
	if not taxa[query][6] == 's__':
		potential_list.add(query)
	
		# add to bad list if its match has a different species
		if taxa[ref][6] != taxa[query][6]:
			bad_list.add(query)

# print only those in the good list and not the bad list
good_taxa = set()
for key in potential_list:
	if not key in bad_list:
		good_taxa.add('; '.join(taxa[key]))

outf = open(sys.argv[3],'w')
for taxon in sorted(good_taxa):
	outf.write(taxon + '\n')
outf.close()
	
