#!/usr/bin/env python
# removes line wraps in sequences
# usage:
# linearize_fasta.py < input.fasta > output.fasta

import sys

seq = ''
for line in sys.stdin.readlines():
    if line.startswith('>'):
        if len(seq) > 0:
            print seq
        print line.strip()
        seq = ''
    else:
        seq += line.strip()

# print last seq
print seq
