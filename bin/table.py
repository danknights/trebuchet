#!/usr/bin/env python
import sys

counts = {}
for line in sys.stdin.readlines():
    line = line.strip()
    if not counts.has_key(line):
       counts[line] = 1
    else:
       counts[line] += 1


counts = sorted(counts.iteritems(),key=lambda xx: xx[1])[::-1]

for pair in counts:
    print pair[0] + ':', pair[1]



