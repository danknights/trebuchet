import sys

count = 0
prevline = '' 
prevprevline = ''
for line in open(sys.argv[1],'r'):
    line = line.strip()
    if count % 2 == 0 and not line.startswith('>'):
        print("Line ",count, "\n")
        print(prevprevline)
        print(prevline)
        print(line)
        sys.exit(0)
    else:
        count += 1
        prevprevline = prevline
        prevline = line
