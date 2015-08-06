# usage:
# 
# trim_fasta_by_primers.py -i in.fna -f fwd_primer -r rev_primer -o newfasta.fna -n 2
# -n is number of mismatches allowed
import sys, os
import tempfile
from subprocess import Popen, PIPE, STDOUT
from optparse import OptionParser

# note these are not reverse-comlemented
PRIMERS = {}
PRIMERS['V3f-Illumina'] = "CCTACGGGNGGCWGCAG" # http://support.illumina.com/documents/documentation/chemistry_documentation/16s/16s-metagenomic-library-prep-guide-15044223-b.pdf
PRIMERS['V4f-EMP'] = "GTGYCAGCMGCCGCGGTAA" # EMP website
PRIMERS['V4r-EMP'] = "GGACTACNVGGGTWTCTAAT" # EMP website
PRIMERS['V4r-Illumina'] = "GACTACHVGGGTATCTAATCC" # http://support.illumina.com/documents/documentation/chemistry_documentation/16s/16s-metagenomic-library-prep-guide-15044223-b.pdf

def make_option_parser():
    parser = OptionParser(usage="usage: %prog [options] filename",
                          version="%prog 1.0")
    parser.add_option("-i","--input_fasta",
                      default=None,
                      type='string',
                      help="Path to input fasta_file [required]") 
    parser.add_option("-f","--forward_primer",
                      default='GTGYCAGCMGCCGCGGTAA',
                      type='string',
                      help="Forward primer. Default is EMP v4. [default %default]")
    parser.add_option("-r","--reverse_primer",
                      default='ATTAGAWACCCBNGTAGTCC',
                      type='string',
                      help="Reverse primer. Default is EMP v4. Unless --reverse_complement is passed, assumes this is already reverse-complemented. [default %default]")
    parser.add_option("-z","--reverse_complement",
                      action="store_true",
                      default=False,
                      help="Reverse complement the reverse primer before matching (default %default)",)
    parser.add_option("-u","--embalmer_command",
                      default='embalm',
                      type='string',
                      help="Path to embalmer command [required]") 
    parser.add_option("-n","--n_mismatches",
                      default=1,
                      type='int',
                      help="Number of primer mismatches allowed. [default %default]")
    parser.add_option("-p","--padding",
                      default=0,
                      type='int',
                      help="Additional nucleotides to include on either end of the region. If this value is greater than then length of either of the primers the program will simply return up to the end of that primer. [default %default]")
    parser.add_option("-T","--n_threads",
                      default=1,
                      type='int',
                      help="Number of concurrent threads for embalmer process [default %default]")
    parser.add_option("-v","--verbose",
                      action="store_true",
                      default=False,
                      help="Verbose output (default %default)",)
    parser.add_option("-d","--dont_delete",
                      action="store_true",
                      default=False,
                      help="Don't delete temporary files (default %default)",)
    parser.add_option("-o","--output_fp",
                      type="string",
                      default=None,
                      help="Path to output file [default os.path.splitext(os.path.basename(options.query))[0] + '-f<fwprimer>-r<reverseprimer>.fna']",)
    return parser


def run_command(cmd, verbose=False):
    if verbose:
        print cmd
    proc = Popen(cmd,shell=True,universal_newlines=True,stdout=PIPE,stderr=PIPE)
    stdout, stderr = proc.communicate()
    return_value = proc.returncode
    return return_value, stdout, stderr


if __name__ == '__main__':
    parser = make_option_parser()
    (options, args) = parser.parse_args()

    if options.output_fp is None:
        options.output_fp = os.path.splitext(os.path.basename(options.query))[0]
        options.output_fp += '-%s-%s.fna' %(options.forward_primer, options.reverse_primer)
    if options.verbose:
        print "Output file is", options.output_fp


    # write the primers to two different temporary FNA files
    fwd_fp = tempfile.mkstemp()[1]
    rev_fp = tempfile.mkstemp()[1]

    fout = open(fwd_fp,'w')
    fout.write('>Forward\n'+options.forward_primer + '\n')
    fout.close()

    fout = open(rev_fp,'w')
    fout.write('>Forward\n'+options.reverse_primer + '\n')
    fout.close()

    # RC the reverse primer if request (requires QIIME)
    if options.reverse_complement:
        rev_rc_fp =  tempfile.mkstemp()[1]
        cmd = "adjust_seq_orientation.py -i " + rev_fp + " -o " + rev_rc_fp
        ret, so, se = run_command(cmd, options.verbose)
        if ret != 0:
            print so
            print se
            raise ValueError("QIIME adjust_seq_orientation.py failed. Is QIIME installed?")
        cmd = "mv " + rev_rc_fp + " " + rev_fp
        ret, so, se = run_command(cmd, options.verbose)

    # Match each primer against each input sequence
    fwd_hits_fp = tempfile.mkstemp()[1]
    rev_hits_fp = tempfile.mkstemp()[1]
    cmd_parts = [options.embalmer_command, '--forage', options.input_fasta, fwd_fp, fwd_hits_fp]
    cmd_parts += [str(1 - float(options.n_mismatches)/len(options.forward_primer)), str(options.n_threads)]
    if options.verbose:
        print "Aligning forward primer to input sequences..."
    ret, so, se = run_command(' '.join(cmd_parts), options.verbose)
    if ret != 0:
            print so
            print se
            raise ValueError("embalmer alignment of forward primer failed. Is embalmer installed?")

    cmd_parts = [options.embalmer_command, '--forage', options.input_fasta, rev_fp, rev_hits_fp]
    cmd_parts += [str(1 - float(options.n_mismatches)/len(options.reverse_primer)), str(options.n_threads)]
    if options.verbose:
        print "Aligning reverse primer to input sequences..."
    ret, so, se = run_command(' '.join(cmd_parts), options.verbose)
    if ret != 0:
            print so
            print se
            raise ValueError("embalmer alignment of reverse primer failed. Is embalmer installed?")


    # read start/end of each sequence from primer hits
    starts = {} # {seqid:start}
    ends = {} # {seqid:end}
    pos = {}  # {seqid:[start,end]}
    for line in open(fwd_hits_fp,'U'):
        line = line.strip()
        words = line.split('\t')
        starts[words[1]] = int(words[9])

    for line in open(rev_hits_fp,'U'):
        line = line.strip()
        words = line.split('\t')
        ends[words[1]] = int(words[8])

    # get intersection
    for key in starts:
        if ends.has_key(key):
            pos[key] = [starts[key],ends[key]]
    if len(pos) == 0:
        raise ValueError("Error: There are were no valid alignments. Does reverse primer need to be reverse-complemented?")

    # trim fasta
    if options.verbose:
        print "Trimming input fasta file..."
    outf = open(options.output_fp,'w')
    header = ''
    seqid = ''

    # if user requested padding, ensure it doesn't go past the primers
    start_padding = min(len(options.forward_primer), options.padding)
    end_padding = min(len(options.reverse_primer), options.padding)

    for line in open(options.input_fasta,'U'):
        line = line.strip()
        if line.startswith('>'):
            header = line
            seqid = line.split()[0][1:]
        else:
            if pos.has_key(seqid):
                outf.write(header + '\n')
                startix = pos[seqid][0]-1
                endix = pos[seqid][1]-1
                startix = max(startix - start_padding, 0)
                endix = min(endix + end_padding, len(line))
                outf.write(line[startix:endix] + '\n')
    outf.close()

    to_remove = [fwd_fp, fwd_hits_fp, rev_fp, rev_hits_fp]
    if not options.dont_delete:
        for f in to_remove:
            os.remove(f)
