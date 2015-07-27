# launches separate usearch threads for all UDB
# Usage:
# parallel_usearch.py -q query.fasta -r ref_db_directory -o outfile
import sys, os
from optparse import OptionParser
import multiprocessing
from subprocess import Popen, PIPE, STDOUT
from collections import OrderedDict

def make_option_parser():
    parser = OptionParser(usage="usage: %prog [options] filename",
                          version="%prog 1.0")
    parser.add_option("-q","--query",
                      default=None,
                      type='string',
                      help="Path to query fasta [required]") 
    parser.add_option("-r","--ref",
                      default=None,
                      type='string',
                      help="Path to directory containing one or more reference udb files, or comma-separated list of udb files [required]")
    parser.add_option("-n","--nthreads",
                      default=1,
                      type='int',
                      help="Number of concurrent threads for each usearch process [default %default]")
    parser.add_option("-A","--max_accepts",
                      default=2,
                      type='int',
                      help="Number of hits to attempt to find for each query seq [default %default]")
    parser.add_option("-R","--max_rejects",
                      default=32,
                      type='int',
                      help="Number of rejects before calling a query a failure [default %default]")
    parser.add_option("-I","--pct_ID",
                      default=.97,
                      type='float',
                      help="Percent identity for alignment [default %default]")
    parser.add_option("-Q","--query_coverage",
                      default=1.0,
                      type='float',
                      help="Fraction of query seq to align [default %default]")
    parser.add_option("-T","--target_coverage",
                      default=0.0,
                      type='float',
                      help="Fraction of reference seq to align [default %default]")
    parser.add_option("-v","--verbose",
                      action="store_true",
                      default=False,
                      help="Verbose output (default %default)",)
    parser.add_option("-o","--output_fp",
                      type="string",
                      default=None,
                      help="Path to output file [required]",)
    return parser

def run_usearch(query_fp, ref_fp, output_fp, nthreads=1,
                max_accepts=2, max_rejects=32, query_cov=1.0, target_cov=0,
                pct_ID=0.97):
    """thread worker function"""
    cmd_dict = OrderedDict()
    cmd_dict['usearch7.0'] = ''
    cmd_dict['-usearch_global'] = query_fp
    cmd_dict['-db'] = ref_fp
    cmd_dict['--blast6out'] = output_fp
    cmd_dict['-id'] = pct_ID
    cmd_dict['-query_cov'] = query_cov
    cmd_dict['-target_cov'] = target_cov
    cmd_dict['-maxaccepts'] = max_accepts
    cmd_dict['-maxrejects'] = max_rejects
    cmd_dict['-strand'] = 'both'
    cmd_dict['-threads'] = nthreads

    cmd = ' '.join([key + ' ' + str(cmd_dict[key]) for key in cmd_dict])
    print cmd
    proc = Popen(cmd,shell=True,universal_newlines=True,stdout=PIPE,stderr=PIPE)
    stdout, stderr = proc.communicate()
    return_value = proc.returncode
    return return_value, stdout, stderr


if __name__ == '__main__':
    parser = make_option_parser()
    (options, args) = parser.parse_args()

    processes = []

    if os.path.isdir(options.ref):
        ref_fps = []
        for fp in os.listdir(options.ref):
            if fp.endswith('.udb'):
                ref_fps.append(os.path.join(options.ref,fp))
    else:
        ref_fps = options.ref.strip().split(',')

    retvals = [''] * len(ref_fps)
    stdouts = [''] * len(ref_fps)
    stderrs = [''] * len(ref_fps)
    output_fps = [options.output_fp + '%03d' %(i) for i in xrange(len(ref_fps))]
    for i, ref_fp in enumerate(ref_fps):
        print "starting thread",i + 1,'of', len(ref_fps),'for ref db',ref_fp
        retvals[i], stdouts[i], stderrs[i] = \
            run_usearch(options.query, 
                        ref_fp, output_fps[i], 
                        nthreads=options.nthreads,
                        max_accepts=options.max_accepts,
                        max_rejects=options.max_rejects,
                        query_cov=options.query_coverage,
                        target_cov=options.target_coverage,
                        pct_ID=options.pct_ID)
        if retvals[i] != 0:
            sys.stderr.write('Warning: usearch run %d on ref DB %s failed with the following error output:\n' %(i, ref_fp))
            sys.stderr.write(stderrs[i] + '\n')

    # merge output files
    print 'Merging output files to',options.output_fp
    output_fp = open(options.output_fp,'w')
    for output_fp_i in output_fps:
        for line in open(output_fp_i,'U'):
            output_fp.write(line)
        os.remove(output_fp_i)
    output_fp.close()
    
