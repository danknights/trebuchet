# launches separate usearch threads for all UDB
# Usage:
# parallel_usearch.py -q query.fasta -r ref_db_directory -o outfile
import sys, os
from optparse import OptionParser
import multiprocessing
from subprocess import Popen, PIPE, STDOUT
from collections import OrderedDict
import math
import time
from datetime import timedelta

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
    parser.add_option("-u","--usearch_command",
                      default='usearch8.0',
                      type='string',
                      help="Path to usearch command [required]") 
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
    parser.add_option("-z","--reverse_complement",
                      action="store_true",
                      default=False,
                      help="Search both strands (default %default)",)
    parser.add_option("-v","--verbose",
                      action="store_true",
                      default=False,
                      help="Verbose output (default %default)",)
    parser.add_option("-l","--split_lines",
                      default=20000,
                      type='int',
                      help="Number of sequences per query split file (default %default)")
    parser.add_option("-o","--output_fp",
                      type="string",
                      default=None,
                      help="Path to output file [default os.path.splitext(os.path.basename(options.query))[0] + '-usearch-out.txt']",)
    return parser

def run_usearch(query_fp, ref_fp, output_fp, usearch_cmd='usearch8.0',nthreads=1,
                max_accepts=2, max_rejects=32, query_cov=1.0, target_cov=0,
                reverse_complement=True, pct_ID=0.97,verbose=False):
    """thread worker function"""
    cmd_dict = OrderedDict()
    cmd_dict[usearch_cmd] = ''
    cmd_dict['-usearch_global'] = query_fp
    cmd_dict['-db'] = ref_fp
    cmd_dict['--blast6out'] = output_fp
    cmd_dict['-id'] = pct_ID
    cmd_dict['-query_cov'] = query_cov
    cmd_dict['-target_cov'] = target_cov
    cmd_dict['-maxaccepts'] = max_accepts
    cmd_dict['-maxrejects'] = max_rejects
    if reverse_complement:
        cmd_dict['-strand'] = 'both'
    else:
        cmd_dict['-strand'] = 'plus'

    cmd_dict['-threads'] = nthreads

    cmd = ' '.join([key + ' ' + str(cmd_dict[key]) for key in cmd_dict])
    if verbose:
        print cmd
    proc = Popen(cmd,shell=True,universal_newlines=True,stdout=PIPE,stderr=PIPE)
    stdout, stderr = proc.communicate()
    return_value = proc.returncode
    return return_value, stdout, stderr

def run_command(cmd, verbose=False):
    if verbose:
        print cmd
    proc = Popen(cmd,shell=True,universal_newlines=True,stdout=PIPE,stderr=PIPE)
    stdout, stderr = proc.communicate()
    return_value = proc.returncode
    return return_value, stdout, stderr

def check_command_args(options, args):
    required = ['query','ref','output_fp']
    if options.query is None:
        raise ValueError('Requires these parameters: ' + ', '.join(required))
    if options.ref is None:
        raise ValueError('Requires these parameters: ' + ', '.join(required))
    
if __name__ == '__main__':
    parser = make_option_parser()
    (options, args) = parser.parse_args()

    if options.output_fp is None:
        options.output_fp = os.path.splitext(os.path.basename(options.query))[0] + '-usearch-out.txt'
    check_command_args(options, args)

    processes = []

    if os.path.isdir(options.ref):
        ref_fps = []
        for fp in os.listdir(options.ref):
            if fp.endswith('.udb'):
                ref_fps.append(os.path.join(options.ref,fp))
        ref_fps = sorted(ref_fps)
    else:
        ref_fps = options.ref.strip().split(',')

    retvals = [''] * len(ref_fps)
    stdouts = [''] * len(ref_fps)
    stderrs = [''] * len(ref_fps)

    output_fps = [options.output_fp + '%03d' %(i) for i in xrange(len(ref_fps))]
    
    # write query file out in batches; do search on each one.
    nseqs = 0
    for line in open(options.query,'U'):
        if line[0] == '>':
            nseqs += 1
    
    print nseqs,'query sequences and',len(ref_fps),'reference databases.'

    starttime = time.time()
    seqcount = 0
    output_file = open(options.output_fp,'w')
    tmp_query_fp = options.output_fp + '_partial_query.tmp'
    tmp_query_file = open(tmp_query_fp,'w')
    tmp_output_fp = options.output_fp + '_partial_output.tmp'
    header = ''
    totalsearches = math.ceil(nseqs/min(nseqs,options.split_lines)) * len(ref_fps)
    nsearches = 0
    for line in open(options.query,'U'):
        tmp_query_file.write(line)
        if not line[0] == '>':
            seqcount += 1
            if seqcount % options.split_lines == 0 or seqcount == nseqs:
                tmp_query_file.close()

                for i, ref_fp in enumerate(ref_fps):
                    currtime = time.time()
                    elapsedtime = currtime - starttime
                    elapsedtimestr = str(timedelta(seconds=round(elapsedtime)))
                    if nsearches > 0:
                        remtime = elapsedtime / float(nsearches) * (totalsearches - nsearches)
                        remtimestr = str(timedelta(seconds=round(remtime)))
                    else:
                        remtimestr = 'Unknown time'
                    print elapsedtimestr, "searching query sequence",seqcount,'of',nseqs,'against ref', i+1,'of',str(len(ref_fps)) + ';',remtimestr,'remaining.'
                    retvals[i], stdouts[i], stderrs[i] = \
                        run_usearch(tmp_query_fp, 
                                ref_fp, tmp_output_fp, 
                                usearch_cmd=options.usearch_command,
                                nthreads=options.nthreads,
                                max_accepts=options.max_accepts,
                                max_rejects=options.max_rejects,
                                query_cov=options.query_coverage,
                                target_cov=options.target_coverage,
                                reverse_complement=options.reverse_complement,
                                pct_ID=options.pct_ID, verbose=options.verbose)
                    nsearches += 1
                    if retvals[i] != 0:
                        sys.stderr.write('Warning: usearch run %d on ref DB %s failed with the following error output:\n' %(i, ref_fp))
                        sys.stderr.write(stderrs[i] + '\n')
                        raise ValueError('USEARCH error - if USEARCH ran out of memory, decrease --split_lines.')
                    else:
                        for line in open(tmp_output_fp,'U'):
                            output_file.write(line)
                if seqcount < nseqs:
                    # if there are more seqs to search, reopen the temp query file
                    tmp_query_file = open(tmp_query_fp,'w')

    # remove temp output files and query files
    os.remove(tmp_query_fp)
    os.remove(tmp_output_fp)
    output_file.close()
