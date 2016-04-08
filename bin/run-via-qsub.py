# Submit job via qsub
# makes a directory in folder out_basename containing stdout, stderr, pbs file
#
# usage:
# run_via_qsub.py -N jobname --command "command" --wallclock "96:00:00" --out_basename stdout_stderr_basename --queue small --nodes 1 --ppn 8 --mem 32gb
# 
# or for default settings:
# run_via_qsub.py -command "command goes here; another command" 


from __future__ import division
from __future__ import print_function

import sys
import os
from optparse import OptionParser
import random
import string
from subprocess import Popen, PIPE

def make_option_parser():
    parser = OptionParser(usage="usage: %prog [options] filename",
                          version="%prog 1.0")
    parser.add_option("-c", "--command",
                      default=None,
                      type='string',
                      help="command to run (required)")
    parser.add_option("-w", "--wallclock",
                      default='96:00:00',
                      type='string',
                      help="Wallclock HH:MM:SS (default 96:00:00)")
    parser.add_option("-N", "--name",
                      default=None,
                      type='string',
                      help="Name (default random string)")
    parser.add_option("-o", "--out_basedir",
                      default='.',
                      type='string',
                      help="Default directory for output files (default is '.')")
    parser.add_option("-q", "--queue",
                      default="small",
                      type='string',
                      help="queue (default small)")
    parser.add_option("-n", "--nodes",
                      default=1,
                      type='int',
                      help="Number of nodes (default 1)")
    parser.add_option("-p", "--ppn",
                      default=1,
                      type='int',
                      help="Processors per node (default 1)")
    parser.add_option("-m", "--mem",
                      default='16gb',
                      type='string',
                      help="RAM per 'node' (note: on the 'small' queue MSI allows node sharing) (default 16gb)")
    parser.add_option("-e", "--email",
                      default=None,
                      type='string',
                      help="Email address; if provided, sends emails at abort, begin, and end (default None)")    
    parser.add_option("-P", "--print_only",
                      action="store_true",
                      default=False,
                      help="Print output and make PBS file only; do not submit job (default run job)",)
    parser.add_option("-v","--verbose",
                      action="store_true",
                      default=False,
                      help="Verbose output (default %default)",)
    return parser


if __name__ == '__main__':
	# make option parser and parse command line flags
    parser = make_option_parser()
    (options, args) = parser.parse_args()

    if options.name is None:
        options.name = ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(10))
        
    pbs_fp = os.path.join(options.out_basedir,options.name + '.pbs')
    stdout_fp = os.path.join(options.out_basedir,options.name + '_stdout.txt')
    stderr_fp = os.path.join(options.out_basedir,options.name + '_stderr.txt')
    
    # make output dir if necessary
    if not os.path.exists(options.out_basedir):
        if options.verbose:
            print('Output directory ' + options.out_basedir + ' not found. Creating it...')
        os.makedirs(options.out_basedir)

    pbs_lines = []
    pbs_lines.append('#!/bin/bash -l')
    pbs_lines.append('#PBS -l walltime=%s,nodes=%d:ppn=%d,mem=%s' %(options.wallclock,options.nodes,options.ppn,options.mem))
    pbs_lines.append('#PBS -m abe')

    if options.email is not None:
        pbs_lines.append('#PBS -M %s' %(options.email))
    pbs_lines.append('#PBS -N %s' %(options.name))
    pbs_lines.append('#PBS -o %s' %(stdout_fp))
    pbs_lines.append('#PBS -e %s' %(stderr_fp))
    pbs_lines.append('#PBS -q %s' %(options.queue))
    pbs_lines.append(options.command)

    pbs_f = open(pbs_fp,'w')
    pbs_f.write('\n'.join(pbs_lines) + '\n')
    pbs_f.close()
    
    cmd = 'qsub ' + pbs_fp
    print('qsub command is: ' + cmd)
    if not options.print_only:
        print('Running qsub...')
    
        proc = Popen(cmd,shell=True,universal_newlines=True,stdout=PIPE,stderr=PIPE)
        stdout, stderr = proc.communicate()
        print(stdout)
        if len(stderr) > 0:
            print('\nStandard err was:\n' + stderr)

