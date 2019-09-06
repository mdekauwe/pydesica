#!/usr/bin/env python

"""
Wrapper script to send off all the qsub scripts

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (06.09.2019)"
__email__ = "mdekauwe@gmail.com"

import subprocess
import sys
import os

def make_qsub_file(qsub_fname, pft, node):

    s = """
#!/bin/bash

#PBS -m ae
#PBS -P w35
#PBS -q normal
#PBS -M mdekauwe@gmail.com
#PBS -l mem=16GB
#PBS -l ncpus=16
#PBS -l walltime=03:00:00
#PBS -l wd
#PBS -j oe
#PBS -l other=gdata1

module load dot
source activate sci

python src/run_sensitivity_exp.py %s 16 %d
        """ % (pft, node)

    f = open(qsub_fname, 'w')
    f.write(s)
    f.close()


for pft in ["rf", "wsf", "dsf", "grw", "saw"]:
    count = 0
    node = 1
    total_exp = (5**6) * (3**2)  # 5 steps ** 6 vars * 3 steps x 2 vars
    while count < total_exp:

        qsub_fn = "qsub_scripts/run_sensitivity_exp_%s_%s.sh" % (pft, node)
        make_qsub_file(qsub_fn, pft, node)

        qs_cmd = "qsub %s" % (qsub_fn)
        error = subprocess.call(qs_cmd, shell=True)
        if error is 1:
            print("Job failed to submit")

        node += 1
        count += 7000
