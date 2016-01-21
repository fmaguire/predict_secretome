#!/usr/bin/env python
"""
Input multi-protein fasta and generate a conservative or permissive
predicted secretome as output.
"""
from __future__ import print_function
import sys
import subprocess
import shutil
import warnings
import argparse
import os
import hashlib
from Bio import SeqIO


def get_parser():
    """
    Generate parser obj with required options
    """

    parser = argparse.ArgumentParser(description='Take in a fasta and '
                                                 'identify the predicted '
                                                 'secretome')

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument('--fasta', '-f',
                        action='store',
                        type=str,
                        dest='input_file',
                        help='Input fasta file for secretome prediction')

    group.add_argument('--check',
                        action='store_true',
                        dest='check',
                        default=False,
                        help='Check dependencies of %(prog)s then quit '
                             '(default: %(default)s)')

    parser.add_argument('--run_name', '-n',
                        action='store',
                        type=str,
                        dest='run_name',
                        default=False,
                        help='Prefix for output (default is input filename)')

    euk_type_group = parser.add_mutually_exclusive_group(required=True)

    euk_type_group.add_argument('--plant',
                                action='store_true',
                                dest='plant',
                                help='Sequences are from plants')

    euk_type_group.add_argument('--fungi',
                                action='store_true',
                                dest='fungi',
                                help='Sequences are from fungi')

    euk_type_group.add_argument('--animal',
                                action='store_true',
                                dest='animal',
                                help='Sequences are from animal')


    #parser.add_argument('--verbose', '-v',
    #                    action='store_true',
    #                    dest='verbose',
    #                    default=False,
    #                    help='Print verbose output '
    #                         '(default: %(default)s)')

    parser.add_argument('--force',
                        action='store_true',
                        dest='force',
                        default=False,
                        help='Force overwriting of existing directory'
                             '(default: %(default)s)')

    parser.add_argument('--permissive', '-p',
                        action='store_true',
                        dest='permissive',
                        default=False,
                        help='Get permissive secretome prediction '
                             '(default: %(default)s)')

    parser.add_argument('--transporter_threshold', '-t',
                        action='store',
                        dest='transporters',
                        default=False,
                        type=int,
                        help='Minimum number of tm domains in mature sequence'
                             'required to consider a protein as a transporter')

    return parser

def which(bin_path, program):
    """
    Check dependency exists and is executable
    """
    fpath = os.path.join(bin_path, program)
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def print_verbose(string, v_flag=False):
    """
    Print string only if v_flag is true

    input: string - string to print
           v_flag - bool to say whether to print
    output: none
    """

    if v_flag:
        print(string)


def check_dependencies(path='predict_secretome/dependencies/bin',
                       check_run=False,
                       verbose=False):
    """
    Check the dependencies are in the bin_path and are executable
    then output a summary if doing a dry 'check run' or in verbose mode
    if doing a dry 'check run' exit after

    input: bin_path - path to dependencies
           argv - args from command line
    output: True if it doesn't exit beforehand
    """

    # if we are running a check run then we want verbose output for this
    # regardless of other settings
    if check_run:
        verbose = True

    dependency_list = ["signalp", "tmhmm", "targetp",
                       "runWolfPsortSummary"]

    check_execs = {dependency: which(path, dependency) for dependency in\
                        dependency_list}

    print_verbose("\n##Checking Depdencies##", v_flag=verbose)
    for dependency in check_execs:
        print_verbose("{0} is present? {1}".format(dependency,
                                                   check_execs[dependency]),
                      v_flag=verbose)

    # if any dependency is missing quit
    if not all(check_execs.values()):
        print("Dependency absent")
        sys.exit(0)

    # if running in check mode quite after checking
    if check_run:
        print("Finished dependency check run")
        sys.exit(0)

    return True

def write_seqs_from_accessions(acc_list,
                               input_file,
                               output_file,
                               append=False):
    """
    Parse an input fasta file and write sequences
    to the output if their accession is in the acc_list
    """
    seqs_written = 0

    if append:
        write_flag = 'a'
    else:
        write_flag = 'w'

    with open(output_file, write_flag) as out_fh:
        fasta_out = SeqIO.FastaIO.FastaWriter(out_fh, wrap=None)
        fasta_out.write_header()
        for record in SeqIO.parse(input_file, 'fasta'):
            if record.description in acc_list:
                seqs_written += 1
                fasta_out.write_record(record)
        if seqs_written > 0:
            fasta_out.write_footer()

