#!/usr/bin/env python3

import shutil
import sys
import subprocess
import warnings
import argparse
import os
import glob
import hashlib
from Bio import SeqIO
from Bio.SeqIO import FastaIO

"""
This is a script desiged to take a multiple protien fasta as input and provide
the predicted secretome as output. The program requires that signalp, tmhmm,
targetp, chlorop, faSomeRecords, wolfpsort, and fasta_formatter be in the
users PATH. See those programs documentation for installation instructions.
Settings are for fungi change inputs for programs as necessary (targetp,
and wolfpsort)
"""



def get_parser():

    parser = argparse.ArgumentParser(description='Take in a fasta and '
                                                 'identify the predicted '
                                                 'secretome')
    parser.add_argument('--fasta', '-f',
                        action='store',
                        type=str,
                        dest='input_file',
                        required=True,
                        help='Input fasta file for secretome prediction '
                             '(REQUIRED)')

    parser.add_argument('--run_name', '-n',
                        action='store',
                        type=str,
                        dest='run_name',
                        default=False,
                        help='Prefix for output (default is input filename)')

    parser.add_argument('--no_cleanup', '-j',
                        action='store_true',
                        dest='no_cleanup'
                        default=False,
                        help="Don't remove intermediate outputs")

    parser.add_argument('--check', '-x',
                        action='store_true',
                        dest='check',
                        default=False,
                        help='Check dependencies of %(prog)s then quit '
                             '(default: %(default)s)')

    parser.add_argument('--verbose', '-v',
                        action='store_true',
                        dest='verbose',
                        default=False,
                        help='Print verbose output '
                             '(default: %(default)s)')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--conservative', '-c',
                       action='store_true',
                       dest='conservative',
                       default=False,
                       help='Get conservative secretome prediction '
                            '(default: %(default)s)')

    group.add_argument('--permissive', '-p',
                       action='store_true',
                       dest='permissive',
                       default=False,
                       help='Get permissive secretome prediction '
                            '(default: %(default)s)')

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


def check_dependencies(path='dependencies/bin',
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
                       "chlorop", "runWolfPsortSummary",
                       "fasta_formatter"]

    check_execs = {dependency: which(bin_path, dependency) for dependency in\
                        dependency_list}

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


def format_fasta(input_file,
                 output_file,
                 path='dependencies/bin',
                 verbose=False):
    """
    Format input fasta and write to output:
    Specifically replace accessions with hash of name to guarantee under 20 chars
    input: input_file - filename
           output_file - filename
           path - path to dependencies bins
           verbose 
    output: rename_mappings - a dict of {new_acc: original_acc}
    """

    print_verbose("Formatting input fasta: {0}".format(input_file),
                  v_flag=verbose)


    # as some dependencies (tmhmm) don't support accessions over 20 chars
    # replace accessions with a 19-char truncated hash of original accession
    # keeping track in a dict where key is the hash and val is original accession 
    in_handle = open(input_file, "rU")
    out_handle = open(output_file, 'w')

    # SeqIO doesn't expose the wrap attribute so need to call FastaIO directly
    # to use single line format
    fasta_out = FastaIO.FastaWriter(out_handle, wrap=None)

    rename_mappings = {}

    for record in SeqIO.parse(in_handle, 'fasta'):
        truncated_md5 = hashlib.md5(\
              record.description.encode('utf-8')).hexdigest()[:19]

        rename_mappings.update({truncated_md5: record.description})

        record.id=truncated_md5
        record.name=''
        record.description=''

        fasta_out.write_record(record)

    in_handle.close()
    out_handle.close()

    print_verbose("Input fasta formatted: {0}".format(input_file),
                  v_flag=verbose)

    return rename_mappings


def signalp(input_file,
            tmp_dir,
            path='dependencies/bin',
            verbose=False):
    """
    Use signalp to identify sequences with signal peptides then extract
    the mature sequences with these peptides cleaved
    input: input_file - filename (reformatted_fasta)
           tmp_dir - path to temporary intermiedate output
           path - path to dependency bins
           verbose
    output: seqs_with_sigpep_removed - filename 
    """


    seqs_with_sigpep_removed = os.path.join(tmp_dir, 'signalp_removed_sigpep.fasta')

    sigp = os.path.join(path, 'signalp')
    sigp_cmd= "{0} -t euk -f short -m {0} {1}".format(sigp,
                                                      seqs_with_sigpep_removed,
                                                      input_file)

    print_verbose("Running SignalP", v_flag=verbose)
    with open(os.devnull) as null:
        sigp_retcode = subprocess.call(command.split(), stdout=null)
    print_verbose("SignalP Complete", v_flag=verbose)

    print_verbose("Creating SignalP protein list", v_flag=verbose)
    # get list of all accessions with signal peptides
    accessions = []
    with open(seqs_without_sigpep, 'rU') as fh:
        for line in fh.readlines()
            if line.starswith('>'):
                line.lstrip('>')
                accessions.append(line.split(' ; ')[0])

    # use this list to get sequences with signal peptides from
    # formatted input sequences
    full_sequences_with_sigpep = output_file
    for record in SeqIO.parse(input_file, 'fasta'):
        if record.description in accessions:
            SeqIO.write(record, full_sequences_with_sigpep, 'fasta')

    return seqs_with_sigpep_removed, full_sequences_with_sigpep


def tmhmm(input_file,
          tmp_dir,
          path='dependencies/bin',
          verbose=False):
    """
    Run tmhmm on the sequences with signal peptides to check for tm
    domains
    This only checks mature sequences so ignores TM in the signal
    peptide
    input: input_file - fasta file without signal peptides
           tmp_dir - working temporary output directory
           path - path to dependencies
           verbose
    output: acc_w_no_tm_domains - file containing list of accessions without tm
                                  domains
    """

    acc_w_no_tm_domains = os.path.join(tmp_dir, 'tmhmm_acc_with_no_tm_domain.txt')

    tmhmm = os.path.join(path, 'tmhmm')
    tmhmm_cmd= "{0} {1}".format(tmhmm,
                                input_file)
    raw_output = os.path.join(tmp_dir, acc_w_no_tm_domains + '_tmp')

    print_verbose("Running TMHMM on mature signalp sequences only", v_flag=verbose)
    with open(raw_output, 'w') as fh:
        tmhmm_retcode = subprocess.call(tmhmm_cmd.split(), stdout=fh)
    print_verbose("TMHMM complete", v_flag=verbose)

    print_verbose("Identifying sequences without tm regions.", v_flag=verbose)
    # Parse tmhmm raw output and write acc without tm domains
    # in non signal peptide sequence to output_file
    with open(raw_output, 'r') as raw_fh:
        with open(acc_with_no_tm_domains, 'w') as out_fh:
            for line in raw_fh:
                if "\tPredHel=0\t" in line:
                    out_fh.write(line.split('\t')[0] + '\n')


    os.remove(raw_output)

    return acc_w_no_tm_domains


def targetp(INPUTS SIGNALP_PASS.fasta
            tmp_dir,
            path='dependencies/bin',
            plant=False,
            verbose=False):
    """
    Runs targetp to identify seqs with 'secreted' desigation in signal peptide
    """
    if plant:
        targetp_flag = "-P"
    else:
        targetp_flag = "-N"

    targetp = os.path.join(path, 'targetp')

    targetp_cmd = "{0} {1} {2}".format(targetp, targetp_flag, input_file)

    raw_targetp_out = os.path.join(tmp_dir, 'raw_targetp_output')

    print_verbose("Running TargetP on SignalP pass seqeunces only", v_flag=verbose)
    with open(raw_targetp_out, 'w') as fh:
        targetp_cmd = subprocess.call(command.split(), stdout=fh)
    print_verbose("TargetP complete", v_flag=verbose)

    print_verbose("Identifying sequences that are secreted", v_flag=verbose)
    targetp_secreted_acc = os.path.join(tmp_dir, 'targetp_secreted_acc.txt')

    # remove header and tail cruft in targetp output
    with open(raw_targetp_out, 'r') as raw_fh:
        raw_output = [line for line in raw_fh.readlines()]
        raw_output = raw_output[8:-2]

    with open(targetp_secreted_acc, 'w') as out_fh:
        for line in raw_output:
            if "S" in line:
                out_fh.write(line.split(' ')[0]+'\n')

    os.remove(raw_targetp_out)

    return targetp_secreted_acc


def wolfpsort(input_file,
              tmp_dir,
              path='dependencies/bin'
              fungi_flag=True,
              verbose=False):
    """
    runs wolfPsort with the fungi setting. change as necessary for you usage.
    """
    wolfpsort = os.path.join(path, 'runWolPsortSummary')

    if fungi_flag:
        wps_opt = "fungi"

    wolfpsort_cmd = "{0} {1} < {2}".format(wolfpsort, wps_opt, input_file)
    raw_output = os.path.join(tmp_dir, 'wolfpsort_raw_output')

    print_verbose("Running WoLFPSORT", v_flag=verbose)
    with open(raw_output, 'w') as raw_out_fh:
        wolfpsort_retcode = subprocess.call(wolfpsort_cmd.split(), stdout = raw_out_fh, shell=True)
    print_verbose("WoLFPSORT complete", v_flag=verbose)

    # Removes header from output file
    with open(raw_output, 'r') as raw_out_fh:
        raw_output = [line for lines in raw_out_fh.readlines()]
        raw_output = raw_output[1:]

    wolfpsort_ec_acc = os.path.join(tmp_dir, 'wolfpsort_extracellular_acc.txt')

    # Places fastas with extracellualr location into good list
    with open(wolfpsort_ec_acc, 'w') as out_fh:
        for line in raw_output:
            if line.split()[1] == "extr":
                out_fh.write(line.split(' ')[0] + '\n')

    os.remove(raw_output)

    return wolfpsort_ec_acc


def secretome(signalp_acc_with_sigpep,
              tmhmm_acc_with_no_tm,
              targetp_secreted_acc,
              wolfpsort_ec_acc,
              tmp_dir,
              conservative=True,
              verbose=False):

    sig_peptides  = set(acc.strip() for acc in open(signalp_acc_with_sigpep))
    no_tm_domains = set(acc.strip() for acc in open(tmhmm_acc_with_no_tm))
    secreted      = set(acc.strip() for acc in open(targetp_secreted_acc))
    extracellular = set(acc.strip() for acc in open(wolfpsort_extracellular_acc))

    if conservative:
        out_flag = 'conservative'
        predicted_acc = set.intersection(sig_peptides, 
                                         no_tm_domains, 
                                         secreted, 
                                         extracellular)
    else:
        out_flag = 'permissive'
        predicted_acc = set.union(sig_peptides, 
                                  no_tm_domains, 
                                  secreted, 
                                  extracellular)

    predicted_secretome_acc = os.path.join(tmp_dir, 
                                           "{0}_predicted_secretome_"
                                           "accessions.txt".format(out_flag))


    with open(predicted_secretome_acc, 'w') as out_fh:
        for acc in predicted_acc:
            print(acc)
            out_fh.write(acc + '\n')


def assemble_secretome(formatted_fasta,
                       pred_secretome_acc, 
                       output_file,
                       rename_mappings, 
                       verbose=False):

    in_handle = formatted_fasta
    out_handle = output_file

    print_verbose("Retreving Secretome fasta", v_flag=verbose)
    for record in SeqIO.parse(in_handle, 'fasta'):
        if record.description in secretome_list:
            record.id=rename_mappings[truncated_md5]
            record.name=''
            record.description=''
            SeqIO.write(record, out_handle, 'fasta')

    in_handle.close()
    out_handle.close()
    print_verbose("Secretome identification Complete", v_flag=verbose)
