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
Finlay Maguire 2014
Adapted from a similar script by Ian Misner 

Input multi-protein fasta and generate a conservative or permissive 
predicted secretome as output.  
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
                        dest='nocleanup',
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

    parser.add_argument('--permissive', '-p',
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
                       "runWolfPsortSummary"]

    check_execs = {dependency: which(path, dependency) for dependency in\
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

    return True

def format_fasta(input_file,
                 tmp_dir,
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

    formatted_fasta = os.path.join(tmp_dir, "formatted_input.fasta")

    # as some dependencies (tmhmm) don't support accessions over 20 chars
    # replace accessions with a 19-char truncated hash of original accession
    # keeping track in a dict where key is the hash and val is original accession 
    in_handle = open(input_file, "r")
    out_handle = open(formatted_fasta, 'w')

    # SeqIO doesn't expose the wrap attribute so need to call FastaIO directly
    # to use single line format
    fasta_out = FastaIO.FastaWriter(out_handle, wrap=None)
    fasta_out.write_header()
    rename_mappings = {}

    for record in SeqIO.parse(in_handle, 'fasta'):
        truncated_md5 = hashlib.md5(\
              record.description.encode('utf-8')).hexdigest()[:19]
              # hash collision calculation fun time:
              # md5 is 128bits (hexdigest 32 chars) so a 19 character truncation 
              # is 76 bits.  Therefore can calculate collision rate using birthday 
              # attack 
              # 

        rename_mappings.update({truncated_md5: record.description})

        record.id=truncated_md5
        record.name=''
        record.description=''

        fasta_out.write_record(record)
   
    fasta_out.write_footer()

    in_handle.close()
    out_handle.close()

    print_verbose("Input fasta formatted: {0}".format(input_file),
                  v_flag=verbose)

    return rename_mappings, formatted_fasta


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
    output: seqs_with_sigpep_removed - filename containg seqs without their sigpeps 
            acc_with_sigpeps - filename containing acc of seqs with sig peptides
            full_sequences_with_sigpep - filname containing seqs with sig peptides including
                                         the sigpep
    """


    seqs_sigpep_removed = os.path.join(tmp_dir, 'signalp_removed_sigpep.fasta')

    sigp = os.path.join(path, 'signalp')
    sigp_cmd= "{0} -t euk -f short -m {1} {2}".format(sigp,
                                                      seqs_sigpep_removed,
                                                      input_file)

    print_verbose("Running SignalP", v_flag=verbose)
    with open(os.devnull) as null:
        sigp_retcode = subprocess.call(sigp_cmd.split(), stdout=null)
    print_verbose("SignalP Complete", v_flag=verbose)

    print_verbose("Creating SignalP protein list", v_flag=verbose)
    # get list of all accessions with signal peptides
    accessions = []
    with open(seqs_sigpep_removed, 'r') as fh:
        for line in fh.readlines():
            if line.startswith('>'):
                strip_line = line.lstrip('>')
                accessions.append(strip_line.split(' ; ')[0])

    acc_with_sigpeps = os.path.join(tmp_dir, 'signalp_acc_with_sigpep')
    with open(acc_with_sigpeps, 'w') as acc_w_sigpep_fh:
        for acc in accessions:
            acc_w_sigpep_fh.write(acc + '\n') 

 
    # use this list to get sequences with signal peptides from
    # formatted input sequences
    full_sequences_with_sigpep = os.path.join(tmp_dir, 'signalp_seqs_with_sigpep_including_sigpeps.fasta')
    with open(full_sequences_with_sigpep, 'w') as sigpep_seqs_fh:
        for record in SeqIO.parse(input_file, 'fasta'):
            if record.description in accessions:
                SeqIO.write(record, sigpep_seqs_fh, 'fasta')

    return seqs_sigpep_removed, acc_with_sigpeps, full_sequences_with_sigpep


def tmhmm(seqs_with_sigpep_removed,
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

    tmhmm = os.path.join(path, 'tmhmm')
    tmhmm_cmd= "{0} {1}".format(tmhmm,
                                seqs_with_sigpep_removed)
    

    print_verbose("Running TMHMM on mature signalp sequences only", v_flag=verbose)
    tmhmm_output = subprocess.check_output(tmhmm_cmd.split())
    tmhmm_output = tmhmm_output.decode('ascii').split('\n')
    print_verbose("TMHMM complete", v_flag=verbose)

    print_verbose("Identifying sequences without tm regions.", v_flag=verbose)
    # Parse tmhmm raw output and write acc without tm domains
    # in non signal peptide sequence to output_file
    parsed_output = [line.split('\t')[0] \
                        for line in tmhmm_output if "\tPredHel=0\t" in line]

    return parsed_output


def targetp(full_seqs_with_sigpeps,
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
    targetp_cmd = "{0} {1} {2}".format(targetp, targetp_flag, full_seqs_with_sigpeps)

    print_verbose("Running TargetP on sequences "
                  "identified by SignalP as having a sig_pep", v_flag=verbose)
    targetp_output = subprocess.check_output(targetp_cmd.split())
    targetp_output = targetp_output.decode('ascii').split('\n')
    print_verbose("TargetP complete", v_flag=verbose)


    print_verbose("Identifying sequences that are secreted", v_flag=verbose)
    # remove header and tail cruft in targetp output
    # and get those that have S as top predicted target 
    parsed_output = [line.split(' ')[0] \
                          for line in targetp_output[8:-3] if "S" in line]

    return parsed_output


def wolfpsort(input_file,
              tmp_dir,
              path='dependencies/bin',
              fungi_flag=True,
              verbose=False):
    """
    runs wolfPsort with the fungi setting. change as necessary for you usage.
    """
    wolfpsort = os.path.join(path, 'runWolfPsortSummary')

    if fungi_flag:
        wps_opt = "fungi"

    wolfpsort_cmd = "{0} {1} < {2}".format(wolfpsort, wps_opt, input_file)
    raw_output = os.path.join(tmp_dir, 'wolfpsort_raw_output')

    print_verbose("Running WoLFPSORT", v_flag=verbose)
    wolfpsort_output = subprocess.check_output(wolfpsort_cmd, shell=True)
    wolfpsort_output = wolfpsort_output.decode('ascii').split('\n')
    print_verbose("WoLFPSORT complete", v_flag=verbose)


    print_verbose("Parsing WoLFPSort output", v_flag=verbose)
    # Removes header from output 
    parsed_output = [line.split(' ')[0] \
                        for line in wolfpsort_output[1:-1] if "extr" in line.split()[1]]

    return parsed_output


def strip_and_read_to_set(fpath):

    with open(fpath, 'r') as fh:
        output = set(line.strip() for line in fh.readlines())

    return output


def secretome(signalp_pred_with_sigpep_fn,
              tmhmm_pred_with_no_tm_fn,
              targetp_pred_secreted_fn,
              wolfpsort_pred_ec_fn,
              tmp_dir,
              conservative=True,
              verbose=False):
    
    sig_peptides  = strip_and_read_to_set(signalp_pred_with_sigpep_fn)
    no_tm_domains = strip_and_read_to_set(tmhmm_pred_with_no_tm_fn)
    secreted      = strip_and_read_to_set(targetp_pred_secreted_fn)
    extracellular = strip_and_read_to_set(wolfpsort_pred_ec_fn)

    if conservative:
        out_flag = 'conservative'
        predicted_acc_list = set.intersection(sig_peptides, 
                                              no_tm_domains, 
                                              secreted, 
                                              extracellular)
    else:
        out_flag = 'permissive'
        predicted_acc_list = set.union(sig_peptides, 
                                       no_tm_domains, 
                                       secreted, 
                                       extracellular)

    prediction_output_fpath = os.path.join(tmp_dir, 
                                           "{0}_predicted_secretome_"
                                           "accessions.txt".format(out_flag))

    if len(predicted_acc_list) is 0:  
        warnings.warn("No secreted proteins found using {0} setting".format(out_flag))

    with open(prediction_output_fpath, 'w') as out_fh:
        for acc in predicted_acc_list:
            out_fh.write(acc + '\n')


    return predicted_acc_list


def generate_output(formatted_fasta,
                    predicted_secretome_acc_list, 
                    rename_mappings,
                    run_name,
                    conservative=True,
                    verbose=False):

    in_handle = open(formatted_fasta, 'r')

    if conservative:
        out_flag = 'conservative'
    else:
        out_flag = 'permissive'

    output = "{0}_{1}_predicted_secretome.fasta".format(run_name, out_flag)
    out_handle = open(output, 'w')

    fasta_out = FastaIO.FastaWriter(out_handle, wrap=None)
    fasta_out.write_header()

    print_verbose("Retreving Secretome fasta", v_flag=verbose)
    for record in SeqIO.parse(in_handle, 'fasta'):
        if record.description in predicted_secretome_acc_list:
            record.id=rename_mappings[record.id]
            record.name=''
            record.description=''
            fasta_out.write_record(record)

    in_handle.close()
    fasta_out.write_footer()
    out_handle.close()
    print_verbose("Secretome identification Complete", v_flag=verbose)
