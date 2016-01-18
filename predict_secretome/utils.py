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

    #parser.add_argument('--transporter_threshold', '-t',
    #                    action='store',
    #                    dest='trans',
    #                    default=0,
    #                    type=int,
    #                    help='Minimum number of tm domains in mature sequence'
    #                         'required to consider a protein as a transporter')

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




def detect_and_output_transporters(all_sequences_fp,
                                   rename_mappings,
                                   tmp_dir,
                                   run_name,
                                   path='dependencies/bin',
                                   mature_seqs=None,
                                   verbose=False,
                                   tm_threshold=2):
    """
    Find sequences which num tm domains > tm_threshold in
    a) mature sequences if they have a signal peptide
    b) full sequences for all other proteins
    input: formatted_fasta_fp- all fasta sequences
           rename_mappings- mappings to original seq names
           mature_sequences_fp- the mature sequences of those seqs with sigpeps
           tmp_dir- temporary file output directory
           tm_threshold- integer minimum number of tm domains to search for
    output: None
    """
    # need to combine sequences that have signalpeptides and thus mature seqs
    # with the full sequences of all the sequences without them
    if mature_seqs:

        all_input_seqs_fh = open(all_sequences_fp, 'r')

        with open(mature_seqs, 'r') as mature_seqs_fh:
            mature_seqs = {record.id: record for \
                            record in SeqIO.parse(mature_seqs_fh, 'fasta')}

        transporter_search_seqs_fp = os.path.join(tmp_dir,
                                                  'mature_if_sigpep_full_otherwise.fas')

        transporter_search_seqs_fh = open(transporter_search_seqs_fp, 'w')

        # if record is in the mature seqs (i.e. has signal peptide) then write
        # that otherwise write the full sequence
        for record in SeqIO.parse(all_input_seqs_fh, 'fasta'):
            if record.id in mature_seqs.keys():
                SeqIO.write(mature_seqs[record.id], transporter_search_seqs_fh, 'fasta')
            else:
                SeqIO.write(record, transporter_search_seqs_fh, 'fasta')

        all_input_seqs_fh.close()
        transporter_search_seqs_fh.close()
    else:
        transporter_search_seqs_fp = all_sequences_fp


    tmhmm_fp = os.path.join(path, 'tmhmm')
    tmhmm_cmd = "{0} {1}".format(tmhmm_fp,
                                 transporter_search_seqs_fp)

    print_verbose("\n##Search for Putative Transporters##", v_flag=verbose)
    tmhmm_output = subprocess.check_output(tmhmm_cmd.split())
    tmhmm_output = tmhmm_output.decode('ascii').split('\n')
    print_verbose("Search complete", v_flag=verbose)

    print_verbose("Parsing results", v_flag=verbose)
    # Parse tmhmm raw output and write acc with more than
    # in non signal peptide sequence to output_file
    putative_transporter_acc = []
    for line in tmhmm_output[:-1]:
        line = line.split('\t')
        predhel = int(line[4].lstrip('PredHel='))
        if predhel >= tm_threshold:
            putative_transporter_acc.append(line[0])

    if putative_transporter_acc:
        transporter_out_fh = open(run_name+'_predicted_transporters.fas', 'w')
        transporter_out = SeqIO.FastaIO.FastaWriter(transporter_out_fh)
        transporter_out.write_header()

        all_input_seqs_fh = open(all_sequences_fp, 'r')
        for record in SeqIO.parse(all_input_seqs_fh, 'fasta'):
            if record.id in putative_transporter_acc:
                record.id = rename_mappings[record.id]
                record.description = ''
                record.name = ''
                transporter_out.write_record(record)

        transporter_out.write_footer()
        transporter_out_fh.close()
        all_input_seqs_fh.close()



def secretome(accessions_with_sig_pep,
              accesions_no_tm_in_mature_seq,
              secreted_accessions,
              extracellular_accessions,
              conservative=True,
              verbose=False):
    """
    Combined predicted accessions
    input: formatted_fasta_fp - formatted fasta file of all input seqs
           fungi_flag - boolean to use fungi settings or not for targetp
           path - path to dependencies
           verbose
    output: extracellular_accessions - list of accessions with signalpeps
                                        predicted to be extracellular

    """

    print_verbose("\n##Combining predictions##", v_flag=verbose)

    sig_peptides = set(accessions_with_sig_pep)
    no_tm_domains = set(accesions_no_tm_in_mature_seq)
    secreted = set(secreted_accessions)
    extracellular = set(extracellular_accessions)


    # if conservative only get those accessions predicted as
    # A) having a signal peptide (signalp)
    # B) not having any TM domains outside of this sigpep
    #     (tmhmm)
    # C) a signal peptide predicted as being for secretion
    #     (targetp)
    # D) a sequence predicted as belonging to the extracellular
    #     compartment (wolFPSort)
    # if permissive get all accessions that belong to any of these
    # categories
    if conservative:
        out_flag = 'conservative'
        predicted_acc_list = set.intersection(sig_peptides,
                                              no_tm_domains,
                                              secreted,
                                              extracellular)
    else:
        out_flag = 'permissive'
        # maybe remove no_tm_domains from this
        # as plenty of things don't have tm domains
        # that aren't secreted
        predicted_acc_list = set.union(sig_peptides,
                                       no_tm_domains,
                                       secreted,
                                       extracellular)

    if len(predicted_acc_list) is 0:
        warnings.warn("No secreted proteins found using {0} setting".format(out_flag))

    return predicted_acc_list



