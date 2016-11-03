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
                        default=False,
                        help='Print verbose output '
                             '(default: %(default)s)')

    parser.add_argument('--permissive', '-p',
                        action='store_true',
                        dest='permissive',
                        default=False,
                        help='Get permissive secretome prediction '
                             '(default: %(default)s)')

    parser.add_argument('--transporter_threshold', '-t',
                        action='store',
                        dest='trans',
                        default=0,
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

def format_fasta(input_file_fp,
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

    print_verbose("\n##Formatting input fasta: {0}##".format(input_file_fp),
                  v_flag=verbose)

    formatted_fasta_fp = os.path.join(tmp_dir, "formatted_input.fasta")

    # as some dependencies (tmhmm) don't support accessions over 20 chars
    # replace accessions with a 19-char truncated hash of original accession
    # keeping track in a dict where key is the hash and val is original accession
    raw_fasta_in_fh = open(input_file_fp, "r")
    formatted_out_fh = open(formatted_fasta_fp, 'w')

    # SeqIO doesn't expose the wrap attribute so need to call FastaIO directly
    # to use single line format
    fasta_out = SeqIO.FastaIO.FastaWriter(formatted_out_fh, wrap=None)
    fasta_out.write_header()
    rename_mappings = {}

    for record in SeqIO.parse(raw_fasta_in_fh, 'fasta'):
        truncated_md5 = hashlib.md5(\
              record.description.encode('utf-8')).hexdigest()[:19]
              # hash collision calculation fun time:
              # md5 is 128bits (hexdigest 32 chars) so a 19 character truncation
              # is 76 bits.  Therefore can calculate collision rate using birthday
              # attack
              #

        rename_mappings.update({truncated_md5: record.description})

        record.id = truncated_md5
        record.name = ''
        record.description = ''

        fasta_out.write_record(record)

    fasta_out.write_footer()

    raw_fasta_in_fh.close()
    formatted_out_fh.close()

    print_verbose("Input fasta formatted: {0}".format(formatted_fasta_fp),
                  v_flag=verbose)

    return rename_mappings, formatted_fasta_fp


def signalp(input_file,
            tmp_dir,
            rename_mappings,
            run_name,
            trans=0,
            path='dependencies/bin',
            verbose=False):
    """
    Use signalp to identify sequences with signal peptides and create
    mature sequences with these signal peptides cleaved out
    input: input_file - filename (reformatted_fasta)
           tmp_dir - path to temporary intermiedate output
           path - path to dependency bins
           verbose
    output: mature_seqs_fp - filename containg mature seqs (sigpeps cleaved)
            accessions_with_sig_pep - list of accessions with sigpeps
            full_sequences_with_sigpep_fp - filname containing full seqs of seqs
                                            id'd as having a sigpep
    """


    # mature seqs are those with signal peptides removed
    mature_seqs_fp = os.path.join(tmp_dir, 'signalp_mature_seqs.fasta')

    sigp = os.path.join(path, 'signalp')
    sigp_cmd = "{0} -t euk -f short -m {1} {2}".format(sigp,
                                                       mature_seqs_fp,
                                                       input_file)


    # we only care about the mature sequences output as they are the sequences
    # predicted as having signal peptides (which have then been subsequently
    # trimmed)

    print_verbose("\n##Detecing signal peptides##", v_flag=verbose)
    sigp_proc = subprocess.Popen(sigp_cmd.split(),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
    _, stderr = sigp_proc.communicate()

    print_verbose("Search Complete", v_flag=verbose)

    if '# No sequences predicted with a signal peptide' in stderr.decode('ascii'):
        print("\nNo signal peptides detected therefore no secretome is predicted\n")
        # if there are no signal peptides can still calculate and return tm if specified
        if trans:
            detect_and_output_transporters(input_file,
                                           rename_mappings,
                                           tmp_dir,
                                           run_name,
                                           path=path,
                                           verbose=verbose,
                                           tm_threshold=trans)
        shutil.rmtree(tmp_dir)
        sys.exit()


    print_verbose("Compiling accessions with signal peptides", v_flag=verbose)
    # get list of all accessions with signal peptides from mature
    # sequences file
    accessions_with_sig_pep = []
    with open(mature_seqs_fp, 'r') as mature_fh:
        for line in mature_fh.readlines():
            if line.startswith('>'):
                strip_line = line.lstrip('>')
                accessions_with_sig_pep.append(strip_line.split(' ; ')[0])


    # use this list to get full sequences (i.e. not mature) with signal peptides from
    # formatted input sequences
    full_sequences_with_sigpep_fp = os.path.join(tmp_dir,
                                                 'signalp_full_seqs_with_sigpep.fasta')

    print_verbose("Assembling full sequences with signal peptides", v_flag=verbose)
    with open(full_sequences_with_sigpep_fp, 'w') as sigpep_seqs_fh:
        fasta_out = SeqIO.FastaIO.FastaWriter(sigpep_seqs_fh, wrap=None)
        fasta_out.write_header()
        for record in SeqIO.parse(input_file, 'fasta'):
            if record.description in accessions_with_sig_pep:
                fasta_out.write_record(record)

        fasta_out.write_footer()


    return mature_seqs_fp, accessions_with_sig_pep, full_sequences_with_sigpep_fp


def tmhmm(mature_seqs_fp,
          path='dependencies/bin',
          verbose=False):
    """
    Run tmhmm on the mature sequences from signalp to ensure they
    don't have transmembrane domains outwith of their signal peptide
    input: mature_seqs_fp - fasta file without signal peptides
           tmp_dir - working temporary output directory
           path - path to dependencies
           verbose
    output: acc_without_tm_in_mature_seq - list of accessions
                                           without TM domains
                                           in their mature
                                           sequences
    """

    tmhmm_fp = os.path.join(path, 'tmhmm')
    tmhmm_cmd = "{0} {1}".format(tmhmm_fp,
                                 mature_seqs_fp)


    print_verbose("\n##Search for TM domains in mature seqs##", v_flag=verbose)
    tmhmm_output = subprocess.check_output(tmhmm_cmd.split())
    tmhmm_output = tmhmm_output.decode('ascii').split('\n')
    print_verbose("Search complete", v_flag=verbose)

    print_verbose("Parsing results", v_flag=verbose)
    # Parse tmhmm raw output and write acc without tm domains
    # in non signal peptide sequence to output_file
    acc_without_tm_in_mature_seq = [line.split('\t')[0] \
                                      for line in tmhmm_output \
                                        if "\tPredHel=0\t" in line]
    return acc_without_tm_in_mature_seq


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


def targetp(full_sequences_with_sigpep_fp,
            path='dependencies/bin',
            plant=False,
            verbose=False):
    """
    Runs targetp on full seqs of those identified as having
    signal peptides by signalp to identify those designated
    as targeted 'secreted'
    input: mature_seqs_fp - fasta file without signal peptides
           tmp_dir - working temporary output directory
           plant - boolean to use plant settings or not for targetp
           path - path to dependencies
           verbose
    /utput: secreted_accessions - list of accessions with signalpeps
                                  predicted to be secreted

    """

    if plant:
        targetp_flag = "-P"
    else:
        targetp_flag = "-N"

    # the fortran 'How' ANN classifier dependency bundled with targetp
    # can't handle paths longer than 80 chars so as a hack fix I'm creating
    # a symlink to /home/user and removing it after execution

    home_targetp = os.path.join(os.path.expanduser('~'),
                                'targetp-1.1')

    if not os.path.exists(home_targetp):
        shutil.copytree('dependencies/targetp-1.1', home_targetp)

    targetp_path = os.path.join(path, 'targetp')
    targetp_cmd = "{0} {1} {2}".format(targetp_path,
                                       targetp_flag,
                                       full_sequences_with_sigpep_fp)

    print_verbose("\n##Identifying sequences with 'secreted' sigpeps##",
                  v_flag=verbose)
    targetp_output = subprocess.check_output(targetp_cmd.split())
    targetp_output = targetp_output.decode('ascii').split('\n')
    print_verbose("Search complete", v_flag=verbose)

    # remove symlink or file
    try:
        shutil.rmtree(home_targetp)
    except OSError:
        os.remove(home_targetp)

    print_verbose("Parsing results", v_flag=verbose)
    # remove header and tail cruft in targetp output
    # and get those that have S as top predicted target
    secreted_accessions = [line.split(' ')[0] \
                           for line in targetp_output[8:-3] if "S" in line]

    return secreted_accessions


def wolfpsort(formatted_fasta_fp,
              path="dependencies/bin",
              fungi_flag=True,
              verbose=False):
    """
    Run wolfPsort on all formatted sequences using fungi setting
    to get a predicted list of 'extracellular' accessions
    input: formatted_fasta_fp - formatted fasta file of all input seqs
           tmp_dir - working temporary output directory
           fungi_flag - boolean to use fungi settings or not for targetp
           path - path to dependencies
           verbose
    output: extracellular_accessions - list of accessions with signalpeps
                                        predicted to be extracellular

    """

    wolfpsort_path = os.path.join(path, 'runWolfPsortSummary')

    if fungi_flag:
        wps_opt = "fungi"

    wolfpsort_cmd = "{0} {1} < {2}".format(wolfpsort_path, wps_opt, formatted_fasta_fp)

    print_verbose("\n##Identifying sequences belonging to 'extracellular' compartment##",
                  v_flag=verbose)
    wolfpsort_output = subprocess.check_output(wolfpsort_cmd, shell=True)
    wolfpsort_output = wolfpsort_output.decode('ascii').split('\n')
    print_verbose("Search complete", v_flag=verbose)


    print_verbose("Parsing results", v_flag=verbose)
    # Removes header from output
    extracellular_accessions = [line.split(' ')[0] \
                                  for line in wolfpsort_output[1:-1] \
                                       if "extr" in line.split()[1]]

    return extracellular_accessions


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


def generate_output(formatted_fasta_fp,
                    predicted_secretome_accessions,
                    rename_mappings,
                    run_name,
                    conservative=True,
                    verbose=False):
    """
    Generate predicted secretome output with original accessions
    """

    print_verbose("\n##Writing predicted secretome fasta file##", v_flag=verbose)

    in_handle = open(formatted_fasta_fp, 'r')

    if conservative:
        out_flag = 'conservative'
    else:
        out_flag = 'permissive'

    output = "{0}_{1}_predicted_secretome.fasta".format(run_name, out_flag)
    out_handle = open(output, 'w')

    fasta_out = SeqIO.FastaIO.FastaWriter(out_handle, wrap=None)
    fasta_out.write_header()

    print_verbose("Retreving Secretome fasta", v_flag=verbose)
    for record in SeqIO.parse(in_handle, 'fasta'):
        if record.description in predicted_secretome_accessions:
            record.id = rename_mappings[record.id]
            record.name = ''
            record.description = ''
            fasta_out.write_record(record)

    in_handle.close()
    fasta_out.write_footer()
    out_handle.close()
    print_verbose("Secretome identification Complete", v_flag=verbose)



