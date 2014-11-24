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


def check_dependencies(bin_path, check_run=False, verbose=False):
    """
    Check the dependencies are in the bin_path and are executable
    output summary if checking dependencies or running in verbose mode
    quit after checking if running check mode

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


def format_fasta(input_file, output_file, argv):
    """
    Format input fasta and write to output
    Specifically replace accessions with hash of name to guarantee under 20 chars
    input: input_file - filename
           output_file - filename
           argv - args
    output: rename_mappings - a dict where keys are new accessions and vals
                              are the original accessions
    """

    print_verbose("Formatting input fasta: {0}".format(input_file),
                  v_flag=argv.verbose)

    intermediate_output = output_file + '_tmp'

    # as some dependencies (tmhmm) don't support accessions over 20 chars
    # make key-value lookup where key is 19-char truncated hash of name
    # and value is original accession for renaming later
    in_handle = open(input_file, "rU")
    out_handle = open(intermediate_output, 'w')

    rename_mappings = {}

    for record in SeqIO.parse(in_handle, 'fasta'):
        truncated_md5 = hashlib.md5(record.description.encode('utf-8')).hexdigest()[:19]

        rename_mappings.update({truncated_md5: record.description})
        record.id=truncated_md5
        record.name=''
        record.description=''

        SeqIO.write(record, out_handle, 'fasta')

    in_handle.close()
    out_handle.close()


    command = ("fasta_formatter -i {0} -o {1} -w 0".format(intermediate_output,
                                                           output_file))

    ret_code = subprocess.call(command.split())

    if ret_code is not 0:
        raise OSError('Process failure: {0}'.format(command))


    os.remove(intermediate_output)

    print_verbose("Input fasta formatted: {0}".format(input_file),
                  v_flag=argv.verbose)

    return rename_mappings


def signalp(argv, verbose=False):
    """
    use Signalp to identify sequences with signal peptides
    """
    command = "signalp -f short -m {0} {1} > {2}".format("removed_SigPep.fasta",
                                                         "singleline.fasta",
                                                         "signalpOUT.txt")

    print_verbose("Running SignalP", v_flag=verbose)
    signalpRUN = subprocess.call(command.split(), shell=True)
    print_verbose("SignalP Complete", v_flag=verbose)


    # Generate the list of sequences with siganal peptides using the mature sequences
    print_verbose("Creating SignalP protein list", v_flag=verbose)

    command2 = "fasta_formatter -i {0} -o {1} -t".format("removed_SigPep.fasta",
                                                         "removed_SigPep_tab.fasta.txt")

    tab = subprocess.call(command2.split(), stdout=file_out2, shell=True)

    # This removes the sequence colulmn from the tab fasta file created above
    command3 = "cut -f1,1 {0}".format("removed_SigPep_tab.fasta.txt")
    file_out3 = open(path + "listaftercut.txt", "w")
    file_out4 = open(path + "goodlistSigP.txt", "w")
    listGood = subprocess.call(command3.split(), stdout=file_out3, shell=True)
    file_out3.close()

    openfile = open(path + "listaftercut.txt", 'r')
    for line in openfile:
            goodname = line.partition(' ')[0] + '\n'
            file_out4.write(goodname)
    file_out2.close()
    file_out4.close()
    openfile.close()

def sigpFasta(verbose=False):
    """
    This function creates a fasta file containing the complete sequences with signal peptides
    """
    command4 = "faSomeRecords {0} {1} {2}".format("singleline.fasta",
                                                  "goodlistSigP.txt",
                                                  "signalP_pass.fasta")
    print_verobse("Retreving SignalP fasta", v_flag=verbose)
    fastaRUN = subprocess.call(command4.split(), shell=True)


def tmhmm(verbose=False):
    """
    This function runs tmhmm on the sequences with signal peptides inorder to check for transmembramne domains.

    NOTE this uses the mature sequences only so any TM regions in the signal peptide will be avoided/ignored


        NOTE in the file tmhmmformat.pl set the follwoing lines as shown:

            $opt_html = 0;       # Make HTML output
            $opt_short = 1;      # Make short output format
            $opt_plot = 0;       # Make plots
            $opt_v1 = 0;         # Use old model (version 1)
            $opt_signal = 10;    # Cut-off used for producing a warning of possible
                                 # Signal sequence
            $opt_Nterm = 0;      # Number of bases to consider for signal peptides
                                 # in the Nterm

    """
    command = "tmhmm {0}".format("removed_SigPep.fasta")

    file_out = open(path + "tmhmmOUT.txt", "w")
    print_verbose("Running tmhmm on mature signalp sequences only", v_flag=verbose)
    tmhmmRUN = subprocess.call(command.split(), stdout=file_out, shell=True)
    file_out.close()
    print_verbose("tmhmm complete", v_flag=verbose)
    print_verbose("Identifying sequences without tm regions.", v_flag=verbose)

    # This section of code parses the output from tmhmm and collects fastas with no TM regions
    openfile = open(path + "tmhmmOUT.txt", "r")
    file_out2 = open(path + "tmhmmGoodlist.txt", "a")
    for line in openfile:
            if "\tPredHel=0\t" in line:
                    goodname = line.partition('\t')[0] + '\n'
                    file_out2.write(goodname)
    openfile.close()
    file_out2.close()

def targetp(verbose=False):
    """
    This function uses targetp to verify the destination of the signal peptide
    NOTE for plant networks use -P over -N
    """

    command = "targetp -N {0}".format("signalP_pass.fasta")
    file_out = open(path + "targetpOUT.txt", "w")
    print_verbose("Running TargetP on SignalP pass seqeunces only", v_flag=verbose)
    targetpRUN = subprocess.check_call(command.split(), stdout=file_out, shell=True)
    print_verbose("TargetP complete", v_flag=verbose)
    file_out.close()
    print_verbose("Identifying sequences that are secreated.", v_flag=verbose)

    # Removes the leader info from the out file
    lines = open(path + 'targetpOUT.txt').readlines()
    open(path + 'targetpOUT_parse.txt', 'w').writelines(lines[8:-2])

    # Puts fastas identified as secreted "S" into goodlist
    openfile = open(path + "targetpOUT_parse.txt", "r")
    file_out2 = open(path + "targetpGoodlist.txt", "a")
    for line in openfile:
            if "S" in line:
                    goodname = line.partition(' ')[0] + '\n'
                    file_out2.write(goodname)
    openfile.close()
    file_out2.close()

def wolfpsort(verbose=False):
    """
    runs wolfPsort with the fungi setting. change as necessary for you usage.
    """
    command = "runWolfPsortSummary fungi < {0}".format("singleline.fasta")
    file_out = open(path + "wolfPsortOUT.txt", "w")
    file_out2 = open(path + "wolfPsortErrorLog.txt", "w")
    print_verbose("Running WoLFPSORT", v_flag=verbose)
    wolfRUN = subprocess.check_call([command], stdout = file_out, stderr=file_out2, shell=True)
    file_out.close()
    file_out2.close()
    print_verbose("WoLFPSORT complete", v_flag=verbose)

    # Removes header from output file
    lines = open(path + 'wolfPsortOUT.txt').readlines()
    open(path + 'wolfPsortOUT_parse.txt', 'w').writelines(lines[1:])

    file_out2 = open(path + "wolfPsortGoodlist.txt", "a")
    # Places fastas with extracellualr location into good list
    searchValue = "extr"
    f = open(path + "wolfPsortOUT_parse.txt", "r+b")
    for line in f:
            if line.split()[1] == searchValue:
                    goodname = line.partition(' ')[0] + '\n'
                    file_out2.write(goodname)
    f.close()


def secretome(signalp_acc_with_sigpep,
              tmhmm_acc_with_no_tm,
              targetp_secreted_acc,
              wolfpsort_extracellular_acc,
              output_file):

    sig_peptides  = set(line.strip() for line in open(signalp_acc_with_sigpep))
    no_tm_domains = set(line.strip() for line in open(tmhmm_acc_with_no_tm))
    secreted      = set(line.strip() for line in open(targetp_secreted_acc))
    extracellular = set(line.strip() for line in open(wolfpsort_extracellular_acc))

    with open(output_file, 'w') as out_fh:
        for line in file1 & file2 & file3 & file4:
                print(line)
                if line:
                        print(line)
                        out_fh.write(line + '\n')

def generate_ouput(input_file, output_file,
                   rename_mappings, secretome_list,
                   verbose=False):

    in_handle = input_file
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


def main(argv):
    """
    Main execution of the program in the proper order
    input: argv from arg parser
    """

    bin_path = os.path.abspath('dependencies/bin')
    check_dependencies(bin_path)

    input_file = argv.input_file


    if argv.run_name is False:
        basename = os.path.basename(input_file)
        run_name = os.path.splitext(basename)[0]


    pwd = os.getcwd()
    tmp_dir = os.path.join(pwd, 'intermediate_outputs_'+run_name)
    try:
        os.makedirs(tmp_dir)
    except OSError:
        if os.path.exists(tmp_dir):
            warnings.warn('intermediate output dir: {0} exists '
                          'overwriting contents'.format(tmp_dir))
        else:
            raise OSError('Error creating intermediate '
                          'output dir: {0}'.format(tmp_dir))

    formatted_fasta = os.path.join(tmp_dir, "formatted_input.fasta")

    renaming_mappings = format_fasta(input_file, formatted_fasta, argv)

    #signalp()
    #sigpFasta()
    #tmhmm()
    #targetp()
    #wolfpsort()
    #secretome()
    #generate_ouput(renaming_mappings)

if __name__=='__main__':

    parser = get_parser()
    argv = parser.parse_args()

    main(argv)

