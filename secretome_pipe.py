#!/usr/bin/env python3

import predicted_secretome.utils as utils
import os
import warnings

"""
This is a script desiged to take a multiple protien fasta as input and provide
the predicted secretome as output. The program requires that signalp, tmhmm,
targetp, chlorop, faSomeRecords, wolfpsort, and fasta_formatter be in the
users PATH. See those programs documentation for installation instructions.
Settings are for fungi change inputs for programs as necessary (targetp,
and wolfpsort)
"""


def main(argv):
    """
    Main execution of the program in the proper order
    input: argv from arg parser
    """

    bin_path = os.path.abspath('dependencies/bin')
    utils.check_dependencies(path=bin_path,
                             check_run=argv.check,
                             verbose=argv.verbose)

    input_file = argv.input_file

    if argv.run_name is False:
        basename = os.path.basename(input_file)
        run_name = os.path.splitext(basename)[0]

    tmp_dir = os.path.join(os.getcwd(), 'intermediate_outputs_'+run_name)

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

    renaming_mappings = utils.format_fasta(input_file,
                                           formatted_fasta,
                                           verbose=argv.verbose)


    seqs_with_sigpep_remove = utils.signalp(formatted_fasta,
                                            seqs_with_sigpep_removed
                                            tmp_dir,
                                            path=bin_path
                                            verbose=argv.verbose)

    acc_with_no_tm_domains = utils.tmhmm(seqs_with_sigpep_remove,
                                         tmp_dir,
                                         path=bin_path
                                         verbose=argv.verbose)

    targetp_secreted_acc = utils.targetp(
                                     tmp_dir,
                                   path=bin_path
                                   verbose=argv.verbose)


    wolfpsort_extracellular_acc = utils.wolfpsort(formatted_fasta,
                                            tmp_dir,
                                            path=bin_path
                                            verbose=argv.verbose)


    secretome_acc = utils.secretome(
                              acc_with_no_tm_domains,
                              targetp_secreted_acc,
                              wolfpsort_extracelluar_acc,
                              tmp_dir,
                              verbose=argv.verbose)

    predicted_secretome = utils.generate_ouput(formatted_fasta,
                                         secretome_acc,
                                         renaming_mappings,
                                         verbose=argv.verbose)

    if not argv.nocleaup:
        os.rmdir(tmp_dir)

if __name__=='__main__':

    parser = utils.get_parser()
    argv = parser.parse_args()

    utils.main(argv)

