#!/usr/bin/env python3

import predict_secretome.utils as utils
import os
import warnings
import shutil

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

    bin_path = os.path.abspath(os.path.join('dependencies','bin'))

    utils.check_dependencies(path=bin_path,
                             check_run=argv.check,
                             verbose=argv.verbose)

    input_file = argv.input_file

    if argv.run_name is False:
        basename = os.path.basename(input_file)
        run_name = os.path.splitext(basename)[0]
    else:
        run_name = argv.run_name 

    tmp_dir = os.path.join(os.getcwd(), 'intermediate_outputs_'+run_name)

    try:
        os.makedirs(tmp_dir)
    except OSError:
        if os.path.exists(tmp_dir):
            warnings.warn('intermediate output dir: {0} exists overwriting contents'.format(tmp_dir))
        else:
            raise OSError('Error creating intermediate '
                          'output dir: {0}'.format(tmp_dir))



    renaming_mappings, formatted_fasta = utils.format_fasta(input_file,
                                                            tmp_dir,
                                                            verbose=argv.verbose)


    seqs_sigpep_removed, \
    acc_with_sigpeps, \
    full_seqs_with_sigpeps = utils.signalp(formatted_fasta,
                                           tmp_dir,
                                           path=bin_path,
                                           verbose=argv.verbose)

    acc_with_no_tm_domains = utils.tmhmm(seqs_sigpep_removed,
                                         tmp_dir,
                                         path=bin_path,
                                         verbose=argv.verbose)

    targetp_secreted_acc = utils.targetp(full_seqs_with_sigpeps,
                                         tmp_dir,
                                         path=bin_path,
                                         verbose=argv.verbose)


    wolfpsort_extracellular_acc = utils.wolfpsort(formatted_fasta,
                                                  tmp_dir,
                                                  path=bin_path,
                                                  verbose=argv.verbose)

    if argv.permissive:
        conservative_flag = False
    else:
        conservative_flag = True

    secretome_acc = utils.secretome(acc_with_sigpeps,
                                    acc_with_no_tm_domains,
                                    targetp_secreted_acc,
                                    wolfpsort_extracellular_acc,
                                    tmp_dir,
                                    conservative=conservative_flag,
                                    verbose=argv.verbose)

    predicted_secretome = utils.generate_output(formatted_fasta,
                                               secretome_acc,
                                               renaming_mappings,
                                               run_name,
                                               conservative=conservative_flag,
                                               verbose=argv.verbose)

    if not argv.nocleanup:
        shutil.rmtree(tmp_dir)

if __name__=='__main__':

    parser = utils.get_parser()
    argv = parser.parse_args()

    main(argv)
