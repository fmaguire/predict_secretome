#!/usr/bin/env python3

import predict_secretome.utils as utils
import os
import warnings
import shutil

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



    rename_mappings, formatted_fasta_fp = utils.format_fasta(input_file,
                                                             tmp_dir,
                                                             verbose=argv.verbose)


    mature_seqs_fp, \
    accessions_with_sig_pep, \
    full_sequences_with_sigpep_fp = utils.signalp(formatted_fasta_fp,
                                                  tmp_dir,
                                                  path=bin_path,
                                                  verbose=argv.verbose)

    accesions_no_tm_in_mature_seq = utils.tmhmm(mature_seqs_fp,
                                               tmp_dir,
                                               path=bin_path,
                                               verbose=argv.verbose)

    secreted_accessions = utils.targetp(full_sequences_with_sigpep_fp,
                                        tmp_dir,
                                        plant=False,
                                        path=bin_path,
                                        verbose=argv.verbose)


    extracellular_accessions = utils.wolfpsort(formatted_fasta_fp,
                                               tmp_dir,
                                               path=bin_path,
                                               verbose=argv.verbose)


    if argv.permissive:
        conservative_flag = False
    else:
        conservative_flag = True

    secretome_accessions = utils.secretome(accessions_with_sig_pep,
                                           accesions_no_tm_in_mature_seq,
                                           secreted_accessions,
                                           extracellular_accessions,
                                           tmp_dir,
                                           conservative=conservative_flag,
                                           verbose=argv.verbose)

    predicted_secretome = utils.generate_output(formatted_fasta_fp,
                                                secretome_accessions,
                                                rename_mappings,
                                                run_name,
                                                conservative=conservative_flag,
                                                verbose=argv.verbose)

    if not argv.nocleanup:
        shutil.rmtree(tmp_dir)

if __name__=='__main__':

    parser = utils.get_parser()
    argv = parser.parse_args()

    main(argv)
