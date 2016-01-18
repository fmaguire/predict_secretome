import os


class predictSecretome(object):
    """
    Runner class for predicting the secretome
    from a protein fasta input file
    """
    def __init__(self, argv):
        """
        Initialise the namespace mainly
        by converting
        """

        self.verbose = argv.verbose
        self.check_run = argv.check

        self.verbose = argv.check

        utils.check_dependencies(check_run=self.check,
                                 verbose=self.verbose)


        self.input_file = argv.input_file


        if not argv.run_name:
            basename = os.path.basename(self.input_file)
            self.run_name = os.path.splitext(basename)[0]
        else:
            self.run_name = argv.run_name

        self.tmp_dir = os.path.join(os.getcwd(),
                                    'tmp_output_{}'.format(self.run_name))

        self.initialise_directory_structure()


    def initialise_directory_structure(self):
        """
        Initialise the output_structure
        """
        try:
            os.makedirs(self.tmp_dir)
        except OSError:
            if os.path.exists(self.tmp_dir) and not self.force:
                warnings.warn('Temporary Output Dir: {} overwritten')


    def reformat_input(self):
        """
        Verify and reformat input fasta
        """
        self.renaming = rename_mappings, formatted_fasta_fp = utils.format_fasta(self.input



class prediction():

    #mature_seqs_fp, \
    #accessions_with_sig_pep, \
    #full_sequences_with_sigpep_fp = utils.signalp(formatted_fasta_fp,
    #                                              tmp_dir,
    #                                              rename_mappings,
    #                                              run_name,
    #                                              trans=argv.trans,
    #                                              verbose=argv.verbose)

    #if argv.trans:
    #    utils.detect_and_output_transporters(formatted_fasta_fp,
    #                                         rename_mappings,
    #                                         tmp_dir,
    #                                         run_name,
    #                                         mature_seqs=mature_seqs_fp,
    #                                         verbose=argv.verbose,
    #                                         tm_threshold=argv.trans)


    #accesions_no_tm_in_mature_seq = utils.tmhmm(mature_seqs_fp,
    #                                            verbose=argv.verbose)

    #secreted_accessions = utils.targetp(full_sequences_with_sigpep_fp,
    #                                    plant=False,
    #                                    verbose=argv.verbose)


    #extracellular_accessions = utils.wolfpsort(formatted_fasta_fp,
    #                                           verbose=argv.verbose)

    #if argv.permissive:
    #    conservative_flag = False
    #else:
    #    conservative_flag = True



    def mature_peptides(self):
        """
        All sequences without predicted
        signal peptides and the mature (i.e. cleaved)
        sequence of those with them
        """
        pass



    def non_membrane_bound(self):
        """
        Accessions without TM domains in their
        mature peptides
        """
        pass

        pass



#def main(argv):
#    """
#    Main execution of the program in the proper order
#    input: argv from arg parser
#    """
#
#
#
#input_file = argv.input_file
#
#    if argv.run_name is False:
#        basename = os.path.basename(input_file)
#        run_name = os.path.splitext(basename)[0]
#    else:
#        run_name = argv.run_name
#
#    tmp_dir = os.path.join(os.getcwd(), 'intermediate_outputs_'+run_name)
#
#    try:
#        os.makedirs(tmp_dir)
#    except OSError:
#        if os.path.exists(tmp_dir):
#            warnings.warn('\n\nintermediate output dir: {0} exists '
#                          'overwriting contents\n'.format(tmp_dir))
#        else:
#            raise OSError('Error creating intermediate '
#                          'output dir: {0}'.format(tmp_dir))
#
#
#
#    rename_mappings, formatted_fasta_fp = utils.format_fasta(input_file,
#                                                             tmp_dir,
#                                                             verbose=argv.verbose)
#
#    mature_seqs_fp, \
#    accessions_with_sig_pep, \
#    full_sequences_with_sigpep_fp = utils.signalp(formatted_fasta_fp,
#                                                  tmp_dir,
#                                                  rename_mappings,
#                                                  run_name,
#                                                  trans=argv.trans,
#                                                  verbose=argv.verbose)
#
#    if argv.trans:
#        utils.detect_and_output_transporters(formatted_fasta_fp,
#                                             rename_mappings,
#                                             tmp_dir,
#                                             run_name,
#                                             mature_seqs=mature_seqs_fp,
#                                             verbose=argv.verbose,
#                                             tm_threshold=argv.trans)
#
#
#    accesions_no_tm_in_mature_seq = utils.tmhmm(mature_seqs_fp,
#                                                verbose=argv.verbose)
#
#    secreted_accessions = utils.targetp(full_sequences_with_sigpep_fp,
#                                        plant=False,
#                                        verbose=argv.verbose)
#
#
#    extracellular_accessions = utils.wolfpsort(formatted_fasta_fp,
#                                               verbose=argv.verbose)
#
#    if argv.permissive:
#        conservative_flag = False
#    else:
#        conservative_flag = True
#
#    secretome_accessions = utils.secretome(accessions_with_sig_pep,
#                                           accesions_no_tm_in_mature_seq,
#                                           secreted_accessions,
#                                           extracellular_accessions,
#                                           conservative=conservative_flag,
#                                           verbose=argv.verbose)
#
#    utils.generate_output(formatted_fasta_fp,
#                          secretome_accessions,
#                          rename_mappings,
#                          run_name,
#                          conservative=conservative_flag,
#                          verbose=argv.verbose)
#
#
#    if not argv.nocleanup:
#        shutil.rmtree(tmp_dir)
#
