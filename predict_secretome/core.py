import os
import warnings
import sys
import shutil
from predict_secretome import utils

class predictSecretome(object):
    """
    Runner class for predicting the secretome
    from a protein fasta input file
    """
    def __init__(self, argv):
        """
        Initialise the namespace mainly
        by:
        - converting user input to a self.settings dicts
        - initialising empty dict for predictions from each tools - self.outputs
        - initialise dict containing all the relevant file paths - self.files
        - generating directory structure for run
        """
        self.dep_path = os.path.join('predict_secretome', 'dependencies',
                                     'bin')

        self.settings = {'check_run': argv.check,
                         'verbose': argv.verbose,
                         'force_overwrite': argv.force,
                         'permissive': argv.permissive}


        # parse the typing argument
        if argv.plant:
            euk_type = 'plant'
        elif argv.animal:
            euk_type = 'animal'
        elif argv.fungi:
            euk_type = 'fungi'
        self.settings.update({'euk_type': argv.euk_type})


        # run check dependencies util
        utils.check_dependencies(self.dep_path,
                                 self.settings['check_run'])

        self.files = {'input_fasta': argv.input_file}

        if not argv.run_name:
            basename = os.path.basename(self.files['input_fasta'])
            argv.run_name = os.path.splitext(basename)[0]

        self.run_name = argv.run_name

        self.settings.update({'output_dir': os.path.join(os.getcwd(),
            '{}_secretome_prediction_output'.format(self.run_name))})

        self.settings.update({'work_dir': os.path.join(self.settings['output_dir'],
                                                           "intermediate_outputs")})

        self.outputs = {}

        self._initialise_directory_structure()



    def _initialise_directory_structure(self):
        """
        Initialise the output_structure checking
        if it already exists and whether to overwrite it.
        """
        path_exists = os.path.exists(self.settings['output_dir'])

        if path_exists and not self.settings['force_overwrite']:
            sys.exit('Output directory ({}) already exists.'
                ' Delete or re-run pipeline with --force option to overwrite'.format(self.settings['output_dir']))
        elif path_exists and self.settings['force_overwrite']:
            shutil.rmtree(self.settings['output_dir'])

        # makedirs recursively makes full path
        # i.e. outputdir and outputdir/workdir
        os.makedirs(self.settings['work_dir'])


    def _reformat_input(self):
                """
        Format input fasta and write to output:
        Specifically replace accessions with hash of name to guarantee under 20 chars
        input: input_file - filename
               output_file - filename
               path - path to dependencies bins
        output: rename_mappings - a dict of {new_acc: original_acc}
        """

        print("\n##Formatting input fasta: {0}##".format(self.files['input_fasta']))


        formatted_fasta_fp = os.path.join(self.settings['work_dir'],
                                          "formatted_input.fasta")

        # as some dependencies (tmhmm) don't support accessions over 20 chars
        # replace accessions with a 19-char truncated hash of original accession
        # keeping track in a dict where key is the hash and val is original accession
        raw_fasta_in_fh = open(self.files['input_fasta'], "r")
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

            rename_mappings.update({truncated_md5: record.description})

            record.id = truncated_md5
            record.name = ''
            record.description = ''

            fasta_out.write_record(record)

        fasta_out.write_footer()

        raw_fasta_in_fh.close()
        formatted_out_fh.close()

        print("Input fasta formatted: {0}".format(formatted_fasta_fp))

        # update namespace on completion
        self.files.update({'reformatted_input': formatted_fasta_fp})
        self.rename_mappings = rename_mappings



    def _signalp(self):
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

        # Dependency
        if 'reformatted_input' not in self.files.keys():
            self._reformat_input()


        # mature seqs are those with signal peptides removed
        mature_seqs_fp = os.path.join(self.settings['work_dir'],
                                      'signalp_mature_seqs.fasta')

        sigp = os.path.join(self.dep_path, 'signalp')
        sigp_cmd = "{0} -t euk -f short -m {1} {2}".format(sigp,
                                                           mature_seqs_fp,
                                                           self.files['reformatted_input'])


        # we only care about the mature sequences output as they are the sequences
        # predicted as having signal peptides (which have then been subsequently
        # trimmed)

        print('##Detecing signal peptides##')
        sigp_proc = subprocess.Popen(sigp_cmd.split(),
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
        _, stderr = sigp_proc.communicate()

        print("Search Complete")

        if '# No sequences predicted with a signal peptide' in stderr.decode('ascii'):
            print("No signal peptides detected no secretome will be predicted")
            # if there are no signal peptides can still calculate and return tm if specified
            sys.exit()


        print("Collecting accessions with signal peptides")
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
        full_sequences_with_sigpep_fp = os.path.join(self.settings['work_dir'],
                                                     'signalp_full_seqs_with_sigpep.fasta')

        print("Assembling full sequences with signal peptides")
        with open(full_sequences_with_sigpep_fp, 'w') as sigpep_seqs_fh:
            fasta_out = SeqIO.FastaIO.FastaWriter(sigpep_seqs_fh, wrap=None)
            fasta_out.write_header()
            for record in SeqIO.parse(input_file, 'fasta'):
                if record.description in accessions_with_sig_pep:
                    fasta_out.write_record(record)

            fasta_out.write_footer()

        self.files.update({'mature_sequences': mature_seqs_fp,
                           'sequences_with_sig': full_sequences_with_sigpep_fp})

        self.outputs.update({'acc_with_sig': accessions_with_sig_pep})


    def _tmhmm(self):
        """
        Run tmhmm on the mature sequences from signalp to ensure they
        don't have transmembrane domains outwith of their signal peptide
        input: mature_seqs_fp - fasta file without signal peptides
        output: acc_without_tm_in_mature_seq - list of accessions
                                               without TM domains
                                               in their mature
                                               sequences
        """
        if "mature_seqs" not in self.files.keys():
            self._signalp()

        tmhmm = os.path.join(self.dep_path, 'tmhmm')
        tmhmm_cmd = "{0} {1}".format(tmhmm,
                                     self.files['mature_sequences'])


        print("\n##Search for TM domains in mature seqs##")
        tmhmm_output = subprocess.check_output(tmhmm_cmd.split())
        tmhmm_output = tmhmm_output.decode('ascii').split('\n')
        print("Search complete")

        print("Parsing results")
        # Parse tmhmm raw output and write acc without tm domains
        # in non signal peptide sequence to output_file
        acc_without_tm_in_mature_seq = [line.split('\t')[0] \
                                          for line in tmhmm_output \
                                            if "\tPredHel=0\t" in line]

        self.outputs.update({'non_tm_mature_accs': acc_without_tm_in_mature_seq})


    def _targetp(self):
       """
       Runs targetp on full seqs of those identified as having
       signal peptides by signalp to identify those designated
       as targeted 'secreted'
       input: mature_seqs_fp - fasta file without signal peptides
              plant - boolean to use plant settings or not for targetp
        /utput: secreted_accessions - list of accessions with signalpeps
                                      predicted to be secreted

        """
        if 'sequences_with_sig' not in self.files.keys():
            self._signalp()

        if self.settings['euk_type'] == 'plant':
            targetp_flag = "-P"
        else:
            targetp_flag = "-N"

        # the fortran 'How' ANN classifier dependency bundled with targetp
        # can't handle paths longer than 80 chars so as a hack fix I'm creating
        # a symlink to /home/user and removing it after execution
        username = os.getlogin()
        home_targetp = '/home/{0}/targetp-1.1'.format(username)
        shutil.copytree(os.path.join(self.dep_path, 'targetp-1.1'), home_targetp)

        targetp_path = os.path.join(self.dep_path, 'targetp')
        targetp_cmd = "{0} {1} {2}".format(targetp_path,
                                           targetp_flag,
                                           full_sequences_with_sigpep_fp)

        print("##Identifying sequences with 'secreted' sigpeps##")
        targetp_output = subprocess.check_output(targetp_cmd.split())
        targetp_output = targetp_output.decode('ascii').split('\n')
        print("Search complete")
        shutil.rmtree(home_targetp)

        print("Parsing results")
        # remove header and tail cruft in targetp output
        # and get those that have S as top predicted target
        secreted_accessions = [line.split(' ')[0] \
                               for line in targetp_output[8:-3] if "S" in line]

        self.outputs.update({'secreted_acc', secreted_accessions})

    def _wolfpsort(self)
        """
        Run wolfPsort on all formatted sequences using fungi setting
        to get a predicted list of 'extracellular' accessions
        input: formatted_fasta_fp - formatted fasta file of all input seqs
        output: extracellular_accessions - list of accessions with signalpeps
                                            predicted to be extracellular

        """

        if 'reformatted_input' not in self.files.keys():
            self._reformat_input()

        wolfpsort_path = os.path.join(path, 'runWolfPsortSummary')

        wolfpsort_cmd = "{0} {1} < {2}".format(wolfpsort_path,
                                               self.settings['euk_type'],
                                               self.files['reformatted_input'])

        print("\n##Identifying sequences belonging to 'extracellular' compartment##")

        wolfpsort_output = subprocess.check_output(wolfpsort_cmd, shell=True)
        wolfpsort_output = wolfpsort_output.decode('ascii').split('\n')
        print("Search complete")


        print("Parsing results")
        # Removes header from output
        extracellular_accessions = [line.split(' ')[0] \
                                      for line in wolfpsort_output[1:-1] \
                                           if "extr" in line.split()[1]]

        self.outputs.update({'extracellular_acc': extracellular_accessions})


    def _get_secretome_accessions(self):
        """
        Get secretome accessions from the predictions
        """
        # if predictions not done, then do them
        print("\n##Combining predictions##")

        sig_peptides = set(self.outputs['acc_with_sig'])
        no_tm_domains = set(self.outputs['non_tm_mature_accs'])
        secreted = set(self.outputs['secreted_accessions'])
        extracellular = set(self.outputs['extracellular_acc'])


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

        if not permissive:
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
            print("No secreted proteins found using {0} setting".format(out_flag))

        self.output.update({'predicted_secretome_acc': predicted_acc_list})



    def _output_secretome(self):
        """
        Using the secretome accessions retrieve sequences from
        the formatted renamed input fasta and reconstitute their
        names.

        Output file to output fasta.
        """

        # Dependency
        if self.secretome_accessions is None:
            self.get_secretome_accessions()


        print("\n##Writing predicted secretome fasta file##")


            #if self.settings['permissive']:
            #    self.files.update({'secretome_fasta':



        formatted_input_fh = open(self.files['formatted_fasta'], 'r')


        secretome_output = "{0}_{1}_predicted_secretome.fasta".format(self.run_name, out_flag)
        out_handle = open(output, 'w')



        formatted_input_fh

        fasta_out = SeqIO.FastaIO.FastaWriter(out_handle, wrap=None)
        fasta_out.write_header()

        print("Retreving Secretome fasta")
        for record in SeqIO.parse(in_handle, 'fasta'):
            if record.description in self.secretome_accessions:
                record.id = rename_mappings[record.id]
                record.name = ''
                record.description = ''
                fasta_out.write_record(record)

        in_handle.close()
        fasta_out.write_footer()
        out_handle.close()
        print("Secretome identification Complete")




    def run_predict(self):
        """
        Main runner method for the class
        -
        """
        self._reformat_input()
        self._signalp()
        self._tmhmm()
        self._targetp()
        self._wolfpsort()
        self._
        self._output_secretome()
        return
